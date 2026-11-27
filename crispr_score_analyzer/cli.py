"""Command-line interface for CRISPR Score Analyzer.

Provides reproducible batch gene ranking, cohort-scoped co-dependency scanning,
dataset diagnostics, vector figure export, and Streamlit UI launching without
requiring manual browser interaction.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from . import __version__
from .analysis_core import (
    analyze_dataset,
    default_builtin_matrix_path,
    file_sha256,
    find_disease_column,
    find_lineage_column,
    gene_correlation_matrix,
    get_lineage_frame,
    lineage_preference,
    load_gene_aliases,
    make_analysis_bundle,
    parse_gene_text,
    parse_gene_upload,
    read_builtin_matrix,
    read_score_csv,
    resolve_gene_list,
)
from .codependency import (
    create_codependency_scan_plot,
    genomewide_codependency,
)
from .figure_export import render_figure_bytes
from .plotting import create_correlation_heatmap, create_rank_plot


def _load_matrix(input_path: Path | None) -> tuple[pd.DataFrame, Path, str]:
    """Load a CSV or bundled NPZ score matrix and return (frame, path, sha256)."""
    if input_path is None:
        candidate = default_builtin_matrix_path()
        if not candidate.is_file():
            raise FileNotFoundError(
                f"No --input file was provided and no bundled matrix was found at {candidate}. "
                "Pass --input path/to/crispr_scores.csv (or .npz)."
            )
        input_path = candidate
    if not input_path.is_file():
        raise FileNotFoundError(f"Input score matrix not found: {input_path}")
    if input_path.suffix.lower() == ".npz":
        df = read_builtin_matrix(input_path)
    else:
        df = read_score_csv(input_path)
    return df, input_path, file_sha256(input_path)


def _apply_scope(
    df: pd.DataFrame,
    cohorts: list[str] | None = None,
    diseases: list[str] | None = None,
) -> tuple[pd.DataFrame, str | None, str | None]:
    lineage_col = find_lineage_column(df)
    disease_col = find_disease_column(df)
    scoped = df
    if cohorts:
        if not lineage_col:
            raise ValueError("Cannot filter by --cohort because no lineage column was found in the input matrix.")
        labels = scoped[lineage_col].fillna("Unknown").astype(str).str.strip().replace("", "Unknown")
        wanted = {c.casefold() for c in cohorts}
        scoped = scoped[labels.str.casefold().isin(wanted)]
    if diseases:
        if not disease_col:
            raise ValueError("Cannot filter by --disease because no primary-disease column was found in the input matrix.")
        labels = scoped[disease_col].fillna("Unknown").astype(str).str.strip().replace("", "Unknown")
        wanted = {d.casefold() for d in diseases}
        scoped = scoped[labels.str.casefold().isin(wanted)]
    if scoped.empty:
        raise ValueError("No cell lines matched the requested --cohort / --disease filters.")
    return scoped, lineage_col, disease_col


def _read_genes_argument(genes_text: str | None, genes_file: Path | None) -> dict:
    if genes_file is not None:
        if not genes_file.is_file():
            raise FileNotFoundError(f"Gene list file not found: {genes_file}")
        return parse_gene_upload(genes_file.read_bytes(), genes_file.name)
    return parse_gene_text(genes_text or "")


def _cmd_info(args: argparse.Namespace) -> int:
    aliases = load_gene_aliases()
    payload: dict = {
        "package": "crispr-score-analyzer",
        "version": __version__,
        "hgnc_aliases_loaded": len(aliases),
    }
    target = args.input
    if target is None:
        default_path = default_builtin_matrix_path()
        target = default_path if default_path.is_file() else None
    if target is not None:
        df, path, digest = _load_matrix(target)
        gene_cols = [c.strip() for c in args.gene_columns.split(",") if c.strip()] if args.gene_columns else None
        rankings, diagnostics = analyze_dataset(df, gene_columns=gene_cols)
        lineage_col = find_lineage_column(df)
        disease_col = find_disease_column(df)
        payload["dataset"] = {
            "path": str(path),
            "sha256": digest,
            "n_cell_lines": int(len(df)),
            "n_ranked_genes": int(len(rankings)),
            "n_excluded_columns": int((diagnostics["status"] == "excluded").sum()),
            "lineage_column": lineage_col,
            "n_lineages": int(df[lineage_col].dropna().nunique()) if lineage_col else 0,
            "disease_column": disease_col,
        }
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


def _cmd_rank(args: argparse.Namespace) -> int:
    df_all, path, digest = _load_matrix(args.input)
    df_scope, lineage_col, _ = _apply_scope(df_all, args.cohort, args.disease)
    gene_cols = [c.strip() for c in args.gene_columns.split(",") if c.strip()] if args.gene_columns else None
    rankings, diagnostics = analyze_dataset(df_scope, gene_columns=gene_cols)
    if rankings.empty:
        raise ValueError(
            "No analyzable gene columns were found. DepMap-style 'SYMBOL (Entrez ID)' headers "
            "are detected automatically; pass --gene-columns for plain-symbol matrices."
        )
    aliases = {} if args.no_aliases else load_gene_aliases()
    parsed = _read_genes_argument(args.genes, args.genes_file)
    if parsed["genes"]:
        matched, missing, resolved = resolve_gene_list(rankings, parsed["genes"], aliases)
    else:
        matched = rankings["gene"].head(args.top).tolist()
        missing, resolved = [], []

    if resolved:
        pairs = ", ".join(f"{item['input']} -> {item['gene']}" for item in resolved)
        print(f"Resolved HGNC aliases: {pairs}", file=sys.stderr)
    if missing:
        print(f"Unmatched genes ({len(missing)}): {', '.join(missing)}", file=sys.stderr)
    if not matched:
        raise ValueError("None of the requested genes matched the ranked dataset.")

    selected = rankings.set_index("gene", drop=False).reindex(matched).reset_index(drop=True)
    display = selected[["gene", "mean_score", "rank", "percentile", "n_valid", "n_missing"]]
    print(display.to_string(index=False))

    if args.output_dir is not None:
        out_dir: Path = args.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        selected.drop(columns=["gene_upper"], errors="ignore").to_csv(out_dir / "selected_results.csv", index=False)
        rankings.drop(columns=["gene_upper"], errors="ignore").to_csv(out_dir / "full_rankings.csv", index=False)
        diagnostics.to_csv(out_dir / "diagnostics.csv", index=False)
        plot_data = get_lineage_frame(df_scope, matched, rankings, lineage_col)
        plot_data.to_csv(out_dir / "plot_data.csv", index=False)

        pref_rows = []
        if lineage_col:
            full_rankings, _ = analyze_dataset(df_all, gene_columns=gene_cols)
            full_lineage_df = get_lineage_frame(df_all, matched, full_rankings, lineage_col)
            pref = lineage_preference(full_lineage_df)
            if not pref.empty:
                pref.to_csv(out_dir / "lineage_preference.csv", index=False)
                pref_rows = pref.to_dict(orient="records")

        if args.plot_format:
            fig = create_rank_plot(
                rankings,
                matched,
                essential_gene=args.ref_a,
                nonessential_gene=args.ref_b,
                n_cell_lines=len(df_scope),
                width=args.width,
                height=args.height,
            )
            img_bytes = render_figure_bytes(fig, args.plot_format, width=args.width, height=args.height)
            (out_dir / f"rank_plot.{args.plot_format}").write_bytes(img_bytes)

        if args.bundle:
            metadata = {
                "software": "crispr-score-analyzer",
                "version": __version__,
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_path": str(path),
                "input_sha256": digest,
                "cohort_lineages": args.cohort or ["ALL"],
                "cohort_diseases": args.disease or [],
                "selected_genes": matched,
                "resolved_aliases": resolved,
                "unmatched_genes": missing,
                "lineage_preference": pref_rows,
            }
            zip_bytes = make_analysis_bundle(rankings, matched, plot_data, metadata, diagnostics)
            (out_dir / "analysis_bundle.zip").write_bytes(zip_bytes)

    return 0


def _cmd_codep(args: argparse.Namespace) -> int:
    df_all, _, _ = _load_matrix(args.input)
    df_scope, _, _ = _apply_scope(df_all, args.cohort, args.disease)
    gene_cols = [c.strip() for c in args.gene_columns.split(",") if c.strip()] if args.gene_columns else None
    rankings, _ = analyze_dataset(df_scope, gene_columns=gene_cols)
    if len(rankings) < 2:
        raise ValueError("At least two ranked genes are required for co-dependency analysis.")
    aliases = {} if args.no_aliases else load_gene_aliases()

    if args.query:
        parsed = parse_gene_text(args.query)
        matched, missing, resolved = resolve_gene_list(rankings, parsed["genes"], aliases)
        if resolved:
            pairs = ", ".join(f"{item['input']} -> {item['gene']}" for item in resolved)
            print(f"Resolved HGNC aliases: {pairs}", file=sys.stderr)
        if not matched:
            raise ValueError(f"Query gene not found in dataset: {', '.join(missing or [args.query])}")
        query_gene = matched[0]
        scan_df = genomewide_codependency(
            df_scope, query_gene, rankings, method=args.method, min_shared=args.min_shared,
        )
        if scan_df.empty:
            print(f"No partner genes met the minimum shared cell-line threshold for {query_gene}.")
        else:
            top_pos = scan_df.head(args.top)
            top_neg = scan_df.tail(args.top).iloc[::-1]
            print(f"Query: {query_gene} ({args.method}, {len(scan_df):,} partners scanned)")
            print("\nTop positive co-dependencies:")
            print(top_pos[["gene", "correlation", "n_shared", "mean_score"]].to_string(index=False))
            print("\nTop negative (anti-correlated) partners:")
            print(top_neg[["gene", "correlation", "n_shared", "mean_score"]].to_string(index=False))
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            scan_df.to_csv(args.output, index=False)
        if args.plot is not None and not scan_df.empty:
            fig = create_codependency_scan_plot(scan_df, query_gene, method=args.method, top_n=args.top)
            fmt = args.plot.suffix.lstrip(".").lower() or "svg"
            args.plot.parent.mkdir(parents=True, exist_ok=True)
            height = int(fig.layout.height or 560)
            args.plot.write_bytes(render_figure_bytes(fig, fmt, width=820, height=height))
        return 0

    parsed = _read_genes_argument(args.genes, args.genes_file)
    matched, missing, resolved = resolve_gene_list(rankings, parsed["genes"], aliases)
    if resolved:
        pairs = ", ".join(f"{item['input']} -> {item['gene']}" for item in resolved)
        print(f"Resolved HGNC aliases: {pairs}", file=sys.stderr)
    if missing:
        print(f"Unmatched genes ({len(missing)}): {', '.join(missing)}", file=sys.stderr)
    if len(matched) < 2:
        raise ValueError("Provide at least 2 matched genes for --genes matrix mode (or use --query for a single gene).")
    corr, counts = gene_correlation_matrix(df_scope, matched, rankings, method=args.method)
    print(corr.round(4).to_string())
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        corr.to_csv(args.output)
    if args.plot is not None:
        fig = create_correlation_heatmap(corr, counts, method=args.method)
        fmt = args.plot.suffix.lstrip(".").lower() or "svg"
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        args.plot.write_bytes(render_figure_bytes(
            fig, fmt, width=int(fig.layout.width or 680), height=int(fig.layout.height or 600),
        ))
    return 0


def _cmd_app(args: argparse.Namespace) -> int:
    app_path = Path(__file__).resolve().parents[1] / "app.py"
    if not app_path.is_file():
        raise FileNotFoundError(
            f"Streamlit entrypoint not found at {app_path}. Clone the repository to launch the web UI."
        )
    try:
        import streamlit  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Streamlit is not installed. Install UI extras with: pip install 'crispr-score-analyzer[app]'"
        ) from exc
    cmd = [sys.executable, "-m", "streamlit", "run", str(app_path)]
    if args.port:
        cmd.extend(["--server.port", str(args.port)])
    if args.headless:
        cmd.extend(["--server.headless", "true"])
    if args.extra_args:
        cmd.extend(args.extra_args)
    return subprocess.call(cmd)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="crispr-score-analyzer",
        description="Cohort-scoped CRISPR gene essentiality ranking, co-dependency analysis, and vector export.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_info = sub.add_parser("info", help="Show package metadata and optional dataset summary.")
    p_info.add_argument("-i", "--input", type=Path, default=None, help="Score matrix (.csv or .npz).")
    p_info.add_argument("--gene-columns", type=str, default=None, help="Comma-separated plain-symbol gene columns.")
    p_info.set_defaults(func=_cmd_info)

    p_rank = sub.add_parser("rank", help="Rank genes in a full or cohort-scoped CRISPR score matrix.")
    p_rank.add_argument("-i", "--input", type=Path, default=None, help="Score matrix (.csv or .npz).")
    p_rank.add_argument("-g", "--genes", type=str, default=None, help="Gene symbols (comma/space/newline separated).")
    p_rank.add_argument("--genes-file", type=Path, default=None, help="Path to a .txt/.csv/.tsv gene list.")
    p_rank.add_argument("--top", type=int, default=10, help="Top N essential genes when no --genes are specified.")
    p_rank.add_argument("--cohort", action="append", default=None, help="Filter to cancer lineage (repeatable).")
    p_rank.add_argument("--disease", action="append", default=None, help="Filter to primary disease (repeatable).")
    p_rank.add_argument("--gene-columns", type=str, default=None, help="Comma-separated plain-symbol gene columns.")
    p_rank.add_argument("--no-aliases", action="store_true", help="Disable offline HGNC alias resolution.")
    p_rank.add_argument("-o", "--output-dir", type=Path, default=None, help="Directory to write CSVs/figures/bundle.")
    p_rank.add_argument("--plot-format", choices=("svg", "pdf", "png"), default=None, help="Export rank plot format.")
    p_rank.add_argument("--width", type=int, default=640, help="Export figure width in px.")
    p_rank.add_argument("--height", type=int, default=480, help="Export figure height in px.")
    p_rank.add_argument("--ref-a", type=str, default="MYC", help="Reference gene A for rank plot.")
    p_rank.add_argument("--ref-b", type=str, default="PTEN", help="Reference gene B for rank plot.")
    p_rank.add_argument("--bundle", action="store_true", help="Also write a complete analysis_bundle.zip.")
    p_rank.set_defaults(func=_cmd_rank)

    p_codep = sub.add_parser("codep", help="Run genome-wide single-gene scan or multi-gene correlation matrix.")
    p_codep.add_argument("-i", "--input", type=Path, default=None, help="Score matrix (.csv or .npz).")
    group = p_codep.add_mutually_exclusive_group(required=True)
    group.add_argument("-q", "--query", type=str, default=None, help="Single query gene for genome-wide scan.")
    group.add_argument("-g", "--genes", type=str, default=None, help="2-25 genes for pairwise correlation matrix.")
    group.add_argument("--genes-file", type=Path, default=None, help="Path to gene list for pairwise matrix.")
    p_codep.add_argument("-m", "--method", choices=("pearson", "spearman"), default="pearson")
    p_codep.add_argument("--cohort", action="append", default=None, help="Filter to cancer lineage (repeatable).")
    p_codep.add_argument("--disease", action="append", default=None, help="Filter to primary disease (repeatable).")
    p_codep.add_argument("--gene-columns", type=str, default=None, help="Comma-separated plain-symbol gene columns.")
    p_codep.add_argument("--min-shared", type=int, default=None, help="Minimum shared finite cell lines per pair.")
    p_codep.add_argument("--top", type=int, default=10, help="Top positive/negative partners to print/plot.")
    p_codep.add_argument("--no-aliases", action="store_true", help="Disable offline HGNC alias resolution.")
    p_codep.add_argument("-o", "--output", type=Path, default=None, help="Output CSV path.")
    p_codep.add_argument("--plot", type=Path, default=None, help="Export figure path (.svg, .pdf, or .png).")
    p_codep.set_defaults(func=_cmd_codep)

    p_app = sub.add_parser("app", help="Launch the interactive Streamlit web interface.")
    p_app.add_argument("--port", type=int, default=8501, help="Port for the Streamlit server.")
    p_app.add_argument("--headless", action="store_true", help="Run Streamlit in headless mode.")
    p_app.add_argument("extra_args", nargs=argparse.REMAINDER, help="Extra arguments forwarded to streamlit run.")
    p_app.set_defaults(func=_cmd_app)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
