"""Testable input, ranking and reproducible-export helpers for the CRISPR app."""

from __future__ import annotations

import csv
import io
import json
import re
import zipfile
from collections import Counter
from collections.abc import Iterable, Mapping
from datetime import date, datetime

import numpy as np
import pandas as pd


GENE_COLUMN_PATTERN = re.compile(r"^([A-Za-z0-9][A-Za-z0-9_.-]*)\s*\(([0-9]+)\)$")
SYMBOL_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
RANKING_COLUMNS = [
    "gene_raw", "gene", "gene_upper", "mean_score", "rank", "percentile",
    "n_valid", "n_missing", "missing_fraction",
]
DIAGNOSTIC_COLUMNS = [
    "column", "position", "gene", "status", "reason", "n_valid", "n_missing",
    "n_invalid", "n_infinite",
]
LINEAGE_COLUMNS = ["gene", "lineage", "crispr_score", "cell_line", "cell_line_name"]
SUMMARY_COLUMNS = ["gene", "lineage", "n_valid", "median", "mean", "std", "min", "max"]
CELL_LINE_ID_NAMES = (
    "modelid", "depmap_id", "depmapid", "cell_line", "cell_line_name",
    "celllinename", "strippedcelllinename",
)
CELL_LINE_NAME_NAMES = (
    "cell_line_display_name", "celllinename", "cell_line_name",
    "strippedcelllinename", "stripped_cell_line_name",
)
UPLOAD_HEADERS = {
    "gene", "genes", "symbol", "symbols", "gene_symbol", "gene_symbols",
    "gene symbol", "gene symbols", "genename", "gene_name", "gene name",
    "基因", "基因名", "基因名称", "基因符号",
}


def read_score_csv(source) -> pd.DataFrame:
    """Read a score matrix without pandas silently renaming duplicate headers.

    Accept a filesystem path, UTF-8 bytes, or a file-like object. Duplicate
    headers remain duplicates so ``analyze_dataset`` can exclude every
    ambiguous gene column instead of silently retaining the first occurrence.
    """
    if isinstance(source, bytes):
        stream = io.TextIOWrapper(io.BytesIO(source), encoding="utf-8-sig", newline="")
    elif hasattr(source, "read"):
        position = source.tell()
        source.seek(0)
        content = source.read()
        source.seek(position)
        stream = (io.TextIOWrapper(io.BytesIO(content), encoding="utf-8-sig", newline="")
                  if isinstance(content, bytes) else io.StringIO(content.lstrip("\ufeff")))
    else:
        stream = open(source, encoding="utf-8-sig", newline="")
    try:
        header = next((row for row in csv.reader(stream) if any(field.strip() for field in row)), None)
        if header is None:
            raise pd.errors.EmptyDataError("No columns to parse from file")
        normalized = [column.strip().casefold() for column in header]
        metadata_names = set(CELL_LINE_ID_NAMES + CELL_LINE_NAME_NAMES)
        duplicates = {
            name for name, count in Counter(normalized).items()
            if count > 1 and (name in metadata_names or "lineage" in name)
        }
        if duplicates:
            labels = list(dict.fromkeys(column for column in header if column.strip().casefold() in duplicates))
            raise ValueError(
                "Duplicate metadata columns are not supported: " + ", ".join(labels)
                + ". Rename or remove the repeated lineage or cell-line identifier/name columns."
            )
        stream.seek(0)
        frame = pd.read_csv(stream, low_memory=False)
        if not isinstance(frame.index, pd.RangeIndex):
            raise ValueError("CSV rows contain more fields than the header; check delimiters and column counts")
        if len(header) != len(frame.columns):
            raise ValueError("CSV header count does not match the parsed score matrix")
        frame.columns = header
        return frame
    finally:
        stream.close()


def parse_gene_text(text: str) -> dict:
    """Split common pasted-list delimiters and deduplicate case-insensitively."""
    tokens = [token.upper() for token in re.split(r"[\s,，;；]+", str(text or "").strip()) if token]
    genes, duplicates, seen = [], [], set()
    for token in tokens:
        if token in seen:
            duplicates.append(token)
        else:
            genes.append(token)
            seen.add(token)
    return {"genes": genes, "input_count": len(tokens), "duplicates": duplicates}


def parse_gene_upload(content: bytes, filename: str) -> dict:
    """Read text lists or the first column of CSV/TSV, with an optional header.

    Only recognized gene-list headings are removed; a headerless first gene is
    always preserved. CSV files may contain extra columns, which are ignored.
    """
    text = content.decode("utf-8-sig")
    if not filename.lower().endswith((".csv", ".tsv")):
        lines = text.splitlines()
        nonempty = next((i for i, line in enumerate(lines) if line.strip()), None)
        if nonempty is not None and lines[nonempty].strip().casefold() in UPLOAD_HEADERS:
            lines.pop(nonempty)
        return parse_gene_text("\n".join(lines))

    delimiter = "\t" if filename.lower().endswith(".tsv") else ","
    try:
        delimiter = csv.Sniffer().sniff(text[:8192], delimiters=",\t;，；").delimiter
    except csv.Error:
        pass
    rows = csv.reader(io.StringIO(text), delimiter=delimiter)
    first_column = [row[0].strip() for row in rows if row and row[0].strip()]
    if first_column and first_column[0].casefold() in UPLOAD_HEADERS:
        first_column.pop(0)
    return parse_gene_text("\n".join(first_column))


def _numeric_values(series: pd.Series) -> tuple[pd.Series, dict]:
    numeric = pd.to_numeric(series, errors="coerce")
    numeric_array = numeric.to_numpy(dtype=float, na_value=np.nan)
    finite = np.isfinite(numeric_array)
    if pd.api.types.is_numeric_dtype(series.dtype):
        n_invalid = 0
    else:
        empty = series.isna() | series.astype("string").str.strip().eq("").fillna(False)
        n_invalid = int((numeric.isna() & ~empty).sum())
    n_valid = int(finite.sum())
    counts = {
        "n_valid": n_valid,
        "n_missing": len(series) - n_valid,
        "n_invalid": n_invalid,
        "n_infinite": int(np.isinf(numeric_array).sum()),
    }
    return numeric.where(finite).astype(float), counts


def analyze_dataset(df: pd.DataFrame, gene_columns: Iterable[str] | None = None):
    """Rank finite per-gene means, using schema or explicit column selection.

    Default identification requires ``SYMBOL (numeric Entrez ID)``. Explicit
    selection also accepts plain gene symbols. Ambiguous symbols are excluded
    as a group, even when only one of their columns contains usable scores.
    ``n_missing`` includes blank, invalid and infinite values; the latter two
    are also counted separately in the diagnostics table.
    """
    selected = None if gene_columns is None else set(gene_columns)
    if selected is not None:
        missing = selected.difference(df.columns)
        if missing:
            raise ValueError("Selected gene columns not found: " + ", ".join(sorted(map(str, missing))))

    candidates, records = [], []
    for position, column in enumerate(df.columns):
        raw = str(column).strip()
        match = GENE_COLUMN_PATTERN.fullmatch(raw)
        explicit = selected is not None and column in selected
        symbol = match.group(1) if match else raw if explicit and SYMBOL_PATTERN.fullmatch(raw) else None
        eligible = bool(symbol) and (selected is None or explicit)
        record = {
            "column": column, "position": position, "gene": symbol or "",
            "status": "excluded", "reason": "not_selected" if selected is not None and not explicit else "non_gene_schema",
            "n_valid": None, "n_missing": None, "n_invalid": None, "n_infinite": None,
        }
        records.append(record)
        if eligible:
            candidates.append((position, column, symbol, record))

    counts = Counter(symbol.upper() for _, _, symbol, _ in candidates)
    rows = []
    for position, column, symbol, record in candidates:
        numeric, value_counts = _numeric_values(df.iloc[:, position])
        record.update(value_counts)
        if counts[symbol.upper()] > 1:
            record["reason"] = "ambiguous_gene_symbol"
        elif not value_counts["n_valid"]:
            record["reason"] = "no_finite_values"
        else:
            record.update(status="included", reason="included")
            rows.append({
                "gene_raw": column, "gene": symbol, "gene_upper": symbol.upper(),
                "mean_score": float(numeric.mean()),
                "n_valid": value_counts["n_valid"], "n_missing": value_counts["n_missing"],
                "missing_fraction": value_counts["n_missing"] / len(df) if len(df) else np.nan,
            })
    rankings = pd.DataFrame(rows)
    if rankings.empty:
        rankings = pd.DataFrame(columns=RANKING_COLUMNS)
    else:
        rankings = rankings.sort_values("mean_score", kind="stable").reset_index(drop=True)
        rankings["rank"] = np.arange(1, len(rankings) + 1)
        rankings["percentile"] = rankings["rank"] / len(rankings) * 100
        rankings = rankings[RANKING_COLUMNS]
    return rankings, pd.DataFrame(records, columns=DIAGNOSTIC_COLUMNS)


def match_genes(rankings: pd.DataFrame, genes: Iterable[str]):
    """Return canonical matched symbols and absent symbols in request order."""
    lookup = dict(zip(rankings["gene_upper"], rankings["gene"]))
    matched, not_found, seen = [], [], set()
    for gene in genes:
        key = str(gene).strip().upper()
        if not key or key in seen:
            continue
        seen.add(key)
        if key in lookup:
            matched.append(lookup[key])
        else:
            not_found.append(str(gene).strip())
    return matched, not_found


def find_lineage_column(df: pd.DataFrame):
    """Prefer exact coarse-lineage labels; never select a sub-lineage label."""
    lookup = {str(column).strip().casefold(): column for column in df.columns}
    for name in ("lineage", "oncotreelineage", "oncotree_lineage"):
        if name in lookup:
            return lookup[name]
    return next((column for column in df.columns if "lineage" in str(column).casefold() and "sub" not in str(column).casefold()), None)


def _cell_line_ids(df: pd.DataFrame) -> pd.Series:
    lookup = {str(column).strip().casefold(): column for column in df.columns}
    for name in CELL_LINE_ID_NAMES:
        if name in lookup:
            values = df[lookup[name]].astype("string")
            fallback = pd.Series(df.index.astype(str), index=df.index, dtype="string")
            return values.mask(values.isna() | values.str.strip().eq(""), fallback)
    return pd.Series(df.index.astype(str), index=df.index, dtype="string")


def _cell_line_names(df: pd.DataFrame, identifiers: pd.Series) -> pd.Series:
    lookup = {str(column).strip().casefold(): column for column in df.columns}
    for name in CELL_LINE_NAME_NAMES:
        if name in lookup:
            values = df[lookup[name]].astype("string").str.strip()
            return values.mask(values.isna() | values.eq(""), identifiers)
    return identifiers


def get_lineage_frame(df: pd.DataFrame, genes: Iterable[str], rankings: pd.DataFrame, lineage_col=None):
    """Build finite long-form plotting data without losing cell-line identity."""
    lineage_col = find_lineage_column(df) if lineage_col is None else lineage_col
    matched, _ = match_genes(rankings, genes)
    lookup = rankings.set_index("gene_upper")["gene_raw"].to_dict()
    lineage = (df[lineage_col].astype("string").str.strip().replace("", pd.NA).fillna("Unknown")
               if lineage_col is not None and lineage_col in df.columns
               else pd.Series("Unknown", index=df.index, dtype="string"))
    identifiers = _cell_line_ids(df)
    names = _cell_line_names(df, identifiers)
    frames = []
    for gene in matched:
        scores, _ = _numeric_values(df[lookup[gene.upper()]])
        frame = pd.DataFrame({
            "gene": gene, "lineage": lineage.to_numpy(),
            "crispr_score": scores.to_numpy(), "cell_line": identifiers.to_numpy(),
            "cell_line_name": names.to_numpy(),
        })
        frames.append(frame.loc[frame["crispr_score"].notna()])
    return pd.concat(frames, ignore_index=True)[LINEAGE_COLUMNS] if frames else pd.DataFrame(columns=LINEAGE_COLUMNS)


def lineage_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize finite observations for each gene and cancer lineage."""
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    valid = frame.copy()
    valid["crispr_score"], _ = _numeric_values(valid["crispr_score"])
    valid = valid.dropna(subset=["crispr_score"])
    valid["lineage"] = valid["lineage"].astype("string").str.strip().replace("", pd.NA).fillna("Unknown")
    return valid.groupby(["gene", "lineage"], sort=False, observed=True)["crispr_score"].agg(
        n_valid="count", median="median", mean="mean", std="std", min="min", max="max",
    ).reset_index()[SUMMARY_COLUMNS]


def _json_default(value):
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, (datetime, date, pd.Timestamp)):
        return value.isoformat()
    if isinstance(value, (np.ndarray, pd.Index)):
        return value.tolist()
    if value is pd.NA:
        return None
    raise TypeError(f"Unsupported metadata type: {type(value).__name__}")


def make_analysis_bundle(rankings: pd.DataFrame, selected_genes: Iterable[str], plot_data: pd.DataFrame,
                         metadata: Mapping, diagnostics=None) -> bytes:
    """Package complete rankings, selected results, plot data and provenance."""
    selected_genes = list(selected_genes)
    matched, not_found = match_genes(rankings, selected_genes)
    indexed = rankings.set_index("gene_upper", drop=False)
    selected = indexed.loc[[gene.upper() for gene in matched]].reset_index(drop=True)
    totals = (rankings["n_valid"] + rankings["n_missing"]).dropna().unique() if not rankings.empty else []
    provenance = dict(metadata)
    provenance["analysis"] = {
        "method": "Arithmetic mean of finite CRISPR scores per gene; lower scores rank first.",
        "ranking": "Ascending mean; sequential 1-based rank; ties preserve dataset column order.",
        "percentile": "100 * rank / number of ranked genes; lower indicates stronger dependency.",
        "missing_values": "Blank, nonnumeric and infinite values excluded from means and plotted scores.",
        "gene_identification": "Strict SYMBOL (numeric Entrez ID) schema or explicit user-selected gene columns; ambiguous symbols excluded.",
        "n_ranked_genes": len(rankings),
        "n_cell_lines": int(totals[0]) if len(totals) == 1 else None,
        "n_selected_genes": len(selected),
        "n_plot_rows": len(plot_data),
        "selected_genes_requested": selected_genes,
        "selected_genes_matched": matched,
        "selected_genes_unmatched": not_found,
    }
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, frame in (("selected_results.csv", selected), ("full_rankings.csv", rankings), ("plot_data.csv", plot_data)):
            archive.writestr(name, frame.to_csv(index=False).encode("utf-8-sig"))
        archive.writestr("metadata.json", json.dumps(provenance, ensure_ascii=False, indent=2, default=_json_default))
        if diagnostics is not None:
            archive.writestr("diagnostics.csv", pd.DataFrame(diagnostics).to_csv(index=False).encode("utf-8-sig"))
    return buffer.getvalue()
