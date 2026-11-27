"""Testable input, ranking and reproducible-export helpers for the CRISPR app."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import zipfile
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from datetime import date, datetime
from importlib import resources
from pathlib import Path

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


_ENTREZ_SUFFIX_PATTERN = re.compile(r"(?<=[^\s,，;；])\s+\([0-9]+\)(?=$|[\s,，;；])")


def parse_gene_text(text: str) -> dict:
    """Split common pasted-list delimiters and deduplicate case-insensitively.

    Also strips DepMap-style space-separated Entrez suffixes such as
    ``MYC (4609)`` while leaving inline parenthesized aliases like ``CYSLT(1)``
    intact.
    """
    cleaned = _ENTREZ_SUFFIX_PATTERN.sub("", str(text or "").strip())
    tokens = [token.upper() for token in re.split(r"[\s,，;；]+", cleaned) if token]
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


def _numeric_block_stats(df: pd.DataFrame, positions: list[int]):
    """Finite-value counts and means for numeric columns, reduced in chunks.

    Means are float64 sums of finite values. Infinite values are excluded from
    the mean and counted separately, matching ``_numeric_values``.
    """
    width = len(positions)
    means = np.full(width, np.nan, dtype=float)
    n_valid = np.zeros(width, dtype=int)
    n_missing = np.zeros(width, dtype=int)
    n_infinite = np.zeros(width, dtype=int)
    contiguous = width > 0 and positions[-1] - positions[0] + 1 == width and all(
        positions[i] + 1 == positions[i + 1] for i in range(width - 1)
    )
    for start in range(0, width, 2048):
        stop = min(start + 2048, width)
        if contiguous:
            col_start = positions[start]
            values = df.iloc[:, col_start:col_start + (stop - start)].to_numpy(dtype=np.float64)
        else:
            chunk = positions[start:stop]
            values = df.iloc[:, chunk].to_numpy(dtype=np.float64)
        finite = np.isfinite(values)
        valid = finite.sum(axis=0)
        infinite = np.isinf(values).sum(axis=0)
        summed = np.where(finite, values, 0.0).sum(axis=0)
        n_valid[start:stop] = valid
        n_infinite[start:stop] = infinite
        n_missing[start:stop] = values.shape[0] - valid
        np.divide(summed, valid, out=means[start:stop], where=valid > 0)
    return means, n_valid, n_missing, n_infinite


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

    n_cols = len(df.columns)
    diag_column = list(df.columns)
    diag_position = list(range(n_cols))
    diag_gene = [""] * n_cols
    diag_status = ["excluded"] * n_cols
    diag_reason = ["non_gene_schema"] * n_cols
    diag_valid: list = [None] * n_cols
    diag_missing: list = [None] * n_cols
    diag_invalid: list = [None] * n_cols
    diag_infinite: list = [None] * n_cols

    candidates = []
    for position, column in enumerate(df.columns):
        raw = str(column).strip()
        match = GENE_COLUMN_PATTERN.fullmatch(raw)
        explicit = selected is not None and column in selected
        symbol = match.group(1) if match else raw if explicit and SYMBOL_PATTERN.fullmatch(raw) else None
        if symbol:
            diag_gene[position] = symbol
        if selected is not None and not explicit:
            diag_reason[position] = "not_selected"
        if symbol and (selected is None or explicit):
            candidates.append((position, column, symbol))

    counts = Counter(symbol.upper() for _, _, symbol in candidates)
    candidate_positions = [position for position, _, _ in candidates]
    if candidates:
        unique_dtypes = set(df.dtypes.iloc[candidate_positions])
        use_fast = all(
            pd.api.types.is_numeric_dtype(dt) and not pd.api.types.is_bool_dtype(dt)
            for dt in unique_dtypes
        )
    else:
        use_fast = False

    if use_fast:
        means, valid_counts, missing_counts, infinite_counts = _numeric_block_stats(df, candidate_positions)

    n_rows = len(df)
    row_raw, row_gene, row_upper = [], [], []
    row_mean, row_valid, row_missing, row_missing_frac = [], [], [], []
    for index, (position, column, symbol) in enumerate(candidates):
        if use_fast:
            n_val = int(valid_counts[index])
            n_mis = int(missing_counts[index])
            n_inv = 0
            n_inf = int(infinite_counts[index])
            mean_score = float(means[index]) if n_val else np.nan
        else:
            numeric, value_counts = _numeric_values(df.iloc[:, position])
            n_val = value_counts["n_valid"]
            n_mis = value_counts["n_missing"]
            n_inv = value_counts["n_invalid"]
            n_inf = value_counts["n_infinite"]
            mean_score = float(numeric.mean()) if n_val else np.nan
        diag_valid[position] = n_val
        diag_missing[position] = n_mis
        diag_invalid[position] = n_inv
        diag_infinite[position] = n_inf
        symbol_upper = symbol.upper()
        if counts[symbol_upper] > 1:
            diag_reason[position] = "ambiguous_gene_symbol"
        elif not n_val:
            diag_reason[position] = "no_finite_values"
        else:
            diag_status[position] = "included"
            diag_reason[position] = "included"
            row_raw.append(column)
            row_gene.append(symbol)
            row_upper.append(symbol_upper)
            row_mean.append(mean_score)
            row_valid.append(n_val)
            row_missing.append(n_mis)
            row_missing_frac.append(n_mis / n_rows if n_rows else np.nan)

    diagnostics = pd.DataFrame({
        "column": diag_column,
        "position": diag_position,
        "gene": diag_gene,
        "status": diag_status,
        "reason": diag_reason,
        "n_valid": diag_valid,
        "n_missing": diag_missing,
        "n_invalid": diag_invalid,
        "n_infinite": diag_infinite,
    }, columns=DIAGNOSTIC_COLUMNS)

    if not row_gene:
        rankings = pd.DataFrame(columns=RANKING_COLUMNS)
    else:
        rankings = pd.DataFrame({
            "gene_raw": row_raw,
            "gene": row_gene,
            "gene_upper": row_upper,
            "mean_score": row_mean,
            "n_valid": row_valid,
            "n_missing": row_missing,
            "missing_fraction": row_missing_frac,
        })
        rankings = rankings.sort_values("mean_score", kind="stable").reset_index(drop=True)
        rankings["rank"] = np.arange(1, len(rankings) + 1)
        rankings["percentile"] = rankings["rank"] / len(rankings) * 100
        rankings = rankings[RANKING_COLUMNS]
    return rankings, diagnostics


def file_sha256(path) -> str:
    """SHA-256 of a file, read in 1 MiB chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def builtin_provenance_path(path) -> Path:
    """Provenance JSON stored beside a bundled ``.npz`` matrix."""
    return Path(path).with_suffix(".provenance.json")


def read_builtin_matrix(path) -> pd.DataFrame:
    """Load the bundled Chronos matrix and the model columns the app expects.

    ``scores`` is float32 or float64 with NaN for unobserved values. Gene
    labels keep the DepMap ``SYMBOL (Entrez ID)`` form so the normal column
    rules apply. The returned frame is ``depmap_id``, ``cell_line_display_name``,
    ``lineage``, optional ``primary_disease`` when bundled, then one column per gene.
    """
    with np.load(path, allow_pickle=False) as archive:
        scores = np.asarray(archive["scores"])
        genes = archive["genes"].tolist()
        depmap_id = archive["depmap_id"].tolist()
        names = archive["cell_line_display_name"].tolist()
        lineage = archive["lineage"].tolist()
        primary_disease = (
            archive["primary_disease"].tolist()
            if "primary_disease" in archive.files else None
        )
    if scores.ndim != 2 or scores.dtype not in (np.float32, np.float64):
        raise ValueError("Built-in score matrix must be a float32 or float64 array")
    if scores.shape != (len(depmap_id), len(genes)):
        raise ValueError("Built-in score matrix does not match its model or gene labels")
    if not (len(depmap_id) == len(names) == len(lineage)):
        raise ValueError("Built-in model metadata lengths do not match")
    if primary_disease is not None and len(primary_disease) != len(depmap_id):
        raise ValueError("Built-in primary_disease metadata length does not match")
    meta_dict = {
        "depmap_id": pd.Series(depmap_id, dtype="string"),
        "cell_line_display_name": pd.Series(names, dtype="string"),
        "lineage": pd.Series(lineage, dtype="string"),
    }
    if primary_disease is not None:
        meta_dict["primary_disease"] = pd.Series(primary_disease, dtype="string")
    meta_df = pd.DataFrame(meta_dict)
    score_df = pd.DataFrame(scores, columns=pd.Index([str(gene) for gene in genes], dtype="object"))
    return pd.concat([meta_df, score_df], axis=1, copy=False)


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


def resolve_gene_list(rankings: pd.DataFrame, genes: Iterable[str], aliases: Mapping | None = None):
    """Match current symbols, then unique previous or alias symbols.

    ``aliases`` maps an uppercase alias to one uppercase approved symbol.
    A direct symbol match wins. An alias is used only when that approved
    symbol is in the ranking. The same approved symbol is returned once.
    """
    aliases = {str(key).upper(): str(value).upper() for key, value in (aliases or {}).items()}
    lookup = dict(zip(rankings["gene_upper"], rankings["gene"]))
    matched, missing, resolved, seen, matched_seen = [], [], [], set(), set()
    for gene in genes:
        key = str(gene).strip().upper()
        if not key or key in seen:
            continue
        seen.add(key)
        approved = key if key in lookup else aliases.get(key)
        if approved in lookup:
            canonical = lookup[approved]
            if approved != key:
                resolved.append({"input": key, "gene": canonical})
            if canonical not in matched_seen:
                matched.append(canonical)
                matched_seen.add(canonical)
        else:
            missing.append(str(gene).strip())
    return matched, missing, resolved


def cohort_mean_contrast(genes: Iterable[str], cohort_rankings: pd.DataFrame, panel_rankings: pd.DataFrame) -> pd.DataFrame:
    """Selected-scope mean minus the all-lineage mean for each requested gene.

    The difference is descriptive. It is not a test of cancer-type selectivity.
    """
    columns = ["gene", "cohort_mean", "panel_mean", "difference", "cohort_n", "panel_n"]
    cohort = cohort_rankings.set_index("gene", drop=False)
    panel = panel_rankings.set_index("gene", drop=False)
    rows = []
    for gene in genes:
        if gene not in cohort.index or gene not in panel.index:
            continue
        cohort_row = cohort.loc[gene]
        panel_row = panel.loc[gene]
        if isinstance(cohort_row, pd.DataFrame):
            cohort_row = cohort_row.iloc[0]
        if isinstance(panel_row, pd.DataFrame):
            panel_row = panel_row.iloc[0]
        cohort_mean = float(cohort_row["mean_score"])
        panel_mean = float(panel_row["mean_score"])
        rows.append({
            "gene": gene,
            "cohort_mean": cohort_mean,
            "panel_mean": panel_mean,
            "difference": cohort_mean - panel_mean,
            "cohort_n": int(cohort_row["n_valid"]),
            "panel_n": int(panel_row["n_valid"]),
        })
    return pd.DataFrame(rows, columns=columns)


def lineage_preference(frame: pd.DataFrame, threshold: float = -0.5) -> pd.DataFrame:
    """Per cancer type, mean/median score minus that gene's all-lineage average.

    ``frame`` is long-form finite scores with ``gene``, ``lineage`` and
    ``crispr_score``. A negative difference is a lower score in that type.
    ``dependent_fraction`` reports the share of finite cell lines below
    ``threshold``. All metrics are descriptive and not a selectivity test.
    """
    columns = [
        "gene", "lineage", "n_valid", "mean", "panel_mean", "difference",
        "median", "panel_median", "median_difference", "dependent_fraction",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    valid = frame.copy()
    valid["crispr_score"], _ = _numeric_values(valid["crispr_score"])
    valid = valid.dropna(subset=["crispr_score"])
    if valid.empty:
        return pd.DataFrame(columns=columns)
    valid["lineage"] = valid["lineage"].astype("string").str.strip().replace("", pd.NA).fillna("Unknown")
    valid["_dep"] = valid["crispr_score"] < threshold
    grouped = valid.groupby(["gene", "lineage"], sort=False, observed=True).agg(
        n_valid=("crispr_score", "count"),
        mean=("crispr_score", "mean"),
        median=("crispr_score", "median"),
        dependent_fraction=("_dep", "mean"),
    ).reset_index()
    panel_mean = valid.groupby("gene", sort=False)["crispr_score"].mean()
    panel_median = valid.groupby("gene", sort=False)["crispr_score"].median()
    grouped["panel_mean"] = grouped["gene"].map(panel_mean).astype(float)
    grouped["difference"] = grouped["mean"] - grouped["panel_mean"]
    grouped["panel_median"] = grouped["gene"].map(panel_median).astype(float)
    grouped["median_difference"] = grouped["median"] - grouped["panel_median"]
    return grouped[columns]


def find_lineage_column(df: pd.DataFrame):
    """Prefer exact coarse-lineage labels; never select a sub-lineage label."""
    lookup = {str(column).strip().casefold(): column for column in df.columns}
    for name in ("lineage", "oncotreelineage", "oncotree_lineage"):
        if name in lookup:
            return lookup[name]
    return next((column for column in df.columns if "lineage" in str(column).casefold() and "sub" not in str(column).casefold()), None)


def find_disease_column(df: pd.DataFrame):
    """Find an optional primary-disease or subtype metadata column if present."""
    lookup = {str(column).strip().casefold(): column for column in df.columns}
    for name in (
        "primary_disease", "oncotreeprimarydisease", "oncotree_primary_disease",
        "oncotreesubtype", "oncotree_subtype", "disease",
    ):
        if name in lookup:
            return lookup[name]
    return None


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


TOP_CELL_LINE_COLUMNS = [
    "gene", "direction", "rank_in_group", "cell_line_name",
    "cell_line", "lineage", "crispr_score",
]


def top_cell_lines(frame: pd.DataFrame, genes: Iterable[str] | None = None, top_n: int = 10) -> pd.DataFrame:
    """Return the most and least dependent cell lines per gene for wet-lab follow-up."""
    if frame.empty or top_n <= 0:
        return pd.DataFrame(columns=TOP_CELL_LINE_COLUMNS)
    valid = frame.copy()
    valid["crispr_score"], _ = _numeric_values(valid["crispr_score"])
    valid = valid.dropna(subset=["crispr_score"])
    if valid.empty:
        return pd.DataFrame(columns=TOP_CELL_LINE_COLUMNS)
    ordered_genes = list(dict.fromkeys(genes)) if genes is not None else list(dict.fromkeys(valid["gene"]))
    parts = []
    for gene in ordered_genes:
        sub = valid[valid["gene"] == gene]
        if sub.empty:
            continue
        most = sub.sort_values(["crispr_score", "cell_line_name"], ascending=[True, True]).head(top_n).copy()
        most["direction"] = "most_dependent"
        most["rank_in_group"] = np.arange(1, len(most) + 1)
        least = sub.sort_values(["crispr_score", "cell_line_name"], ascending=[False, True]).head(top_n).copy()
        least["direction"] = "least_dependent"
        least["rank_in_group"] = np.arange(1, len(least) + 1)
        parts.extend((most, least))
    if not parts:
        return pd.DataFrame(columns=TOP_CELL_LINE_COLUMNS)
    return pd.concat(parts, ignore_index=True)[TOP_CELL_LINE_COLUMNS]


def gene_correlation_matrix(df: pd.DataFrame, genes: Iterable[str], rankings: pd.DataFrame,
                            method: str = "pearson") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute pairwise gene co-dependency correlation and shared valid sample counts."""
    method = str(method).lower()
    if method not in {"pearson", "spearman"}:
        raise ValueError("Correlation method must be 'pearson' or 'spearman'")
    matched, _ = match_genes(rankings, genes)
    if len(matched) < 2:
        empty = pd.DataFrame(index=matched, columns=matched, dtype=float)
        return empty, empty.astype(int)
    lookup = rankings.set_index("gene_upper")["gene_raw"].to_dict()
    cols = {}
    for gene in matched:
        numeric, _ = _numeric_values(df[lookup[gene.upper()]])
        cols[gene] = numeric.to_numpy(dtype=float)
    matrix = pd.DataFrame(cols)
    corr = matrix.corr(method=method)
    finite_mask = matrix.notna().to_numpy(dtype=int)
    pair_counts = pd.DataFrame(finite_mask.T @ finite_mask, index=matched, columns=matched, dtype=int)
    return corr, pair_counts


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
        preference_rows = provenance.get("lineage_preference")
        if isinstance(preference_rows, list) and preference_rows:
            archive.writestr(
                "lineage_preference.csv",
                pd.DataFrame(preference_rows).to_csv(index=False).encode("utf-8-sig"),
            )
        archive.writestr("metadata.json", json.dumps(provenance, ensure_ascii=False, indent=2, default=_json_default))
        if diagnostics is not None:
            archive.writestr("diagnostics.csv", pd.DataFrame(diagnostics).to_csv(index=False).encode("utf-8-sig"))
    return buffer.getvalue()


def default_builtin_matrix_path() -> Path:
    """Return the path to the optional bundled DepMap matrix (or CRISPR_BUILTIN_MATRIX override)."""
    override = os.environ.get("CRISPR_BUILTIN_MATRIX")
    if override:
        return Path(override)
    return Path(__file__).resolve().parents[1] / "data" / "depmap_chronos.npz"


def load_gene_aliases(path: str | Path | None = None) -> dict[str, str]:
    """Load the unambiguous HGNC previous/alias-to-approved symbol mapping.

    Checks an explicit ``path`` when provided; otherwise reads the packaged
    ``crispr_score_analyzer/data/gene_aliases.json`` resource and falls back to
    ``data/gene_aliases.json`` in a repository checkout. Returns ``{}`` if no
    alias table is available.
    """
    candidates: list[Path] = []
    if path is not None:
        candidates.append(Path(path))
    else:
        try:
            pkg_file = resources.files("crispr_score_analyzer").joinpath("data", "gene_aliases.json")
            if pkg_file.is_file():
                payload = json.loads(pkg_file.read_text(encoding="utf-8"))
                return {str(k).upper(): str(v).upper() for k, v in payload.get("map", {}).items()}
        except (FileNotFoundError, ModuleNotFoundError, OSError):
            pass
        candidates.append(Path(__file__).resolve().parents[1] / "data" / "gene_aliases.json")
    for candidate in candidates:
        if candidate.is_file():
            payload = json.loads(candidate.read_text(encoding="utf-8"))
            return {str(k).upper(): str(v).upper() for k, v in payload.get("map", {}).items()}
    return {}


def filter_gene_suggestions(
    query: str,
    symbols: Sequence[str],
    aliases: Mapping[str, str] | None = None,
    limit: int = 24,
) -> list[dict[str, str | None]]:
    """Deterministically rank symbol and HGNC-alias suggestions for interactive search.

    Prioritizes exact approved symbol matches, then case-insensitive prefix
    matches, then exact/prefix HGNC alias resolutions, and finally substring
    matches. Designed for native text inputs without fragile client-side
    fuzzy combobox filtering.
    """
    if limit <= 0:
        return []
    cleaned = _ENTREZ_SUFFIX_PATTERN.sub("", str(query or "").strip())
    q = cleaned.upper()
    if not q:
        return []

    symbol_by_upper: dict[str, str] = {}
    ordered_symbols: list[tuple[str, str]] = []
    for sym in symbols:
        text = str(sym).strip()
        if not text:
            continue
        upper = text.upper()
        if upper not in symbol_by_upper:
            symbol_by_upper[upper] = text
            ordered_symbols.append((upper, text))

    results: list[dict[str, str | None]] = []
    seen: set[str] = set()

    def add(upper_gene: str, matched_via: str | None = None) -> bool:
        if upper_gene in seen or upper_gene not in symbol_by_upper:
            return False
        seen.add(upper_gene)
        results.append({"gene": symbol_by_upper[upper_gene], "matched_via": matched_via})
        return len(results) >= limit

    # 1. Exact approved symbol match.
    if q in symbol_by_upper and add(q):
        return results

    # 2. Exact HGNC alias match.
    if aliases:
        target = aliases.get(q)
        if target and add(str(target).upper(), matched_via=q):
            return results

    # 3. Prefix matches on approved symbols.
    for upper, _ in ordered_symbols:
        if upper.startswith(q) and add(upper):
            return results

    # 4. Prefix matches on HGNC aliases (sorted for deterministic output).
    if aliases and len(q) >= 2:
        for alias_key in sorted(aliases):
            if alias_key.startswith(q):
                target = str(aliases[alias_key]).upper()
                if add(target, matched_via=alias_key):
                    return results

    # 5. Substring matches on approved symbols.
    for upper, _ in ordered_symbols:
        if q in upper and add(upper):
            return results

    return results

