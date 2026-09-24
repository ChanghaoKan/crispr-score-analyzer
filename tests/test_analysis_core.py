import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import numpy as np
import pandas as pd

from analysis_core import (
    analyze_dataset, find_lineage_column, get_lineage_frame, lineage_summary,
    make_analysis_bundle, match_genes, parse_gene_text, parse_gene_upload, read_score_csv,
)


class GeneInputTests(unittest.TestCase):
    def test_delimiters_case_and_order(self):
        parsed = parse_gene_text("myc\tTP53，KIF18A;MYC\n brca1；tp53")
        self.assertEqual(parsed, {"genes": ["MYC", "TP53", "KIF18A", "BRCA1"], "input_count": 6, "duplicates": ["MYC", "TP53"]})

    def test_csv_does_not_discard_first_gene(self):
        self.assertEqual(parse_gene_upload(b"MYC\nTP53\n", "genes.csv")["genes"], ["MYC", "TP53"])
        self.assertEqual(parse_gene_upload(b"MYC,1\nTP53,2\n", "genes.csv")["genes"], ["MYC", "TP53"])

    def test_csv_header_bom_quotes_and_tab(self):
        self.assertEqual(parse_gene_upload('\ufeffGene Symbol,score\n"MYC",1\nTP53,2'.encode(), "genes.csv")["genes"], ["MYC", "TP53"])
        self.assertEqual(parse_gene_upload(b"symbol\tother\nMYC\t1\nTP53\t2", "genes.tsv")["genes"], ["MYC", "TP53"])

    def test_text_and_empty(self):
        self.assertEqual(parse_gene_upload(b"genes\nmyc tp53", "genes.txt")["genes"], ["MYC", "TP53"])
        self.assertEqual(parse_gene_text(""), {"genes": [], "input_count": 0, "duplicates": []})


class RankingTests(unittest.TestCase):
    def test_extra_fields_cannot_silently_shift_scores(self):
        shifted = b"ModelID,lineage,MYC (4609)\nACH-1,Liver,-1.2,-0.2\nACH-2,Lung,-1.4,0.3\n"
        with self.assertRaisesRegex(ValueError, "more fields than the header"):
            read_score_csv(shifted)

    def test_csv_loader_preserves_duplicate_headers_for_exclusion(self):
        source = b"ModelID,age,MYC (4609),MYC (4609),TP53 (7157)\nACH-1,22,-1,-2,0\n"
        for representation in (source, io.BytesIO(source)):
            frame = read_score_csv(representation)
            self.assertEqual(frame.columns.tolist().count("MYC (4609)"), 2)
            ranks, diagnostics = analyze_dataset(frame)
            self.assertEqual(ranks.gene.tolist(), ["TP53"])
            self.assertEqual((diagnostics.reason == "ambiguous_gene_symbol").sum(), 2)
            self.assertEqual(diagnostics.set_index("column").loc["age", "reason"], "non_gene_schema")

    def test_csv_loader_paths_bom_and_no_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scores.csv"
            path.write_bytes('\ufeffModelID,"MYC (4609)"\nACH-1,-1\n'.encode())
            frame = read_score_csv(path)
            self.assertEqual(frame.columns.tolist(), ["ModelID", "MYC (4609)"])
            self.assertEqual(frame.iloc[0, 1], -1)
        self.assertTrue(read_score_csv(b"ModelID,MYC (4609)\n").empty)
        with self.assertRaises(pd.errors.EmptyDataError):
            read_score_csv(b"")

    def test_csv_loader_rejects_duplicate_lineage_metadata(self):
        for header in ("lineage,lineage", "OncotreeLineage,oncotreelineage", "lineage_1,lineage_1"):
            with self.assertRaisesRegex(ValueError, "Duplicate metadata columns"):
                read_score_csv(f"{header},MYC (4609)\nLiver,Lung,-1\n".encode())

    def test_csv_loader_rejects_duplicate_cell_line_ids_and_names(self):
        for header in ("ModelID,ModelID", "depmap_id,depmap_id", "cell_line_display_name,cell_line_display_name"):
            with self.assertRaisesRegex(ValueError, "Duplicate metadata columns"):
                read_score_csv(f"{header},MYC (4609)\nACH-1,ACH-2,-1\n".encode())
        frame = read_score_csv(b"MYC,MYC,age,age\n-1,-2,22,23\n")
        ranks, diagnostics = analyze_dataset(frame, ["MYC"])
        self.assertTrue(ranks.empty)
        self.assertEqual((diagnostics.reason == "ambiguous_gene_symbol").sum(), 2)

    def test_small_and_constant_genes_retained_metadata_excluded(self):
        df = pd.DataFrame({"ModelID": ["A", "B"], "age": [20, 30], "MYC (4609)": [-8, -8], "TP53 (7157)": [3, 3], "PLAIN": [-1, -2]})
        ranks, diagnostics = analyze_dataset(df)
        self.assertEqual(ranks["gene"].tolist(), ["MYC", "TP53"])
        self.assertEqual(ranks["n_valid"].tolist(), [2, 2])
        self.assertEqual(ranks["mean_score"].tolist(), [-8, 3])
        self.assertEqual(ranks["percentile"].tolist(), [50, 100])
        self.assertEqual(diagnostics.set_index("column").loc["age", "reason"], "non_gene_schema")

    def test_invalid_infinite_missing_and_no_finite_values(self):
        df = pd.DataFrame({"MYC (4609)": ["-1", "bad", np.inf, None, ""], "EMPTY (1)": [np.nan] * 5})
        ranks, diagnostics = analyze_dataset(df)
        self.assertEqual(ranks.loc[0, "mean_score"], -1)
        self.assertEqual(ranks.loc[0, "n_valid"], 1)
        self.assertEqual(ranks.loc[0, "n_missing"], 4)
        self.assertEqual(ranks.loc[0, "missing_fraction"], .8)
        self.assertEqual(diagnostics.loc[0, "n_invalid"], 1)
        self.assertEqual(diagnostics.loc[0, "n_infinite"], 1)
        self.assertEqual(diagnostics.loc[1, "reason"], "no_finite_values")

    def test_duplicate_symbols_never_silently_overwrite(self):
        df = pd.DataFrame({"MYC (4609)": [-1], "myc (999)": [-2], "TP53 (7157)": [0]})
        ranks, diagnostics = analyze_dataset(df)
        self.assertEqual(ranks["gene"].tolist(), ["TP53"])
        self.assertEqual(diagnostics["reason"].tolist()[:2], ["ambiguous_gene_symbol"] * 2)
        duplicate_names = pd.DataFrame([[-1, -2]], columns=["MYC (4609)", "MYC (4609)"])
        self.assertTrue(analyze_dataset(duplicate_names)[0].empty)

    def test_explicit_symbols_and_unknown_selection(self):
        df = pd.DataFrame({"MYC": [0], "TP53": [-1], "age": [50]})
        ranks, _ = analyze_dataset(df, ["MYC", "TP53"])
        self.assertEqual(ranks["gene"].tolist(), ["TP53", "MYC"])
        with self.assertRaises(ValueError):
            analyze_dataset(df, ["absent"])

    def test_schema_is_strict_and_empty_is_supported(self):
        df = pd.DataFrame({"MYC (text)": [1], "MYC (4609) extra": [2], "something": [3]})
        self.assertTrue(analyze_dataset(df)[0].empty)
        self.assertTrue(analyze_dataset(df.iloc[:0])[0].empty)

    def test_matching_preserves_input_order(self):
        ranks, _ = analyze_dataset(pd.DataFrame({"MYC (4609)": [-1], "TP53 (7157)": [-2]}))
        self.assertEqual(match_genes(ranks, ["myc", "unknown", "tp53", "MYC"]), (["MYC", "TP53"], ["unknown"]))


class LineageAndExportTests(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            "ModelID": ["ACH-1", "ACH-2", "ACH-3"],
            "OncotreeSubtype": ["a", "b", "c"], "OncotreeLineage": ["Lung", "Lung", None],
            "MYC (4609)": [-1, -3, np.inf], "TP53 (7157)": [0, np.nan, 2],
        })
        self.ranks, self.diagnostics = analyze_dataset(self.df)

    def test_lineage_lookup_exact_preference_and_sublineage_exclusion(self):
        self.assertEqual(find_lineage_column(self.df), "OncotreeLineage")
        self.assertEqual(find_lineage_column(self.df.assign(lineage="Liver")), "lineage")
        self.assertIsNone(find_lineage_column(pd.DataFrame({"sublineage": []})))

    def test_long_data_keeps_ids_unknown_and_only_finite(self):
        frame = get_lineage_frame(self.df, ["TP53", "MYC"], self.ranks)
        self.assertEqual(frame["gene"].tolist(), ["TP53", "TP53", "MYC", "MYC"])
        self.assertEqual(frame["cell_line"].tolist(), ["ACH-1", "ACH-3", "ACH-1", "ACH-2"])
        self.assertEqual(frame["cell_line_name"].tolist(), frame["cell_line"].tolist())
        self.assertEqual(frame.loc[1, "lineage"], "Unknown")
        summary = lineage_summary(frame)
        myc = summary[summary["gene"] == "MYC"].iloc[0]
        self.assertEqual(myc["n_valid"], 2)
        self.assertEqual(myc["median"], -2)

    def test_cell_line_names_kept_alongside_ids_with_fallback(self):
        data = self.df.assign(cell_line_display_name=["HepG2", " ", None])
        frame = get_lineage_frame(data, ["TP53", "MYC"], self.ranks)
        self.assertEqual(frame["cell_line_name"].tolist(), ["HepG2", "ACH-3", "HepG2", "ACH-2"])
        self.assertEqual(frame["cell_line"].tolist(), ["ACH-1", "ACH-3", "ACH-1", "ACH-2"])

    def test_bundle_contains_complete_data_and_reproducibility(self):
        plot = get_lineage_frame(self.df, ["TP53", "MYC"], self.ranks)
        bundle = make_analysis_bundle(self.ranks, ["TP53", "MYC", "ABSENT"], plot, {"dataset_version": "test", "parameters": {"color": "red"}}, self.diagnostics)
        with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
            self.assertEqual(set(archive.namelist()), {"selected_results.csv", "full_rankings.csv", "plot_data.csv", "metadata.json", "diagnostics.csv"})
            selected = pd.read_csv(archive.open("selected_results.csv"))
            self.assertEqual(selected["gene"].tolist(), ["TP53", "MYC"])
            full = pd.read_csv(archive.open("full_rankings.csv"))
            self.assertEqual(len(full), len(self.ranks))
            metadata = json.loads(archive.read("metadata.json"))
            self.assertEqual(metadata["dataset_version"], "test")
            self.assertEqual(metadata["analysis"]["n_cell_lines"], 3)
            self.assertEqual(metadata["analysis"]["n_plot_rows"], 4)
            self.assertEqual(metadata["analysis"]["selected_genes_unmatched"], ["ABSENT"])

    def test_empty_bundle_and_absent_lineage(self):
        ranks, diag = analyze_dataset(pd.DataFrame())
        self.assertTrue(get_lineage_frame(pd.DataFrame(), [], ranks).empty)
        self.assertTrue(lineage_summary(pd.DataFrame()).empty)
        with zipfile.ZipFile(io.BytesIO(make_analysis_bundle(ranks, [], pd.DataFrame(), {}, diag))) as archive:
            self.assertEqual(json.loads(archive.read("metadata.json"))["analysis"]["n_ranked_genes"], 0)

    def test_no_lineage_still_exports_selected_gene_scores(self):
        frame = get_lineage_frame(self.df.drop(columns="OncotreeLineage"), ["MYC"], self.ranks)
        self.assertEqual(frame["lineage"].tolist(), ["Unknown", "Unknown"])
        self.assertEqual(frame["crispr_score"].tolist(), [-1, -3])


if __name__ == "__main__":
    unittest.main()
