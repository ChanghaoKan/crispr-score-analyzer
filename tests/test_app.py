"""Offline integration checks for submitted Streamlit results and data changes."""

import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import streamlit as st
from streamlit.testing.v1 import AppTest


APP_PATH = Path(__file__).resolve().parents[1] / "app.py"
GENES = [
    "MYC", "PTEN", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6",
    "E2F7", "E2F8", "CDK1", "CDK2", "CCNB1", "CCND1", "CCNE1", "PLK1", "AURKA",
]


def score_csv(offset=0.0):
    """Small complete matrix, with cohort-specific means and valid gene headers."""
    content = io.StringIO()
    writer = csv.writer(content)
    writer.writerow(["ModelID", "lineage", *[f"{gene} ({i + 1})" for i, gene in enumerate(GENES)]])
    for row in range(12):
        writer.writerow([
            f"ACH-{row + 1:06d}", "Liver" if row < 6 else "Lung",
            *[round(offset - 1.6 + i * 0.12 + row * 0.01, 4) for i in range(len(GENES))],
        ])
    return content.getvalue().encode("utf-8")


def built_in_score_csv():
    """Mirror the pinned public file: seven metadata columns, then plain symbols."""
    content = io.StringIO()
    writer = csv.writer(content)
    writer.writerow([
        "depmap_id", "cell_line_display_name", "lineage_1", "lineage_2",
        "lineage_3", "lineage_6", "lineage_4", *GENES,
    ])
    rows = list(csv.reader(io.StringIO(score_csv().decode("utf-8"))))
    for row in rows[1:]:
        writer.writerow([row[0], "Cell " + row[0], row[1], "Subtype", "Cancer", "Type6", "Type4", *row[2:]])
    return content.getvalue().encode("utf-8")


class Upload(io.BytesIO):
    def __init__(self, content, name="custom.csv"):
        super().__init__(content)
        self.name = name
        self.type = "text/csv"


class AppSubmissionTests(unittest.TestCase):
    def setUp(self):
        st.cache_data.clear()
        st.cache_resource.clear()
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.data_path = Path(self.tempdir.name) / "scores.csv"
        self.data_path.write_bytes(built_in_score_csv())
        self.current_upload = None
        original_uploader = st.file_uploader

        def upload_fixture(*args, **kwargs):
            if kwargs.get("key") == "dataset_file":
                return self.current_upload
            return original_uploader(*args, **kwargs)

        hf_patch = patch("huggingface_hub.hf_hub_download", return_value=str(self.data_path))
        upload_patch = patch("streamlit.file_uploader", side_effect=upload_fixture)
        hf_patch.start()
        upload_patch.start()
        self.addCleanup(hf_patch.stop)
        self.addCleanup(upload_patch.stop)
        self.app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
        self.assert_healthy()

    def tearDown(self):
        st.cache_data.clear()
        st.cache_resource.clear()

    def assert_healthy(self):
        self.assertEqual(len(self.app.exception), 0, [e.message for e in self.app.exception])

    def result(self, key):
        try:
            return self.app.session_state[key]
        except KeyError:
            return None

    def run_rank(self, genes="MYC"):
        self.app.text_area(key="rank_genes").set_value(genes)
        self.app.button(key="rank_run").click().run()
        self.assert_healthy()
        return self.result("rank_result")

    def switch_view(self, view):
        self.app.segmented_control(key="analysis_tab").set_value(view).run()
        self.assert_healthy()

    def select_custom_source(self):
        self.app.segmented_control(key="dataset_source").set_value("custom").run()
        self.assert_healthy()

    def test_run_commits_inputs_and_draft_preserves_last_result(self):
        first = self.run_rank("MYC")
        self.assertEqual(first["genes"], ["MYC"])
        first_ranking = first["rankings"].copy(deep=True)
        first_signature = first["input_signature"]

        self.app.text_area(key="rank_genes").set_value("PTEN").run()
        self.assert_healthy()
        self.assertEqual(self.result("rank_result")["genes"], ["MYC"])
        self.assertTrue(self.result("rank_result")["rankings"].equals(first_ranking))
        self.assertEqual(self.result("rank_result")["input_signature"], first_signature)
        self.assertTrue(any("Inputs changed" in warning.value for warning in self.app.warning))

        self.app.button(key="rank_run").click().run()
        self.assert_healthy()
        self.assertEqual(self.result("rank_result")["genes"], ["PTEN"])
        self.assertNotEqual(self.result("rank_result")["input_signature"], first_signature)
        self.assertFalse(any("Inputs changed" in warning.value for warning in self.app.warning))

    def test_unmatched_submission_keeps_result_and_warns(self):
        original = self.run_rank("MYC")
        self.run_rank("NOT_A_REAL_GENE")
        self.assertEqual(self.result("rank_result")["genes"], ["MYC"])
        self.assertEqual(self.result("rank_result")["metadata"]["created_at_utc"],
                         original["metadata"]["created_at_utc"])
        self.assertTrue(any("No matched genes" in warning.value for warning in self.app.warning))

    def test_example_and_clear_change_draft_without_replacing_result(self):
        cases = [
            ("rank", {"rank_genes": "MYC"}, {"rank_genes": "E2F1\nE2F2\nE2F3"}),
            ("box", {"box_genes": "MYC"}, {"box_genes": "E2F1\nE2F2"}),
            ("multi", {"multi_bg": "MYC", "multi_hl": "PTEN"},
             {"multi_bg": "CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1", "multi_hl": "PLK1\nAURKA"}),
        ]
        for view, initial, example in cases:
            with self.subTest(view=view):
                self.switch_view(view)
                for key, text in initial.items():
                    self.app.text_area(key=key).set_value(text)
                self.app.button(key=f"{view}_run").click().run()
                self.assert_healthy()
                committed = self.result(f"{view}_result")

                self.app.button(key=f"{view}_example").click().run()
                self.assert_healthy()
                for key, text in example.items():
                    self.assertEqual(self.app.text_area(key=key).value, text)
                self.assertEqual(self.result(f"{view}_result")["genes"], committed["genes"])
                self.assertTrue(any("Inputs changed" in warning.value for warning in self.app.warning))

                self.app.button(key=f"{view}_clear").click().run()
                self.assert_healthy()
                for key in initial:
                    self.assertEqual(self.app.text_area(key=key).value, "")
                self.assertEqual(self.result(f"{view}_result")["genes"], committed["genes"])
                self.assertTrue(any("Inputs changed" in warning.value for warning in self.app.warning))

    def test_submitting_another_tab_does_not_submit_rank_draft(self):
        self.run_rank("MYC")
        self.app.text_area(key="rank_genes").set_value("PTEN").run()
        self.switch_view("multi")
        self.app.text_area(key="multi_bg").set_value("CDK1 CDK2")
        self.app.text_area(key="multi_hl").set_value("MYC")
        self.app.button(key="multi_run").click().run()
        self.assert_healthy()
        self.assertEqual(self.result("rank_result")["genes"], ["MYC"])
        self.assertEqual(set(self.result("multi_result")["genes"]), {"CDK1", "CDK2", "MYC"})
        self.switch_view("rank")
        self.assertEqual(self.app.text_area(key="rank_genes").value, "PTEN")
        self.assertEqual(self.result("rank_result")["genes"], ["MYC"])

    def test_cohort_change_invalidates_result_until_run(self):
        self.run_rank("MYC")
        self.app.multiselect(key="cohort").set_value(["Liver"]).run()
        self.assert_healthy()
        self.assertIsNone(self.result("rank_result"))
        result = self.run_rank("MYC")
        self.assertEqual(result["genes"], ["MYC"])
        # The cohort must change the scientific result, not only a title.
        myc = result["rankings"].loc[result["rankings"]["gene"] == "MYC"].iloc[0]
        self.assertAlmostEqual(float(myc["mean_score"]), -1.575)

    def test_new_upload_invalidates_result_and_uses_new_scores(self):
        self.run_rank("MYC")
        self.select_custom_source()
        self.assertIsNone(self.result("rank_result"))
        self.assertEqual(len(self.app.get("plotly_chart")), 0)
        self.current_upload = Upload(score_csv(offset=0.5))
        self.app.run()
        self.assert_healthy()
        self.assertIsNone(self.result("rank_result"))
        result = self.run_rank("MYC")
        myc = result["rankings"].loc[result["rankings"]["gene"] == "MYC"].iloc[0]
        self.assertAlmostEqual(float(myc["mean_score"]), -1.045)

    def test_plain_symbol_upload_requires_explicit_column_mapping(self):
        self.run_rank("MYC")
        self.current_upload = Upload(built_in_score_csv(), "plain_symbols.csv")
        self.select_custom_source()
        self.assert_healthy()
        self.assertIsNone(self.result("rank_result"))
        self.assertGreater(len(self.app.warning), 0)
        self.app.checkbox(key="custom_mapping").set_value(True).run()
        self.app.multiselect(key="gene_columns").set_value(["MYC", "PTEN"]).run()
        self.assert_healthy()
        result = self.run_rank("MYC")
        self.assertEqual(set(result["rankings"]["gene"]), {"MYC", "PTEN"})

    def test_invalid_upload_shows_error_without_retaining_results(self):
        self.run_rank("MYC")
        self.select_custom_source()
        invalid_inputs = {
            "no_gene_columns": b"ModelID,lineage,notes\nACH-1,Liver,hello\n",
            "malformed_csv": b"ModelID,lineage,MYC (4609)\nACH-1,Liver,-1\nACH-2,Lung,-0.5,extra\n",
            "shifted_columns": b"ModelID,lineage,MYC (4609)\nACH-1,Liver,-1.2,-0.2\nACH-2,Lung,-1.4,0.3\n",
            "invalid_encoding": b"\xff\xfe\x00unreadable\x00",
        }
        for name, content in invalid_inputs.items():
            with self.subTest(name=name):
                self.current_upload = Upload(content, name + ".csv")
                self.app.run()
                self.assert_healthy()
                messages = self.app.error if name != "no_gene_columns" else self.app.warning
                self.assertGreater(len(messages), 0)
                self.assertIsNone(self.result("rank_result"))

    def test_language_switch_keeps_committed_result_and_localizes_run(self):
        original = self.run_rank("MYC")
        self.app.selectbox(key="lang_select").set_value("zh").run()
        self.assert_healthy()
        self.assertEqual(self.result("rank_result")["genes"], ["MYC"])
        self.assertEqual(len(self.app.warning), 0)
        self.app.selectbox(key="theme_select").set_value("dark").run()
        self.assert_healthy()
        self.assertEqual(self.result("rank_result")["input_signature"], original["input_signature"])
        self.assertEqual(len(self.app.warning), 0)
        for view, key in [("rank", "rank_run"), ("box", "box_run"), ("multi", "multi_run")]:
            self.switch_view(view)
            self.assertIn("生成", self.app.button(key=key).label)
        self.app.selectbox(key="lang_select").set_value("en").run()
        self.assert_healthy()
        self.switch_view("rank")
        self.assertIn("Run", self.app.button(key="rank_run").label)

    def test_more_than_eight_box_genes_requires_explicit_selection(self):
        requested = ["PTEN", "MYC", "E2F8", "E2F7", "E2F6", "E2F5", "E2F4", "E2F3", "E2F2", "E2F1"]
        self.switch_view("box")
        self.app.text_area(key="box_genes").set_value("\n".join(requested))
        self.app.button(key="box_run").click().run()
        self.assert_healthy()
        self.assertIsNone(self.result("box_result"))
        selector = self.app.multiselect(key="box_display_genes")
        self.assertEqual(selector.options, requested)
        selected = requested[-3:]
        selector.set_value(selected)
        self.app.button(key="box_selection_run").click().run()
        self.assert_healthy()
        self.assertEqual(self.result("box_result")["genes"], selected)

    def test_editing_box_draft_discards_pending_selection_and_keeps_saved_plot(self):
        self.switch_view("box")
        self.app.text_area(key="box_genes").set_value("MYC")
        self.app.button(key="box_run").click().run()
        self.assert_healthy()
        committed = self.result("box_result")

        self.app.text_area(key="box_genes").set_value(" ".join(GENES[:10]))
        self.app.button(key="box_run").click().run()
        self.assert_healthy()
        self.assertIsNotNone(self.result("box_candidates"))
        self.assertEqual(self.result("box_result")["genes"], ["MYC"])

        self.app.text_area(key="box_genes").set_value("PTEN").run()
        self.assert_healthy()
        self.assertIsNone(self.result("box_candidates"))
        self.assertEqual(self.result("box_result")["input_signature"], committed["input_signature"])
        self.assertTrue(any("Inputs changed" in warning.value for warning in self.app.warning))

    def test_box_selection_uses_valid_counts_and_preserves_cell_identity(self):
        rows = list(csv.reader(io.StringIO(score_csv().decode("utf-8"))))
        rows[0].insert(1, "cell_line_display_name")
        for i, row in enumerate(rows[1:], 1):
            row.insert(1, f"Cell {i}")
        rows[1][3] = ""  # MYC score missing for ACH-000001.
        rows[2][2] = ""  # ACH-000002 has no lineage label.
        content = io.StringIO()
        csv.writer(content).writerows(rows)
        self.current_upload = Upload(content.getvalue().encode("utf-8"), "missing_scores.csv")
        self.select_custom_source()
        self.assert_healthy()
        self.switch_view("box")
        self.app.text_area(key="box_genes").set_value("MYC\nPTEN")
        self.app.multiselect(key="box_lineages").set_value(["Liver"])
        self.app.checkbox(key="box_points").set_value(True)
        self.app.button(key="box_run").click().run()
        self.assert_healthy()

        frame = self.result("box_result")["plot_data"]
        self.assertEqual(set(frame["lineage"]), {"Liver"})
        self.assertEqual(frame.groupby("gene").size().to_dict(), {"MYC": 4, "PTEN": 5})
        self.assertFalse(frame["crispr_score"].isna().any())
        myc = frame.loc[frame["gene"] == "MYC"]
        self.assertEqual(set(myc["cell_line"]), {f"ACH-{i:06d}" for i in range(3, 7)})
        self.assertEqual(set(myc["cell_line_name"]), {f"Cell {i}" for i in range(3, 7)})

        figure = json.loads(self.app.get("plotly_chart")[0].proto.spec)
        traces = figure["data"]
        self.assertEqual([set(trace["y"]) for trace in traces], [{"Liver (n=4)"}, {"Liver (n=5)"}])
        self.assertTrue(all(trace["boxpoints"] == "all" for trace in traces))
        self.assertEqual(traces[0]["customdata"][0], ["ACH-000003", "Cell 3"])

        self.app.multiselect(key="box_lineages").set_value(["Unknown"])
        self.app.button(key="box_run").click().run()
        self.assert_healthy()
        unknown = self.result("box_result")["plot_data"]
        self.assertEqual(set(unknown["lineage"]), {"Unknown"})
        self.assertEqual(set(unknown["cell_line"]), {"ACH-000002"})
        self.assertEqual(unknown.groupby("gene").size().to_dict(), {"MYC": 1, "PTEN": 1})


if __name__ == "__main__":
    unittest.main()
