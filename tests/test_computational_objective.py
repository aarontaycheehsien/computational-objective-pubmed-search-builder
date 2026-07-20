import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "computational_objective.py"
SPEC = importlib.util.spec_from_file_location("computational_objective", MODULE_PATH)
computational_objective = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(computational_objective)


def record(pmid, title, abstract="", mesh=None, relevant=True):
    return {
        "pmid": str(pmid),
        "title": title,
        "abstract": abstract,
        "mesh_headings": [{"name": name} for name in (mesh or [])],
        "relevant": relevant,
    }


class ComputationalObjectiveTests(unittest.TestCase):
    def test_split_references_uses_half_quarter_quarter(self):
        pmids = [str(i) for i in range(1, 9)]
        result = computational_objective.split_pmids(pmids, seed="fixed")

        self.assertEqual(result["development_count"], 4)
        self.assertEqual(result["validation_count"], 2)
        self.assertEqual(result["unseen_count"], 2)
        self.assertEqual(
            set(result["development_pmids"]) | set(result["validation_pmids"]) | set(result["unseen_pmids"]),
            set(pmids),
        )

    def test_term_stats_uses_manual_category_map_and_background(self):
        records = [
            record(1, "Cancer therapy randomized trial", mesh=["Neoplasms"]),
            record(2, "Cancer therapy cohort study", mesh=["Neoplasms"]),
            record(3, "Tumor surgery cohort study", mesh=["Surgical Procedures, Operative"]),
        ]
        split = {"development": ["1", "2", "3"], "validation": [], "unseen": []}
        category_map = {
            "cancer": "condition",
            "tumor": "condition",
            "therapy": "treatment",
            "surgery": "treatment",
            "randomized": "study_type",
            "cohort": "study_type",
        }
        term_rank = {
            "pubmed_total_estimate": 1000,
            "ranked_terms": [
                {"term": "cancer", "field": "tiab", "background_count": 10},
                {"term": "therapy", "field": "tiab", "background_count": 20},
            ],
        }

        result = computational_objective.build_term_statistics(
            records,
            split,
            category_map=category_map,
            term_rank_payload=term_rank,
            require_background=False,
        )

        by_term = {row["term"]: row for row in result["terms"]}
        self.assertEqual(by_term["cancer"]["category"], "condition")
        self.assertEqual(by_term["therapy"]["category"], "treatment")
        self.assertEqual(by_term["cohort"]["category"], "study_type")
        self.assertEqual(by_term["cancer"]["development_df"], 2)
        self.assertEqual(by_term["cancer"]["background_prevalence"], 0.01)

    def test_candidate_filter_and_prune_remove_redundant_terms(self):
        terms = [
            {
                "term": "cancer",
                "field": "tiab",
                "category": "condition",
                "development_df": 3,
                "development_prevalence": 1.0,
                "development_pmids": ["1", "2", "3"],
                "background_prevalence": 0.01,
                "suggested_layer": "cancer[tiab]",
            },
            {
                "term": "tumor",
                "field": "tiab",
                "category": "condition",
                "development_df": 2,
                "development_prevalence": 0.667,
                "development_pmids": ["1", "2"],
                "background_prevalence": 0.01,
                "suggested_layer": "tumor[tiab]",
            },
            {
                "term": "therapy",
                "field": "tiab",
                "category": "treatment",
                "development_df": 3,
                "development_prevalence": 1.0,
                "development_pmids": ["1", "2", "3"],
                "background_prevalence": 0.02,
                "suggested_layer": "therapy[tiab]",
            },
        ]
        term_stats = {"terms": terms}

        candidates = computational_objective.candidate_terms_for_grid(
            term_stats,
            development_threshold=0.20,
            population_threshold=0.02,
            mesh_count=0,
        )
        pruned = computational_objective.prune_terms(candidates, ["1", "2", "3"])

        kept = {row["term"] for row in pruned["kept_terms"]}
        removed = {row["term"] for row in pruned["removed_terms"]}
        self.assertEqual(kept, {"cancer", "therapy"})
        self.assertEqual(removed, {"tumor"})
        self.assertIn("AND", computational_objective.build_query(pruned["kept_terms"]))

    def test_derive_query_returns_selected_parameters_and_query(self):
        records = [
            record(1, "Cancer therapy alpha"),
            record(2, "Cancer therapy beta"),
            record(3, "Tumor surgery gamma"),
            record(4, "Tumor surgery delta"),
            record(5, "Cancer therapy epsilon"),
            record(6, "Tumor surgery zeta"),
            record(7, "Cancer therapy eta"),
            record(8, "Tumor surgery theta"),
        ]
        split = {"development": ["1", "2", "3", "4"], "validation": ["5", "6"], "unseen": ["7", "8"]}
        category_map = {
            "cancer": "condition",
            "tumor": "condition",
            "therapy": "treatment",
            "surgery": "treatment",
        }

        result = computational_objective.derive_query(
            records,
            split,
            category_map=category_map,
            term_rank_payload=None,
            require_background=False,
            development_thresholds=[0.25, 0.50],
            population_thresholds=[0.02],
            mesh_counts=[0],
            objective="recall",
        )

        self.assertEqual(result["selected_parameters"]["population_threshold"], 0.02)
        self.assertGreater(result["kept_count"], 0)
        self.assertIn("[tiab]", result["query"])
        self.assertIsNotNone(result["unseen_evaluation"])

    def test_cli_split_and_derive_write_outputs(self):
        records = [
            record(1, "Cancer therapy alpha"),
            record(2, "Cancer therapy beta"),
            record(3, "Tumor surgery gamma"),
            record(4, "Tumor surgery delta"),
        ]
        split = {"development_pmids": ["1", "2"], "validation_pmids": ["3"], "unseen_pmids": ["4"]}
        category_map = {"condition": ["cancer", "tumor"], "treatment": ["therapy", "surgery"]}

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            records_path = tmp_path / "records.json"
            split_path = tmp_path / "split.json"
            map_path = tmp_path / "map.json"
            output_path = tmp_path / "derived.json"
            records_path.write_text(json.dumps({"records": records}), encoding="utf-8")
            split_path.write_text(json.dumps(split), encoding="utf-8")
            map_path.write_text(json.dumps(category_map), encoding="utf-8")

            rc = computational_objective.main(
                [
                    "derive-query",
                    "--records-json",
                    str(records_path),
                    "--split-json",
                    str(split_path),
                    "--category-map-json",
                    str(map_path),
                    "--development-thresholds",
                    "0.5",
                    "--population-thresholds",
                    "0.02",
                    "--mesh-counts",
                    "0",
                    "--output",
                    str(output_path),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(output_path.exists())
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["operation"], "derive-query")


if __name__ == "__main__":
    unittest.main()
