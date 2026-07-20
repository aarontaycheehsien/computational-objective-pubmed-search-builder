#!/usr/bin/env python3
"""Computational objective PubMed query generation helpers.

This is a practical, deterministic implementation of the Scells et al. style
workflow for offline artifacts. PubMed fetching, term-rank background counts,
and final query execution remain delegated to pubmed_tool.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


DEFAULT_SEED = "computational-objective-pubmed-v1"
DEFAULT_DEVELOPMENT_THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
DEFAULT_POPULATION_THRESHOLDS = [0.001, 0.01, 0.02, 0.05, 0.10, 0.20]
DEFAULT_MESH_COUNTS = [0, 1, 5, 10, 15, 20, 25]
VALID_CATEGORIES = ("condition", "treatment", "study_type")
TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9-]{2,}")
SPLIT_TOKEN_RE = re.compile(r"[\s,;]+")

STOPWORDS = {
    "about",
    "above",
    "after",
    "again",
    "against",
    "also",
    "among",
    "and",
    "any",
    "are",
    "because",
    "been",
    "before",
    "being",
    "between",
    "both",
    "but",
    "can",
    "could",
    "during",
    "each",
    "few",
    "for",
    "from",
    "further",
    "had",
    "has",
    "have",
    "having",
    "here",
    "how",
    "into",
    "its",
    "itself",
    "may",
    "more",
    "most",
    "not",
    "our",
    "out",
    "over",
    "own",
    "patient",
    "patients",
    "result",
    "results",
    "same",
    "should",
    "some",
    "such",
    "than",
    "that",
    "the",
    "their",
    "them",
    "then",
    "there",
    "these",
    "they",
    "this",
    "those",
    "through",
    "under",
    "until",
    "using",
    "very",
    "was",
    "were",
    "what",
    "when",
    "where",
    "which",
    "while",
    "who",
    "with",
}

CONDITION_HINTS = {
    "abnormality",
    "abnormalities",
    "adenocarcinoma",
    "cancer",
    "carcinoma",
    "condition",
    "disease",
    "disorder",
    "infection",
    "injury",
    "lesion",
    "malignancy",
    "neoplasm",
    "neoplasms",
    "pain",
    "syndrome",
    "tumor",
    "tumour",
}
TREATMENT_HINTS = {
    "assay",
    "device",
    "diagnosis",
    "diagnostic",
    "drug",
    "implant",
    "intervention",
    "procedure",
    "radiotherapy",
    "screening",
    "surgery",
    "test",
    "testing",
    "tests",
    "therapy",
    "treatment",
}
STUDY_TYPE_HINTS = {
    "accuracy",
    "case-control",
    "cohort",
    "controlled",
    "cross-sectional",
    "diagnostic",
    "meta-analysis",
    "prospective",
    "randomized",
    "retrospective",
    "review",
    "sensitivity",
    "specificity",
    "trial",
}


class ComputationalObjectiveError(Exception):
    pass


def read_json(path: str | Path) -> object:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except OSError as exc:
        raise ComputationalObjectiveError(f"Could not read JSON file: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ComputationalObjectiveError(f"Could not parse JSON file: {path}: {exc}") from exc


def write_json(payload: object, output: str | None = None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if output:
        Path(output).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)


def stable_digest(seed: str, value: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def pmid_tokens(values: Iterable[str]) -> Iterable[str]:
    for value in values:
        for token in SPLIT_TOKEN_RE.split(str(value).strip()):
            if token:
                yield token


def normalize_pmids(values: Iterable[str]) -> dict[str, object]:
    pmids: list[str] = []
    seen: set[str] = set()
    malformed: list[str] = []
    duplicates: list[str] = []
    for token in pmid_tokens(values):
        if not token.isdigit() or int(token) == 0:
            malformed.append(token)
            continue
        normalized = str(int(token))
        if normalized in seen:
            duplicates.append(normalized)
            continue
        seen.add(normalized)
        pmids.append(normalized)
    return {
        "pmids": pmids,
        "pmid_count": len(pmids),
        "duplicates": duplicates,
        "malformed": malformed,
    }


def split_pmids(pmids: Iterable[str], *, seed: str = DEFAULT_SEED) -> dict[str, object]:
    normalized = normalize_pmids(pmids)
    ordered = list(normalized["pmids"])  # type: ignore[arg-type]
    count = len(ordered)
    ranked = sorted(ordered, key=lambda pmid: stable_digest(seed, pmid))
    if count == 0:
        development: list[str] = []
        validation: list[str] = []
        unseen: list[str] = []
        status = "empty"
    elif count == 1:
        development, validation, unseen = ordered[:], [], []
        status = "insufficient_for_validation"
    elif count == 2:
        development, validation, unseen = [ranked[0]], [ranked[1]], []
        status = "insufficient_for_unseen"
    elif count == 3:
        development, validation, unseen = [ranked[0]], [ranked[1]], [ranked[2]]
        status = "small_sample"
    else:
        development_count = max(1, count // 2)
        validation_count = max(1, (count - development_count) // 2)
        development_set = set(ranked[:development_count])
        validation_set = set(ranked[development_count : development_count + validation_count])
        unseen_set = set(ranked[development_count + validation_count :])
        development = [pmid for pmid in ordered if pmid in development_set]
        validation = [pmid for pmid in ordered if pmid in validation_set]
        unseen = [pmid for pmid in ordered if pmid in unseen_set]
        status = "ok" if development_count >= 25 else "development_below_paper_minimum"
    return {
        "operation": "split-references",
        "seed": seed,
        "status": status,
        "pmid_count": count,
        "development_pmids": development,
        "development_count": len(development),
        "validation_pmids": validation,
        "validation_count": len(validation),
        "unseen_pmids": unseen,
        "unseen_count": len(unseen),
        "duplicates": normalized["duplicates"],
        "malformed": normalized["malformed"],
        "note": "Default split is 50% development, 25% validation, and 25% unseen when enough PMIDs exist.",
    }


def load_records(payload: object) -> list[dict[str, object]]:
    if isinstance(payload, dict):
        records = payload.get("records")
        if records is None and isinstance(payload.get("result"), dict):
            records = payload["result"].get("records")  # type: ignore[index]
    elif isinstance(payload, list):
        records = payload
    else:
        records = None
    if not isinstance(records, list):
        raise ComputationalObjectiveError("Record JSON must be a list or an object with a records list.")
    cleaned: list[dict[str, object]] = []
    for raw in records:
        if not isinstance(raw, dict):
            continue
        pmid = str(raw.get("pmid") or raw.get("uid") or "").strip()
        if not pmid:
            continue
        record = dict(raw)
        record["pmid"] = pmid
        cleaned.append(record)
    if not cleaned:
        raise ComputationalObjectiveError("No records with PMIDs were found.")
    return cleaned


def text_terms(record: dict[str, object]) -> set[str]:
    text = f"{record.get('title') or ''} {record.get('abstract') or ''}"
    terms: set[str] = set()
    for match in TOKEN_RE.finditer(text):
        term = match.group(0).lower().strip("-")
        if len(term) < 3 or term in STOPWORDS:
            continue
        if term.isdigit():
            continue
        terms.add(term)
    return terms


def mesh_terms(record: dict[str, object]) -> set[str]:
    output: set[str] = set()
    for item in record.get("mesh_headings") or []:
        if isinstance(item, dict):
            name = item.get("name") or item.get("term") or item.get("label")
        else:
            name = item
        if name:
            output.add(" ".join(str(name).split()))
    return output


def record_relevance(record: dict[str, object]) -> bool:
    value = record.get("relevant")
    if value is None:
        value = record.get("is_relevant")
    if value is None:
        return True
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "relevant", "included"}


def parse_float_list(value: str | None, defaults: list[float]) -> list[float]:
    if not value:
        return defaults[:]
    result = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not result:
        raise ComputationalObjectiveError("Threshold list cannot be empty.")
    return result


def parse_int_list(value: str | None, defaults: list[int]) -> list[int]:
    if not value:
        return defaults[:]
    result = [int(part.strip()) for part in value.split(",") if part.strip()]
    if not result:
        raise ComputationalObjectiveError("MeSH count list cannot be empty.")
    return result


def load_split(payload: object) -> dict[str, list[str]]:
    if not isinstance(payload, dict):
        raise ComputationalObjectiveError("Split JSON must be an object.")
    return {
        "development": [str(pmid) for pmid in payload.get("development_pmids", [])],
        "validation": [str(pmid) for pmid in payload.get("validation_pmids", [])],
        "unseen": [str(pmid) for pmid in payload.get("unseen_pmids", [])],
    }


def load_background(term_rank_payload: object | None) -> tuple[dict[tuple[str, str], dict[str, object]], int | None]:
    if not isinstance(term_rank_payload, dict):
        return {}, None
    pubmed_total = term_rank_payload.get("pubmed_total_estimate")
    total = int(pubmed_total) if pubmed_total else None
    rows = term_rank_payload.get("ranked_terms")
    background: dict[tuple[str, str], dict[str, object]] = {}
    if not isinstance(rows, list):
        return background, total
    for row in rows:
        if not isinstance(row, dict):
            continue
        field = str(row.get("field") or "").strip().lower()
        term = normalize_term(str(row.get("term") or ""))
        if not field or not term:
            continue
        count = row.get("background_count")
        prevalence = None
        if count is not None and total:
            prevalence = int(count) / total
        background[(field, term)] = {
            "background_count": count,
            "background_prevalence": prevalence,
            "lift": row.get("lift"),
            "suggested_layer": row.get("suggested_layer"),
        }
    return background, total


def normalize_term(term: str) -> str:
    return " ".join(term.lower().split())


def load_category_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    payload = read_json(path)
    mapping: dict[str, str] = {}
    if not isinstance(payload, dict):
        raise ComputationalObjectiveError("Category map must be a JSON object.")
    for key, value in payload.items():
        if isinstance(value, str):
            category = normalize_category(value)
            if category:
                mapping[normalize_term(str(key))] = category
        elif isinstance(value, list):
            category = normalize_category(str(key))
            if not category:
                continue
            for term in value:
                mapping[normalize_term(str(term))] = category
    return mapping


def normalize_category(value: str) -> str | None:
    normalized = value.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "health_condition": "condition",
        "disease": "condition",
        "intervention": "treatment",
        "test": "treatment",
        "diagnostic_test": "treatment",
        "therapy": "treatment",
        "study": "study_type",
        "study_design": "study_type",
        "methodology": "study_type",
        "none": "none",
        "reject": "none",
        "rejected": "none",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized in {*VALID_CATEGORIES, "none"}:
        return normalized
    return None


def heuristic_category(term: str, field: str, mapping: dict[str, str]) -> str | None:
    normalized = normalize_term(term)
    if normalized in mapping:
        return mapping[normalized]
    words = set(re.split(r"[\s,/;-]+", normalized))
    if words & STUDY_TYPE_HINTS:
        return "study_type"
    if words & TREATMENT_HINTS:
        return "treatment"
    if words & CONDITION_HINTS:
        return "condition"
    if field == "mesh":
        if normalized.endswith(" neoplasms") or " diseases" in normalized:
            return "condition"
        if "therapeutics" in normalized or "diagnosis" in normalized:
            return "treatment"
    return None


def build_term_statistics(
    records: list[dict[str, object]],
    split: dict[str, list[str]],
    *,
    category_map: dict[str, str],
    term_rank_payload: object | None = None,
    require_background: bool = False,
) -> dict[str, object]:
    development_pmids = set(split["development"])
    if not development_pmids:
        raise ComputationalObjectiveError("Split must contain development PMIDs.")
    records_by_pmid = {str(record["pmid"]): record for record in records}
    missing = sorted(pmid for pmid in development_pmids if pmid not in records_by_pmid)
    dev_records = [records_by_pmid[pmid] for pmid in split["development"] if pmid in records_by_pmid]
    if not dev_records:
        raise ComputationalObjectiveError("No development records were found in records JSON.")

    background, pubmed_total = load_background(term_rank_payload)
    vector_map: dict[tuple[str, str], set[str]] = defaultdict(set)
    for record in dev_records:
        pmid = str(record["pmid"])
        for term in text_terms(record):
            vector_map[("tiab", term)].add(pmid)
        for term in mesh_terms(record):
            vector_map[("mesh", normalize_term(term))].add(pmid)

    terms: list[dict[str, object]] = []
    for (field, normalized), pmids in sorted(vector_map.items()):
        bg = background.get((field, normalized), {})
        if require_background and field == "tiab" and "background_prevalence" not in bg:
            continue
        category = heuristic_category(normalized, field, category_map)
        background_prevalence = bg.get("background_prevalence")
        terms.append(
            {
                "term": normalized,
                "field": field,
                "category": category or "none",
                "development_df": len(pmids),
                "development_prevalence": round(len(pmids) / len(dev_records), 6),
                "development_pmids": sorted(pmids),
                "background_count": bg.get("background_count"),
                "background_prevalence": background_prevalence,
                "lift": bg.get("lift"),
                "suggested_layer": bg.get("suggested_layer") or pubmed_atom(normalized, field),
            }
        )

    return {
        "operation": "term-stats",
        "development_count": len(dev_records),
        "missing_development_pmids": missing,
        "pubmed_total_estimate": pubmed_total,
        "require_background": require_background,
        "term_count": len(terms),
        "terms": terms,
        "warnings": []
        if pubmed_total
        else ["Population/background prevalence is unavailable for terms without term-rank evidence."],
    }


def pubmed_atom(term: str, field: str) -> str:
    if field == "mesh":
        escaped = " ".join(term.split())
        return f'"{escaped}"[Mesh:noexp]'
    escaped = " ".join(term.split())
    if " " in escaped:
        return f'"{escaped}"[tiab]'
    return f"{escaped}[tiab]"


def candidate_terms_for_grid(
    term_stats: dict[str, object],
    *,
    development_threshold: float,
    population_threshold: float,
    mesh_count: int,
) -> list[dict[str, object]]:
    raw_terms = [term for term in term_stats.get("terms", []) if isinstance(term, dict)]
    tiab_terms: list[dict[str, object]] = []
    mesh_terms_list: list[dict[str, object]] = []
    for term in raw_terms:
        category = term.get("category")
        if category not in VALID_CATEGORIES:
            continue
        field = term.get("field")
        dev_prev = float(term.get("development_prevalence") or 0.0)
        bg_prev_raw = term.get("background_prevalence")
        bg_prev = float(bg_prev_raw) if bg_prev_raw is not None else 0.0
        if field == "tiab":
            if dev_prev >= development_threshold and bg_prev <= population_threshold:
                tiab_terms.append(term)
        elif field == "mesh":
            mesh_terms_list.append(term)
    mesh_terms_list.sort(
        key=lambda row: (
            -int(row.get("development_df") or 0),
            str(row.get("term") or ""),
        )
    )
    return tiab_terms + mesh_terms_list[: max(0, mesh_count)]


def term_vector(term: dict[str, object], pmid_order: list[str]) -> tuple[int, ...]:
    pmids = set(str(pmid) for pmid in term.get("development_pmids", []))
    return tuple(1 if pmid in pmids else 0 for pmid in pmid_order)


def or_vectors(vectors: list[tuple[int, ...]], width: int) -> tuple[int, ...]:
    if not vectors:
        return tuple(1 for _ in range(width))
    return tuple(1 if any(vector[index] for vector in vectors) else 0 for index in range(width))


def and_vectors(vectors: list[tuple[int, ...]], width: int) -> tuple[int, ...]:
    if not vectors:
        return tuple(1 for _ in range(width))
    return tuple(1 if all(vector[index] for vector in vectors) else 0 for index in range(width))


def query_vector(terms: list[dict[str, object]], pmid_order: list[str]) -> tuple[int, ...]:
    by_category: dict[str, list[tuple[int, ...]]] = {category: [] for category in VALID_CATEGORIES}
    for term in terms:
        category = str(term.get("category") or "")
        if category in by_category:
            by_category[category].append(term_vector(term, pmid_order))
    category_vectors = [or_vectors(vectors, len(pmid_order)) for vectors in by_category.values() if vectors]
    return and_vectors(category_vectors, len(pmid_order))


def prune_terms(terms: list[dict[str, object]], development_pmids: list[str]) -> dict[str, object]:
    kept = list(terms)
    removed: list[dict[str, object]] = []
    baseline = query_vector(kept, development_pmids)
    baseline_count = sum(baseline)
    ordered = sorted(
        kept,
        key=lambda row: (
            -int(row.get("development_df") or 0),
            str(row.get("category") or ""),
            str(row.get("term") or ""),
        ),
    )
    for term in ordered:
        if term not in kept:
            continue
        category = str(term.get("category") or "")
        if category in VALID_CATEGORIES:
            category_remaining = [row for row in kept if row.get("category") == category]
            if len(category_remaining) <= 1:
                continue
        trial = [row for row in kept if row is not term]
        trial_vector = query_vector(trial, development_pmids)
        if sum(trial_vector) >= baseline_count:
            kept = trial
            removed.append(
                {
                    "term": term.get("term"),
                    "field": term.get("field"),
                    "category": term.get("category"),
                    "reason": "removal did not reduce development-set retrieval count",
                }
            )
            baseline = trial_vector
            baseline_count = sum(baseline)
    return {
        "kept_terms": kept,
        "removed_terms": removed,
        "development_retrieved_count": baseline_count,
        "development_total": len(development_pmids),
    }


def build_query(terms: list[dict[str, object]]) -> str:
    lines: list[str] = []
    for category in VALID_CATEGORIES:
        atoms = [
            str(term.get("suggested_layer") or pubmed_atom(str(term.get("term") or ""), str(term.get("field") or "tiab")))
            for term in terms
            if term.get("category") == category
        ]
        atoms = sorted(dict.fromkeys(atom for atom in atoms if atom))
        if not atoms:
            continue
        if len(atoms) == 1:
            block = atoms[0]
        else:
            block = "(\n  " + "\n  OR ".join(atoms) + "\n)"
        lines.append(block)
    return "\nAND\n".join(lines)


def beta_score(precision: float | None, recall: float, beta: float) -> float | None:
    if precision is None:
        return None
    if precision == 0 and recall == 0:
        return 0.0
    beta_sq = beta * beta
    return (1 + beta_sq) * precision * recall / ((beta_sq * precision) + recall)


def evaluate_terms_on_records(
    terms: list[dict[str, object]],
    records: list[dict[str, object]],
    pmids: list[str],
) -> dict[str, object]:
    pmid_set = set(pmids)
    universe = [record for record in records if str(record["pmid"]) in pmid_set]
    relevant = [record for record in universe if record_relevance(record)]
    has_nonrelevant = any(not record_relevance(record) for record in universe)
    retrieved_pmids: list[str] = []
    for record in universe:
        category_matches: list[bool] = []
        record_text_terms = text_terms(record)
        record_mesh_terms = {normalize_term(term) for term in mesh_terms(record)}
        for category in VALID_CATEGORIES:
            category_terms = [term for term in terms if term.get("category") == category]
            if not category_terms:
                continue
            matched = False
            for term in category_terms:
                value = str(term.get("term") or "")
                if term.get("field") == "mesh":
                    matched = value in record_mesh_terms
                else:
                    matched = value in record_text_terms
                if matched:
                    break
            category_matches.append(matched)
        if category_matches and all(category_matches):
            retrieved_pmids.append(str(record["pmid"]))
    relevant_pmids = {str(record["pmid"]) for record in relevant}
    retrieved_set = set(retrieved_pmids)
    tp = len(retrieved_set & relevant_pmids)
    recall = tp / len(relevant_pmids) if relevant_pmids else 0.0
    precision = tp / len(retrieved_set) if retrieved_set and has_nonrelevant else None
    return {
        "pmid_count": len(pmids),
        "record_count": len(universe),
        "relevant_count": len(relevant_pmids),
        "retrieved_count": len(retrieved_set),
        "true_positive_count": tp,
        "recall": round(recall, 6),
        "precision": round(precision, 6) if precision is not None else None,
        "f1": round(beta_score(precision, recall, 1.0), 6) if precision is not None else None,
        "f3": round(beta_score(precision, recall, 3.0), 6) if precision is not None else None,
        "retrieved_pmids": sorted(retrieved_set),
        "missed_relevant_pmids": sorted(relevant_pmids - retrieved_set),
        "precision_available": precision is not None,
    }


def objective_value(metrics: dict[str, object], objective: str) -> float:
    value = metrics.get(objective)
    if isinstance(value, (int, float)):
        return float(value)
    return float(metrics.get("recall") or 0.0)


def derive_query(
    records: list[dict[str, object]],
    split: dict[str, list[str]],
    *,
    category_map: dict[str, str],
    term_rank_payload: object | None,
    require_background: bool,
    development_thresholds: list[float],
    population_thresholds: list[float],
    mesh_counts: list[int],
    objective: str,
) -> dict[str, object]:
    term_stats = build_term_statistics(
        records,
        split,
        category_map=category_map,
        term_rank_payload=term_rank_payload,
        require_background=require_background,
    )
    grid_results: list[dict[str, object]] = []
    for dev_threshold in development_thresholds:
        for pop_threshold in population_thresholds:
            for mesh_count in mesh_counts:
                candidates = candidate_terms_for_grid(
                    term_stats,
                    development_threshold=dev_threshold,
                    population_threshold=pop_threshold,
                    mesh_count=mesh_count,
                )
                pruned = prune_terms(candidates, split["development"])
                kept_terms = pruned["kept_terms"]
                validation_metrics = evaluate_terms_on_records(kept_terms, records, split["validation"])
                grid_results.append(
                    {
                        "development_threshold": dev_threshold,
                        "population_threshold": pop_threshold,
                        "mesh_count": mesh_count,
                        "candidate_count": len(candidates),
                        "kept_count": len(kept_terms),
                        "removed_count": len(pruned["removed_terms"]),
                        "development_retrieved_count": pruned["development_retrieved_count"],
                        "validation": validation_metrics,
                        "score": round(objective_value(validation_metrics, objective), 6),
                    }
                )
    grid_results.sort(
        key=lambda row: (
            -float(row["score"]),
            int(row["kept_count"]),
            float(row["development_threshold"]),
            -float(row["population_threshold"]),
            int(row["mesh_count"]),
        )
    )
    if not grid_results:
        raise ComputationalObjectiveError("No grid results were generated.")
    best = grid_results[0]
    candidates = candidate_terms_for_grid(
        term_stats,
        development_threshold=float(best["development_threshold"]),
        population_threshold=float(best["population_threshold"]),
        mesh_count=int(best["mesh_count"]),
    )
    pruned = prune_terms(candidates, split["development"])
    kept_terms = pruned["kept_terms"]
    return {
        "operation": "derive-query",
        "method": "Computational objective PubMed query generation after Scells et al. 2020",
        "objective": objective,
        "objective_note": (
            "When precision is unavailable because only relevant records are supplied, F-score objectives fall back to recall."
        ),
        "split_summary": {
            "development_count": len(split["development"]),
            "validation_count": len(split["validation"]),
            "unseen_count": len(split["unseen"]),
        },
        "term_stats": {
            "development_count": term_stats["development_count"],
            "term_count": term_stats["term_count"],
            "warnings": term_stats["warnings"],
        },
        "selected_parameters": {
            "development_threshold": best["development_threshold"],
            "population_threshold": best["population_threshold"],
            "mesh_count": best["mesh_count"],
            "score": best["score"],
        },
        "selected_validation": best["validation"],
        "unseen_evaluation": evaluate_terms_on_records(kept_terms, records, split["unseen"]) if split["unseen"] else None,
        "candidate_count": len(candidates),
        "kept_count": len(kept_terms),
        "removed_count": len(pruned["removed_terms"]),
        "kept_terms": [
            {
                "term": term.get("term"),
                "field": term.get("field"),
                "category": term.get("category"),
                "development_df": term.get("development_df"),
                "development_prevalence": term.get("development_prevalence"),
                "background_prevalence": term.get("background_prevalence"),
                "suggested_layer": term.get("suggested_layer"),
            }
            for term in kept_terms
        ],
        "removed_terms": pruned["removed_terms"],
        "query": build_query(kept_terms),
        "grid_result_count": len(grid_results),
        "top_grid_results": grid_results[:20],
        "limitations": [
            "Default categorization is heuristic unless a manual category map or UMLS-derived map is supplied.",
            "The script uses unigrams and MeSH headings; phrase-level MetaMap extraction is not built in.",
            "Precision, NNR, F1, and F3 require labelled nonrelevant records; with relevant-only inputs, validation is recall-only.",
        ],
    }


def pmids_from_args(values: list[str] | None, file_path: str | None) -> list[str]:
    tokens: list[str] = []
    if values:
        tokens.extend(values)
    if file_path:
        try:
            tokens.append(Path(file_path).read_text(encoding="utf-8-sig"))
        except OSError as exc:
            raise ComputationalObjectiveError(f"Could not read PMID file: {file_path}") from exc
    return tokens


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Computational objective PubMed query generator.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    split_parser = subparsers.add_parser("split-references", help="Create 50/25/25 development, validation, and unseen PMID sets.")
    split_parser.add_argument("--pmids", nargs="*", help="PMIDs as separate values or comma/space/semicolon-separated strings.")
    split_parser.add_argument("--pmids-file", help="UTF-8 file containing PMIDs.")
    split_parser.add_argument("--seed", default=DEFAULT_SEED)
    split_parser.add_argument("--output")

    term_parser = subparsers.add_parser("term-stats", help="Compute development-set term vectors and statistics.")
    term_parser.add_argument("--records-json", required=True)
    term_parser.add_argument("--split-json", required=True)
    term_parser.add_argument("--term-rank-json")
    term_parser.add_argument("--category-map-json")
    term_parser.add_argument("--require-background", action="store_true")
    term_parser.add_argument("--output")

    derive_parser = subparsers.add_parser("derive-query", help="Grid-search, prune, and emit a computationally derived PubMed query.")
    derive_parser.add_argument("--records-json", required=True)
    derive_parser.add_argument("--split-json", required=True)
    derive_parser.add_argument("--term-rank-json")
    derive_parser.add_argument("--category-map-json")
    derive_parser.add_argument("--require-background", action="store_true")
    derive_parser.add_argument("--development-thresholds", help="Comma-separated thresholds, default 0.05..0.30.")
    derive_parser.add_argument("--population-thresholds", help="Comma-separated thresholds, default Scells-style grid.")
    derive_parser.add_argument("--mesh-counts", help="Comma-separated MeSH counts, default 0,1,5,10,15,20,25.")
    derive_parser.add_argument("--objective", choices=["recall", "precision", "f1", "f3"], default="recall")
    derive_parser.add_argument("--output")

    prune_parser = subparsers.add_parser("prune-query", help="Prune candidate terms from a saved term-stats JSON and explicit thresholds.")
    prune_parser.add_argument("--term-stats-json", required=True)
    prune_parser.add_argument("--split-json", required=True)
    prune_parser.add_argument("--development-threshold", type=float, default=0.20)
    prune_parser.add_argument("--population-threshold", type=float, default=0.02)
    prune_parser.add_argument("--mesh-count", type=int, default=20)
    prune_parser.add_argument("--output")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "split-references":
            result = split_pmids(pmids_from_args(args.pmids, args.pmids_file), seed=args.seed)
        elif args.command == "term-stats":
            records = load_records(read_json(args.records_json))
            split = load_split(read_json(args.split_json))
            term_rank_payload = read_json(args.term_rank_json) if args.term_rank_json else None
            result = build_term_statistics(
                records,
                split,
                category_map=load_category_map(args.category_map_json),
                term_rank_payload=term_rank_payload,
                require_background=args.require_background,
            )
        elif args.command == "derive-query":
            records = load_records(read_json(args.records_json))
            split = load_split(read_json(args.split_json))
            term_rank_payload = read_json(args.term_rank_json) if args.term_rank_json else None
            result = derive_query(
                records,
                split,
                category_map=load_category_map(args.category_map_json),
                term_rank_payload=term_rank_payload,
                require_background=args.require_background,
                development_thresholds=parse_float_list(args.development_thresholds, DEFAULT_DEVELOPMENT_THRESHOLDS),
                population_thresholds=parse_float_list(args.population_thresholds, DEFAULT_POPULATION_THRESHOLDS),
                mesh_counts=parse_int_list(args.mesh_counts, DEFAULT_MESH_COUNTS),
                objective=args.objective,
            )
        elif args.command == "prune-query":
            term_stats = read_json(args.term_stats_json)
            split = load_split(read_json(args.split_json))
            if not isinstance(term_stats, dict):
                raise ComputationalObjectiveError("Term-stats JSON must be an object.")
            candidates = candidate_terms_for_grid(
                term_stats,
                development_threshold=args.development_threshold,
                population_threshold=args.population_threshold,
                mesh_count=args.mesh_count,
            )
            pruned = prune_terms(candidates, split["development"])
            result = {
                "operation": "prune-query",
                "parameters": {
                    "development_threshold": args.development_threshold,
                    "population_threshold": args.population_threshold,
                    "mesh_count": args.mesh_count,
                },
                **pruned,
                "query": build_query(pruned["kept_terms"]),
            }
        else:
            parser.error(f"Unknown command: {args.command}")
    except ComputationalObjectiveError as exc:
        parser.exit(2, f"computational_objective.py: error: {exc}\n")

    write_json(result, getattr(args, "output", None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
