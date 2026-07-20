---
name: computational-objective-pubmed
description: "Generate automated candidate PubMed Boolean queries from known relevant PMIDs using a Scells et al.-style computational adaptation of the Hausner objective method. Use when Codex needs to split references into development/validation/unseen sets, compute objective term statistics, categorize terms into condition/treatment/study-type clauses, run threshold grid search, apply bit-vector pruning, evaluate generated queries, or produce an algorithmic starting query for an information specialist. Do not use as a final librarian-reviewed search; label outputs as computational candidates."
---

# Computational Objective PubMed

## Core Goal

Create an automated PubMed Boolean candidate query from known relevant records. This skill adapts Scells et al. 2020, "A Computational Approach for Objectively Derived Systematic Review Search Strategies" (PMC7148214), which computationalizes Hausner-style objective search development.

Use this skill to produce a reproducible baseline query and evaluation report. Use `objective-pubmed` afterward when a conservative, audit-first, information-specialist workflow is needed.

## Start Here

Read `references/workflow.md` for every build.

Require:

- a plain-language review topic or title
- known relevant PMIDs or a saved PubMed record JSON containing known relevant records

If the user has no known relevant PMIDs, stop and ask for them or use another skill to build a non-computational search. This skill does not discover a gold standard by itself.

## Required References

- `references/workflow.md` - end-to-end build sequence.
- `references/computational-method.md` - Scells-style algorithm and v1 approximations.
- `references/categorization.md` - term category modes and mapping file format.
- `references/evaluation-and-limitations.md` - metrics, validation, and handoff caveats.
- `references/tools.md` - bundled script commands.

## Bundled Scripts

- `scripts/computational_objective.py` - split PMIDs, compute term statistics, grid-search thresholds, prune terms, assemble PubMed query, and evaluate known-set recall.
- `scripts/pubmed_tool.py` - PubMed E-utilities fetch, mine, term-rank, search, sample, validate, and audit scaffold.
- `scripts/mesh_tool.py` - MeSH lookup, tree context, and sweeps.
- `scripts/manifest_tool.py` - provenance ledger.
- `scripts/hooks_tool.py` - final query QA.
- `scripts/audit_markdown.py` - render audit Markdown when structured audit JSON exists.

## Non-Negotiable Rules

- Treat generated queries as computational candidates, not final searches.
- Do not claim paper-equivalent MetaMap/UMLS categorization unless a UMLS/MetaMap-derived mapping file was supplied.
- Do not report precision, NNR, F1, or F3 as meaningful unless labelled nonrelevant records are available. With relevant-only PMIDs, report known-set recall only.
- Keep the unseen set untouched until after validation-set parameter tuning.
- Preserve JSON outputs, final query text, and `run_manifest.json` for reproducibility.
