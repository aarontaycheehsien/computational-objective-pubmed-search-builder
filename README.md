# Computational Objective PubMed Search Builder

A Codex skill and standalone Python toolkit for generating reproducible candidate PubMed Boolean searches from a set of known relevant studies.

The workflow is a practical computational adaptation of the objective search-development method described by Scells et al. (2020). It divides known relevant records into development, validation, and unseen sets; mines candidate title/abstract and MeSH terms; tunes term-selection thresholds; applies bit-vector pruning; and evaluates the resulting query against held-out known studies.

> [!IMPORTANT]
> The generated strategy is an algorithmic starting point, not a final librarian-reviewed search. Known-item recall does not establish that every relevant PubMed record has been retrieved.

## What this repository provides

- A Codex `SKILL.md` that guides the complete workflow.
- Deterministic development/validation/unseen splitting.
- Candidate-term statistics and threshold grid search.
- Condition, treatment/test, and study-type concept categorization.
- Bit-vector term pruning while retaining interpretable concept blocks.
- PubMed E-utilities helpers for record retrieval, background ranking, query counts, sampling, and validation.
- MeSH lookup and tree-context helpers.
- Query QA, provenance manifests, and Markdown audit-report generation.
- A standard-library-only Python implementation with an offline unit-test suite.

## Method overview

1. Supply a review topic and a set of known relevant PMIDs.
2. Fetch or load the corresponding PubMed records.
3. Split the records 50%/25%/25% into development, validation, and unseen sets.
4. Mine title/abstract terms and MeSH headings from the development set.
5. Compare development-set document frequency with PubMed background frequency.
6. Assign terms to condition, treatment/intervention/test, or study-type clauses.
7. Tune frequency thresholds and MeSH counts on the validation set.
8. Prune terms that do not improve development-set retrieval.
9. Evaluate the selected query once against the unseen set.
10. Run PubMed count and syntax QA, then hand the candidate to an information specialist.

The defaults follow the structure of the published method, with documented pragmatic substitutions. See [the computational method](references/computational-method.md) and [evaluation limitations](references/evaluation-and-limitations.md).

## Requirements

- Python 3.10 or later
- Internet access for commands that query PubMed or the MeSH RDF service
- Known relevant PMIDs, or a saved JSON file containing known relevant PubMed records

The core scripts and tests use only the Python standard library.

## Install as a Codex skill

Clone the repository into your personal Codex skills directory:

```powershell
git clone https://github.com/aarontaycheehsien/computational-objective-pubmed-search-builder.git "$HOME/.codex/skills/computational-objective-pubmed"
```

Restart Codex if the skill is not discovered immediately. Invoke it as `$computational-objective-pubmed` and provide a review topic plus known relevant PMIDs.

The scripts can also be used directly without installing the repository as a Codex skill.

## Quick start

Create a text file containing one known relevant PMID per line, then run:

```bash
# Fetch the known records. Replace the example PMIDs with your own.
python scripts/pubmed_tool.py mine \
  --pmids 12345678 23456789 34567890 \
  --output known_relevant_mine.json

# Create deterministic development, validation, and unseen sets.
python scripts/computational_objective.py split-references \
  --pmids-file known_relevant_pmids.txt \
  --output split.json

# Collect PubMed background-frequency evidence using development PMIDs.
python scripts/pubmed_tool.py term-rank \
  --pmids 12345678 23456789 \
  --fields tiab,mesh \
  --max-terms 500 \
  --output term_rank_development.json

# Derive and evaluate a candidate query.
python scripts/computational_objective.py derive-query \
  --records-json known_relevant_mine.json \
  --split-json split.json \
  --term-rank-json term_rank_development.json \
  --category-map-json category_map.json \
  --require-background \
  --objective recall \
  --output derived_query.json
```

`category_map.json` should assign candidate terms to `condition`, `treatment`, `study_type`, or `none`. Two supported formats and their aliases are documented in [categorization](references/categorization.md). If no map is supplied, the tool uses a lightweight heuristic that is not equivalent to the paper's MetaMap/UMLS categorization.

After extracting the generated query to `derived_query.txt`, count and check it:

```bash
python scripts/pubmed_tool.py search \
  --query-file derived_query.txt \
  --retmax 0 \
  --output derived_query_search.json \
  --summary

python scripts/hooks_tool.py final-qa --strategy-file derived_query.txt
```

See the [end-to-end workflow](references/workflow.md) and [command reference](references/tools.md) for more detail.

## Inputs and outputs

The principal inputs are:

- `known_relevant_pmids.txt`: one PMID per line.
- `known_relevant_mine.json`: PubMed records containing `pmid`, `title`, and preferably `abstract` and `mesh_headings`.
- `category_map.json`: an optional manual or UMLS/MetaMap-derived term-to-category mapping.
- `term_rank_development.json`: optional PubMed background-frequency evidence.

The derivation output records the candidate query, selected grid parameters, kept and removed terms, validation metrics, and unseen metrics when available. Preserve these files with `run_manifest.json` to maintain an audit trail.

## Metrics and interpretation

When only relevant records are supplied, the meaningful evaluation is known-set recall on development, validation, and unseen records. A PubMed result count is a workload proxy, not precision.

Precision, NNR, F1, and F3 require labelled nonrelevant records. Even then, reported precision is relative to the labelled collection rather than absolute PubMed precision.

## Important limitations

- The toolkit approximates Scells et al.; it is not a paper-equivalent reproduction.
- MetaMap is not bundled. Default categorization is heuristic unless you provide a manual or UMLS/MetaMap-derived mapping.
- Candidate free-text terms are unigrams, so phrase meaning can be lost.
- MeSH expansion may improve recall while substantially increasing screening workload.
- PubMed translation, indexing lag, incomplete abstracts, and the supplied known-study set affect results.
- Changing the query after inspecting the unseen set invalidates that set as an independent final evaluation.

## Testing

Run the offline test suite from the repository root:

```bash
python -m unittest discover -s tests -v
```

## Repository structure

```text
SKILL.md                         Codex skill instructions
agents/openai.yaml               Skill display metadata
scripts/computational_objective.py  Core query derivation
scripts/pubmed_tool.py           PubMed retrieval and evaluation helpers
scripts/mesh_tool.py             MeSH lookup helpers
scripts/hooks_tool.py            Query QA checks
scripts/manifest_tool.py         Reproducibility manifest
scripts/audit_markdown.py        Audit-report renderer
references/                      Workflow and methodological documentation
tests/                           Offline unit tests
```

## Methodological basis

Scells, H., Zuccon, G., Koopman, B., & Clark, J. (2020). *A Computational Approach for Objectively Derived Systematic Review Search Strategies*. In **Advances in Information Retrieval**, 385–398. [https://doi.org/10.1007/978-3-030-45439-5_26](https://doi.org/10.1007/978-3-030-45439-5_26). [Free full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7148214/).

That work computationally adapts the earlier Hausner objective method for systematic-review search development.

## License

[MIT](LICENSE) © 2024–2026 Aaron Tay.
