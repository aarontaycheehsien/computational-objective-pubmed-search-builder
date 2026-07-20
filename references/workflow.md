# Workflow

This skill produces an algorithmic PubMed candidate query. It is a starting point for expert review, not a final search strategy.

## 1. Intake

Require a review title/topic and known relevant PMIDs. If the user supplies only a topic, ask for known relevant PMIDs or a saved record JSON.

Do not use pasted Boolean syntax as term evidence. Existing Boolean can be recorded as a comparator only.

## 2. Fetch Or Load Records

If the user supplies PMIDs, fetch and mine them:

```bash
python scripts/pubmed_tool.py mine --pmids <PMIDs> --output known_relevant_mine.json
```

Inspect the saved JSON before using it. Exclude missing, malformed, retracted, or clearly out-of-scope records.

If the user already has a saved record JSON, it must be a list of records or an object with `records`. Each record needs `pmid`, `title`, and ideally `abstract`, `mesh_headings`, and optional `relevant`.

## 3. Split References

Create development, validation, and unseen sets:

```bash
python scripts/computational_objective.py split-references --pmids-file known_relevant_pmids.txt --output split.json
```

Default split is 50% development, 25% validation, and 25% unseen. The script warns when the development set is below the paper's minimum of 25 records.

## 4. Collect Background Evidence

Run objective background ranking on the development PMIDs when network access is available:

```bash
python scripts/pubmed_tool.py term-rank --pmids <development PMIDs> --fields tiab,mesh --max-terms 500 --output term_rank_development.json
```

Use `--require-background` in the derive step when you want candidate text terms limited to terms with PubMed background evidence.

## 5. Categorize Terms

Provide a category map when possible. See `categorization.md`.

Without a map, the script uses a lightweight heuristic. That is useful for prototyping but not equivalent to the paper's MetaMap/UMLS semantic-type categorization.

## 6. Derive Candidate Query

Run:

```bash
python scripts/computational_objective.py derive-query --records-json known_relevant_mine.json --split-json split.json --term-rank-json term_rank_development.json --category-map-json category_map.json --require-background --objective recall --output derived_query.json
```

The output includes selected parameters, kept and removed terms, validation metrics, unseen metrics if available, and the generated PubMed query.

## 7. PubMed Count And QA

Save the generated query to `derived_query.txt`, then run:

```bash
python scripts/pubmed_tool.py search --query-file derived_query.txt --retmax 0 --output derived_query_search.json --summary
python scripts/hooks_tool.py final-qa --strategy-file derived_query.txt
```

The count is a workload proxy. It is not precision unless labelled nonrelevant records are available.

## 8. Handoff

Report:

- input relevant-record source
- split sizes
- categorization mode
- parameter grid and selected parameters
- validation and unseen known-set recall
- whether precision/F-score metrics were available
- generated query
- limitations and recommended expert-review steps
