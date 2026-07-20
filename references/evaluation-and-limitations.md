# Evaluation And Limitations

## Metrics

With known relevant PMIDs only, report:

- development known-set retrieval
- validation known-set recall
- unseen known-set recall, if an unseen set exists
- generated PubMed result count as workload proxy

Do not call the workload count precision.

When labelled nonrelevant records are supplied in `records` with `relevant: false`, the script can compute precision and F-scores over those labelled records. These are still collection-relative metrics, not absolute PubMed precision.

## Objective Selection

`derive-query` supports:

- `--objective recall`
- `--objective precision`
- `--objective f1`
- `--objective f3`

If precision is unavailable, F-score and precision objectives fall back to recall. State that fallback in the report.

## Validation Discipline

Use development records to derive terms and prune. Use validation to select thresholds. Use unseen only once for the final estimate.

If the query is changed after unseen evaluation, the unseen result is no longer an independent final estimate.

## Important Limitations

- Generated queries are candidate starting points.
- Heuristic categorization can place terms incorrectly.
- Unigram extraction misses phrase semantics such as "myocardial infarction".
- MeSH additions can increase recall while lowering precision.
- Known-set recall does not prove complete retrieval of all relevant PubMed records.
- PubMed query translation, indexing lag, and missing abstracts can affect results.

## Recommended Handoff Language

Use wording like:

```text
This is a computationally derived PubMed candidate query. It was generated from the supplied known relevant records using development/validation/unseen splitting, threshold tuning, and bit-vector pruning. It should be reviewed by an information specialist before use as a final systematic review search.
```
