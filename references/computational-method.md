# Computational Method

This skill follows the structure of Scells et al. 2020, with pragmatic v1 substitutions where required infrastructure is unavailable.

## Paper-Inspired Steps

1. Start from a test set of known relevant studies.
2. Split into development, validation, and unseen sets.
3. Identify candidate title/abstract terms and MeSH terms from the development set.
4. Filter terms by development-set document frequency and population/background document frequency.
5. Assign terms to query categories: health condition, treatment/intervention/test, and study type.
6. Assemble category `OR` clauses and combine categories with `AND`.
7. Use bit-vector pruning to remove terms that do not reduce development-set retrieval.
8. Tune thresholds on validation data.
9. Evaluate the selected query on unseen data.

## V1 Approximations

- Candidate text terms are unigrams from titles and abstracts.
- MetaMap is not bundled. Categorization is manual-map or heuristic unless the user supplies a UMLS/MetaMap-derived category map.
- Population document frequency comes from `pubmed_tool.py term-rank` when provided. Without that file, background evidence is marked incomplete.
- The pruning implementation keeps at least one term in each active category. This preserves interpretable concept blocks even when deleting the last category term would not reduce development-set recall.
- Precision, NNR, F1, and F3 require labelled nonrelevant records. With only known relevant PMIDs, the script reports recall-style known-item metrics.

## Threshold Grid

Default development thresholds:

```text
0.05, 0.10, 0.15, 0.20, 0.25, 0.30
```

Default population thresholds:

```text
0.001, 0.01, 0.02, 0.05, 0.10, 0.20
```

Default MeSH counts:

```text
0, 1, 5, 10, 15, 20, 25
```

Use custom grids only when a protocol or experiment requires them.

## Bit-Vector Pruning

The script represents each term by the development PMIDs it retrieves. Terms within a category are OR-ed. Active categories are AND-ed. Candidate terms are considered for removal in deterministic order by descending development document frequency.

A term is removed when removal does not reduce the number of development PMIDs retrieved and it is not the final term in its category.
