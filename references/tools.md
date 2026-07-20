# Tools

## Core Script

Split PMIDs:

```bash
python scripts/computational_objective.py split-references --pmids-file known_relevant_pmids.txt --output split.json
```

Compute term statistics:

```bash
python scripts/computational_objective.py term-stats --records-json known_relevant_mine.json --split-json split.json --term-rank-json term_rank_development.json --category-map-json category_map.json --output term_stats.json
```

Derive query:

```bash
python scripts/computational_objective.py derive-query --records-json known_relevant_mine.json --split-json split.json --term-rank-json term_rank_development.json --category-map-json category_map.json --require-background --objective recall --output derived_query.json
```

Prune fixed-threshold candidates:

```bash
python scripts/computational_objective.py prune-query --term-stats-json term_stats.json --split-json split.json --development-threshold 0.20 --population-threshold 0.02 --mesh-count 20 --output pruned_query.json
```

## PubMed Support

Fetch and mine PMIDs:

```bash
python scripts/pubmed_tool.py mine --pmids <PMIDs> --output known_relevant_mine.json
```

Compute PubMed background term evidence:

```bash
python scripts/pubmed_tool.py term-rank --pmids <development PMIDs> --fields tiab,mesh --max-terms 500 --output term_rank_development.json
```

Count generated query:

```bash
python scripts/pubmed_tool.py search --query-file derived_query.txt --retmax 0 --output derived_query_search.json --summary
```

Inspect samples when needed:

```bash
python scripts/pubmed_tool.py sample --query-file derived_query.txt --retmax 10 --output derived_query_sample.json
```

## QA And Manifest

```bash
python scripts/hooks_tool.py final-qa --strategy-file derived_query.txt
python scripts/manifest_tool.py init --manifest run_manifest.json --topic-slug topic
python scripts/manifest_tool.py add --manifest run_manifest.json --kind artifact --command "python scripts/computational_objective.py derive-query ..." --output derived_query.json --label "computational derived query"
python scripts/manifest_tool.py show --manifest run_manifest.json --validate --check-files
```
