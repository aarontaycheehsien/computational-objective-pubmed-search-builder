# Categorization

The paper uses MetaMap and UMLS semantic types/groups to assign terms to Hausner categories. This skill does not bundle MetaMap.

## Modes

### Manual Map

Preferred for serious use. Provide `--category-map-json`.

Supported JSON forms:

```json
{
  "condition": ["cancer", "neoplasm", "tumor"],
  "treatment": ["brachytherapy", "radiotherapy"],
  "study_type": ["randomized", "cohort"]
}
```

or:

```json
{
  "cancer": "condition",
  "brachytherapy": "treatment",
  "randomized": "study_type"
}
```

Accepted categories:

- `condition`
- `treatment`
- `study_type`
- `none`

Aliases such as `health_condition`, `intervention`, `diagnostic_test`, and `study_design` are normalized.

### Heuristic

Default fallback. The script uses small keyword families to guess condition, treatment/intervention/test, and study type. Label all such output as heuristic and not paper-equivalent.

### UMLS/MetaMap-Derived

If the user supplies a mapping generated from MetaMap/UMLS, use it as a manual map. Then the output may be described as using UMLS-derived categorization, but still document the source, UMLS version, MetaMap settings, and any unmapped terms.

## Rejected Terms

Terms mapped to `none` or left uncategorized are not included in the generated query. Report how many terms were uncategorized when handing off the result.
