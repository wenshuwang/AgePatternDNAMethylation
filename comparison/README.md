# Original-run comparison inputs

Large archived checkpoints and matrices are not committed to Git. To resume
the exact original filtering run without reprocessing raw data, place the
files below under `comparison/kevin_old/`:

```text
running_site_list.rds
site_list_by_samples.rds
valid_samples.rds
POST_sites.rds
allTRUE_OG.csv
allFALSE_OG.csv
maleTRUE_OG.csv
maleFALSE_OG.csv
femaleTRUE_OG.csv
femaleFALSE_OG.csv
```

When this complete set is present and local matrices are absent,
`analysis/02_filtering.Rmd` uses it as a fallback and validates the original
checkpoints before analysis. Complete locally generated inputs are preferred.
The expected counts are 393,628 CpGs present in at least 18 of 23 datasets,
256,529 CpGs in the downstream sample-coverage list, and 4,641 valid samples.

Set `Sys.setenv(AGEPATTERN_INPUT_MODE = "local")` before running filtering to
force locally generated `data/matrix/` and `intermediates/POST_sites.rds`
instead. Set the value to `"archived"` to require the archived set.
