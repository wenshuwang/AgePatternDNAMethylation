# Age-Related Patterns of DNA Methylation Changes

Analysis code for the manuscript [Age-Related Patterns of DNA Methylation
Changes](https://www.biorxiv.org/content/10.1101/2024.12.10.627727v2).

The study analyzes 4,641 samples from 23 GEO datasets to characterize
age-associated CpGs, sex differences, epigenetic clocks, and non-linear DNA
methylation trajectories.

## Start here

The live pipeline is in [`analysis/`](analysis/README.md). Files are numbered in
the order they should be run:

1. preprocessing
2. filtering
3. clustering

Each stage has a `_helpers.Rmd` notebook that must be run before its matching
analysis notebook. For example, run `01_preprocessing_helpers.Rmd` before
`01_preprocessing.Rmd`.

For a complete reproduction directly from the 23 raw datasets, the simplest
entry point is:

```r
source("run_preprocessing.R")
```

That command rebuilds the preprocessing checkpoints, 10 horizontal split
files, and all six mean/SD matrices. It does not read anything from
`comparison/kevin_old/`.

## Repository map

```text
.
|-- analysis/             numbered notebooks for the current pipeline
|-- R/                    shared project-path utilities
|-- data/                 local raw, processed, matrix, and reference data
|-- comparison/           optional archived-original analysis inputs
|-- relevant_rds/         generated R checkpoints (not tracked by Git)
|-- fig_repository/       final manuscript figures
|-- archive/              exploratory and superseded analyses
|-- run_preprocessing.R   exact raw-data preprocessing entry point
|-- check_inputs.R        reports missing required inputs
`-- REORGANIZATION_NOTES.md
```

Large data, RDS checkpoints, manuscript drafts, and private working files are
not stored in Git.

## Expected data layout

```text
data/
|-- Raw/
|   `-- GSE<accession>/
|       |-- methylation_data.csv
|       |-- metadata.csv
|       `-- cleaned_metadata.csv     optional
|-- pp_datasets/
|   `-- GSE<accession>_pp.csv
|-- matrix/
|   |-- allTRUE.csv
|   |-- allFALSE.csv
|   |-- maleTRUE.csv
|   |-- maleFALSE.csv
|   |-- femaleTRUE.csv
|   `-- femaleFALSE.csv
`-- ref/                            clock definitions and HM450 manifest
```

The notebooks locate the repository root through `R/project_paths.R`, so their
data paths do not depend on the directory from which they are opened.

## Rebuilding from raw data

`run_preprocessing.R` uses the raw methylation tables and metadata to derive
the eligible 0-80-year-old samples and CpG coverage counts from scratch. It
stops before the expensive split/matrix work unless the raw inputs reproduce
all original checkpoints:

- 486,427 rows in `running_site_list.rds` (including 850 liftover artifacts
  with missing counts)
- 485,577 standard HM450 probes in `site_list_by_samples.rds`
- 393,628 CpGs present in at least 18 of 23 datasets
- 256,529 CpGs measured in at least 4,176 valid samples
- 4,641 valid samples

The 393,628-site manuscript checkpoint is saved as
`intermediates/POST_sites_18of23_datasets.rds`. The downstream matrix
background is the distinct 256,529-site list saved as
`intermediates/POST_sites.rds`.

Reproducing these exact counts still requires the same raw tables, metadata
edits, and HM450 manifest used for the original analysis. A count mismatch is
reported explicitly rather than being hidden or replaced with Kevin's saved
list.

## Starting from archived original results

Raw data are only required to rerun preprocessing. Filtering can instead
resume from the original saved checkpoints and matrices. Put the archived
files in `comparison/kevin_old/` using the layout documented in
[`comparison/README.md`](comparison/README.md).

Filtering selects inputs as follows:

1. If locally generated matrices and both validated `POST_sites` lists exist,
   use them after checking their counts and matrix row identities.
2. Otherwise, if the complete archived-original set exists, validate and use
   it.
3. Stop with a clear missing-file report if neither set is complete.

The archived validation requires the original counts of 393,628 CpGs present
in at least 18 of 23 datasets, 256,529 downstream CpGs, and 4,641 valid
samples. This prevents an original CpG list from being silently combined with
matrices produced from a different preprocessing run.

When archived mode starts, filtering writes these two validated lists to
`intermediates/POST_sites_18of23_datasets.rds` and
`intermediates/POST_sites_4176_samples.rds`. It does not overwrite a locally
rebuilt `intermediates/POST_sites.rds`.

To override automatic selection for the current R session:

```r
Sys.setenv(AGEPATTERN_INPUT_MODE = "archived") # require original saved inputs
Sys.setenv(AGEPATTERN_INPUT_MODE = "local")    # require locally generated inputs
```

## Check inputs

From the repository root, run:

```r
source("check_inputs.R")
```

The checker reports preprocessing, filtering, and clustering readiness
separately. Missing raw files do not block filtering when either a complete
local matrix set or the archived-original input set is available. It does not
download or modify anything.

## Main outputs

- Preprocessing: Supplementary Figures 1, 2, and 5
- Filtering: Supplementary Figures 3 and 4; Table 1; Figures 3 and 4
- Clustering: Figures 2 and 5

## Notes

- `archive/` is retained for provenance and is not part of the live pipeline.
- Additional RDS objects needed by clustering are created during filtering.
- Reference clock files and the HM450 manifest must be available under
  `data/ref/`.

## Contact

[kchen24@stanford.edu](mailto:kchen24@stanford.edu)
