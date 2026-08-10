# Analysis pipeline

Run the notebooks in numerical order. Within each numbered stage, run the
helper notebook first because it loads packages and defines functions used by
the main notebook.

| Order | Helper notebook | Analysis notebook | Main outputs |
|---|---|---|---|
| 1 | `01_preprocessing_helpers.Rmd` | `01_preprocessing.Rmd` | standardized data, matrices, supplementary figures |
| 2 | `02_filtering_helpers.Rmd` | `02_filtering.Rmd` | age-associated CpGs, tables, filtering figures |
| 3 | `03_clustering_helpers.Rmd` | `03_clustering.Rmd` | trajectory clusters, Figures 2 and 5 |

## Resume behavior

The notebooks save intermediate `.rds` files so long analyses can resume from
later chunks. Run earlier chunks at least once to create the objects required
by downstream stages.

If the raw-data preprocessing stage cannot be rerun, filtering can resume from
the validated original checkpoint and matrix files under
`comparison/kevin_old/`. See the root `README.md` and `comparison/README.md`
for the required filenames. The filtering notebook automatically keeps this
archived set separate from locally rebuilt matrices.

## Paths

Every notebook loads `../R/project_paths.R` and then resolves inputs from the
repository root. Do not replace these with computer-specific absolute paths.

## Archived code

Exploratory and superseded notebooks are under `archive/`. They are kept for
history but are not prerequisites for this pipeline.
