# Age-Related Patterns of DNA Methylation Changes

Analysis code accompanying the manuscript [Age-Related Patterns of DNA
Methylation Changes](https://www.biorxiv.org/content/10.1101/2024.12.10.627727v2).

## Study overview

This project analyzes DNA methylation measurements from 4,641 samples across
23 GEO datasets. It compares CpGs used by nine epigenetic clocks, identifies
age-associated CpGs, examines sex-specific differences, and clusters
non-linear methylation trajectories.

## Repository layout

```text
.
|-- README.md
|-- ClusteringHelpers.Rmd     # functions and dependencies for clustering
|-- Clustering.Rmd            # clustering analysis and Figures 2 and 5
|-- scripts/
|   |-- PreprocHelpers.Rmd    # preprocessing functions and dependencies
|   |-- Preproc.Rmd           # preprocessing and supplementary figures
|   |-- FilterHelpers.Rmd     # filtering functions and dependencies
|   `-- Filtering.Rmd         # filtering, tables, and manuscript figures
|-- data/
|   `-- README.md             # expected input-data layout
|-- fig_repository/           # final manuscript figures
`-- extra/                    # archived exploratory and superseded analyses
```

Large datasets and generated intermediate files are intentionally not included
in the repository. See `data/README.md` for the expected local layout.

## Pipeline order

Run notebooks in this order:

1. `scripts/PreprocHelpers.Rmd`, then `scripts/Preproc.Rmd`
2. `scripts/FilterHelpers.Rmd`, then `scripts/Filtering.Rmd`
3. `ClusteringHelpers.Rmd`, then `Clustering.Rmd`

The helper notebook for each stage defines the functions, packages, and shared
objects used by the analysis notebook. Run every chunk in the helper before
running the corresponding analysis notebook. The analysis is deliberately
split into chunks and writes checkpoints so a long run can be resumed.

## Working directories

The preprocessing and filtering notebooks use paths relative to `scripts/`.
Open those notebooks from that directory. The clustering notebooks use paths
relative to the repository root. Changing the working directory will cause
input-file errors.

## Local folders

Create these untracked folders before running the analysis:

```text
data/Raw/                 one folder per GEO accession
data/pp_datasets/         standardized datasets
data/matrix/              age-by-CpG matrices
data/ref/                 manifests and clock reference files
data/splits/              temporary matrix splits
intermediates/            saved R objects and checkpoints
relevant_rds/             filtering results and downstream R objects
fig_repository/           final manuscript figures
```

Every active notebook loads `project_paths.R` first and resolves paths from the
repository root. The code therefore uses the same locations whether a notebook
is opened from the root or from `scripts/`.

Each `data/Raw/GSE<accession>/` folder is expected to contain
`methylation_data.csv` and `metadata.csv`. A `cleaned_metadata.csv` file may
also be present, but preprocessing deliberately starts from `metadata.csv`.

## Reproducibility notes

- Raw GEO and reference data are not bundled in this archive.
- The supplied `relevant_rds` folder contains the three filtering-result objects.
  Later filtering chunks generate additional RDS objects needed by clustering.
- Clock definitions and the HM450 manifest under `data/ref/` were not included
  in the supplied downloads and must be provided separately.
- The notebooks require R plus CRAN and Bioconductor packages loaded in the
  helper notebooks.
- Some steps are computationally intensive and were designed to resume from
  saved `.rds` files.
- `extra/` is retained for provenance, but it is not part of the current
  three-stage pipeline.

## Manuscript outputs

- Preprocessing: Supplementary Figures 1, 2, and 5
- Filtering: Supplementary Figures 3 and 4; Table 1; Figures 3 and 4
- Clustering: Figures 2 and 5

## Contact

Questions about the scientific analysis: [kchen24@stanford.edu](mailto:kchen24@stanford.edu)

Repository documentation reorganized July 2026. Scientific analysis code and
reported results were not altered by this documentation cleanup.
