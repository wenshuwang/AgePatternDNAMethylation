# Data layout

The study data are too large to distribute in this repository. The analysis
notebooks expect the following local directories:

```text
data/
|-- Raw/          # one folder per GEO accession
|-- original_datasets/ # original oGSE beta matrices for the 18-of-23 rebuild
|-- pp_datasets/  # standardized per-dataset CSV files
|-- splits/       # horizontal matrix chunks
|-- vert_splits/  # sample-oriented matrix chunks
|-- matrix/       # combined age-by-CpG mean and SD matrices
`-- ref/          # array manifests and epigenetic-clock CpG definitions
```

Do not commit participant-level or controlled-access data. The `.gitignore`
excludes generated and local data directories while preserving this file.

## Naming conventions

- Raw methylation values: `Raw/GSE<accession>/methylation_data.csv`
- Raw metadata: `Raw/GSE<accession>/metadata.csv`
- Optional cleaned metadata: `Raw/GSE<accession>/cleaned_metadata.csv`
- Processed datasets: `GSE<accession>_pp.csv`
- Combined matrices: `<group><average>.csv`, where group is `all`, `male`, or
  `female`, and average is `TRUE` for means or `FALSE` for standard deviations

Reference files used by the notebooks belong under `data/ref/`, including the
HM450 manifest and the `CpGsToInvestigate/` clock-definition files.

For access to the study's complete data bundle, contact the repository authors.
