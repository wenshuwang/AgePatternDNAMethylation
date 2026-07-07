# Reorganization notes

This cleanup improves navigation without changing statistical methods or
scientific results.

## Changes made

- Rewrote the main README with the actual file locations, execution order,
  expected working directories, data layout, and manuscript outputs.
- Replaced the allow-list `.gitignore`, which hid nearly every new file, with
  explicit rules for local data, checkpoints, rendered notebooks, and R state.
- Fixed the missing closing quote in the `Clustering.Rmd` YAML title.
- Documented the expected data filenames and reference-data locations.
- Clearly marked `extra/` as an archive rather than an active pipeline stage.
- Centralized working-directory handling in `project_paths.R`.
- Updated raw-data links to `data/Raw/GSE<accession>/methylation_data.csv`
  and `metadata.csv`.
- Updated processed-data, matrix, and checkpoint links to `data/pp_datasets`,
  `data/matrix`, and `relevant_rds` respectively.

## Recommended future work

The current notebooks rely on session state created by manually running helper
notebooks. A deeper refactor should extract reusable functions into ordinary R
files, use a project-root path helper such as `here`, and capture package
versions with `renv`. That work should be validated against the original saved
intermediates and manuscript figures before replacing the published workflow.
