# Rebuild only the six manuscript age-summary matrices from existing validated
# splits. Use this recovery entry point when 256,529-row matrices were created
# by an earlier fork version; raw and PP preprocessing do not need to be rerun.

source("R/project_paths.R")
set_project_root()
source("R/preprocessing_core.R")

matrix_dir <- project_file("data", "matrix")
old_matrix_dir <- project_file("data", "matrix_256529")
matrix_names <- as.vector(outer(
  c("all", "male", "female"),
  c("TRUE.csv", "FALSE.csv"),
  paste0
))
current_matrices <- file.path(matrix_dir, matrix_names)

# Preserve the completed stricter-background run for comparison without
# copying several gigabytes. Renaming the folder does not duplicate its data.
if (all(file.exists(current_matrices))) {
  first_ids <- data.table::fread(
    current_matrices[[1L]], select = 1L, showProgress = FALSE
  )[[1L]]
  if (length(first_ids) == 256529L) {
    if (dir.exists(old_matrix_dir)) {
      stop("Archive folder already exists: ", old_matrix_dir)
    }
    if (!file.rename(matrix_dir, old_matrix_dir)) {
      stop("Could not archive the 256,529-row matrix folder.")
    }
    message("Archived the prior matrices under data/matrix_256529/.")
  } else if (length(first_ids) == 393628L) {
    message("The current matrices already have the manuscript's 393,628 rows.")
  } else {
    stop("Current matrices have an unexpected number of rows: ", length(first_ids))
  }
}
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

split_files <- list.files(
  project_file("data", "splits"),
  pattern = "^split_[0-9]+[.]csv$",
  full.names = TRUE
)
split_numbers <- as.integer(sub(
  "^split_([0-9]+)[.]csv$", "\\1", basename(split_files)
))
split_files <- split_files[order(split_numbers)]
if (length(split_files) != 10L || !identical(split_numbers[order(split_numbers)], 1:10)) {
  stop("Expected the validated split_1.csv through split_10.csv files.")
}

manuscript_sites <- readRDS(
  project_file("intermediates", "POST_sites_18of23_datasets.rds")
)
sample_coverage_sites <- readRDS(
  project_file("intermediates", "POST_sites_4176_samples.rds")
)
if (length(manuscript_sites) != 393628L ||
    length(sample_coverage_sites) != 256529L) {
  stop("The preprocessing CpG checkpoints do not have the validated counts.")
}

saveRDS(manuscript_sites, project_file("intermediates", "POST_sites.rds"))
full_metadata <- readRDS(project_file("intermediates", "full_metadata.rds"))

for (gender in c("all", "male", "female")) {
  for (average in c(TRUE, FALSE)) {
    generateMatrix(
      input_folder = project_file("data", "splits"),
      cpg_list = manuscript_sites,
      full_metadata = full_metadata,
      gender = gender,
      average = average,
      output_folder = matrix_dir,
      split_files = split_files
    )
  }
}

for (path in file.path(matrix_dir, matrix_names)) {
  ids <- data.table::fread(path, select = 1L, showProgress = FALSE)[[1L]]
  if (!identical(as.character(ids), as.character(manuscript_sites))) {
    stop("Matrix rows do not match the manuscript CpG background: ", path)
  }
}

cat("\nMANUSCRIPT MATRICES REBUILT AND VALIDATED\n")
cat("- Matrices: 6\n")
cat("- CpG rows per matrix: ", length(manuscript_sites), "\n", sep = "")
cat("- Prior 256,529-row matrices: data/matrix_256529/\n")
