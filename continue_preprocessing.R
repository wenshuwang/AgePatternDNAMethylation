# Continue after run_preprocessing.R has rebuilt and validated all 23 PP files.
# This step creates the large split files and the six age-summary matrices.

source("R/project_paths.R")
set_project_root()
source("R/preprocessing_core.R")

required <- c(
  "intermediates/POST_sites_18of23_datasets.rds",
  "intermediates/POST_sites_4176_samples.rds",
  "intermediates/full_metadata.rds",
  "intermediates/valid_samples.rds"
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Preprocessing checkpoints are missing. Run source('run_preprocessing.R') first:\n- ",
    paste(missing, collapse = "\n- ")
  )
}

pp_files <- list.files(
  project_file("data", "pp_datasets"),
  pattern = "^GSE[0-9]+_pp[.]csv$",
  full.names = TRUE
)
if (length(pp_files) != 23L) {
  stop("Expected 23 validated GSE*_pp.csv files; found ", length(pp_files), ".")
}

full_metadata <- readRDS(project_file("intermediates", "full_metadata.rds"))
valid_samples <- readRDS(project_file("intermediates", "valid_samples.rds"))
all_sites <- loadStandardCpGList()
manuscript_sites <- readRDS(
  project_file("intermediates", "POST_sites_18of23_datasets.rds")
)
if (length(manuscript_sites) != 393628L) {
  stop("The manuscript matrix background must contain exactly 393,628 CpGs.")
}
saveRDS(manuscript_sites, project_file("intermediates", "POST_sites.rds"))

split_files <- makeSplits(
  cpg_list = all_sites,
  full_metadata = full_metadata,
  valid_samples = valid_samples,
  sites_per_list = 50000L,
  output_folder = project_file("data", "splits")
)

for (gender in c("all", "male", "female")) {
  for (average in c(TRUE, FALSE)) {
    generateMatrix(
      input_folder = project_file("data", "splits"),
      cpg_list = manuscript_sites,
      full_metadata = full_metadata,
      gender = gender,
      average = average,
      output_folder = project_file("data", "matrix"),
      split_files = split_files
    )
  }
}

message("Split and matrix generation finished.")
