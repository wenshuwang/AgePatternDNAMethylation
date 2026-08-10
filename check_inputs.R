source("R/project_paths.R")
set_project_root()
source("R/analysis_inputs.R")

check_project_status <- function() {
  gse_ids <- c(
    "32148", "36054", "40279", "50660", "50759", "51057", "53740",
    "61256", "67705", "73103", "80261", "85568", "89253", "90124",
    "94734", "106648", "114134", "124366", "137495", "138279",
    "62924", "51180", "30870"
  )

  raw_files <- unlist(lapply(gse_ids, function(id) {
    file.path(
      "data", "Raw", paste0("GSE", id),
      c("methylation_data.csv", "metadata.csv")
    )
  }))
  missing_raw <- raw_files[!file.exists(raw_files)]

  pp_files <- list.files(
    project_file("data", "pp_datasets"),
    pattern = "_pp\\.csv$",
    full.names = TRUE
  )

  cat("Preprocessing inputs\n")
  cat("- Complete raw dataset pairs:", length(raw_files) - length(missing_raw),
      "of", length(raw_files), "files\n")
  cat("- Preprocessed datasets:", length(pp_files), "of 23 expected files\n\n")

  analysis_inputs <- tryCatch(
    resolve_analysis_inputs(),
    error = function(error) error
  )
  analysis_ready <- !inherits(analysis_inputs, "error")

  cat("Filtering inputs\n")
  if (analysis_ready) {
    cat("- Ready\n")
    cat("- Selected mode:", analysis_inputs$mode, "\n")
    cat("- Downstream CpGs:", length(analysis_inputs$filtered_sites), "\n")
    if (analysis_inputs$mode == "archived") {
      cat("- Dataset-coverage CpGs:",
          length(analysis_inputs$dataset_presence_sites), "\n")
    }
  } else {
    cat("- Not ready:", conditionMessage(analysis_inputs), "\n")
  }

  clustering_files <- project_file(
    "relevant_rds",
    c("all_results.rds", "male_results.rds", "female_results.rds")
  )
  clustering_ready <- all(file.exists(clustering_files))
  cat("\nClustering prerequisites\n")
  cat("- Filtering result checkpoints:",
      sum(file.exists(clustering_files)), "of", length(clustering_files), "\n")

  if (length(missing_raw)) {
    cat(
      "\nRaw preprocessing is incomplete, but this does not block filtering ",
      "when a complete matrix/checkpoint set is available.\n",
      sep = ""
    )
  }

  invisible(list(
    raw_ready = !length(missing_raw),
    preprocessed_ready = length(pp_files) == 23L,
    analysis_ready = analysis_ready,
    analysis_mode = if (analysis_ready) analysis_inputs$mode else NA_character_,
    clustering_ready = clustering_ready,
    missing_raw = missing_raw
  ))
}

project_status <- check_project_status()
