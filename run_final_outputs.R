# Generate the manuscript tables and figures from completed filtering outputs.
#
# Default behavior creates Table 1, Figure 3, Supplementary Figure 4, and the
# clustering-based Figures 2 and 5. To generate only the quick outputs, run:
# options(agepattern.run_clustering = FALSE)
# source("run_final_outputs.R")

if (!file.exists("R/project_paths.R")) {
  stop(
    "Run this file from the AgePatternDNAMethylation project folder.\n",
    "Example:\n",
    "setwd('C:/Users/wensh/Downloads/AgePatternDNAMethylation-main')\n",
    "source('run_final_outputs.R')",
    call. = FALSE
  )
}

source("R/project_paths.R")
set_project_root()
source("R/final_outputs.R")

output_dir <- "final_outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Creating Table 1, Figure 3, and Supplementary Figure 4...")
generate_quick_final_outputs(".", output_dir)

run_clustering <- isTRUE(getOption("agepattern.run_clustering", TRUE))

if (run_clustering) {
  message("Creating clustering Figures 2 and 5. This is the slower step...")

  required_clustering_inputs <- c(
    "analysis/03_clustering_helpers.Rmd",
    "analysis/03_clustering.Rmd",
    "data/matrix/allTRUE.csv",
    "data/matrix/maleTRUE.csv",
    "data/matrix/femaleTRUE.csv",
    "relevant_rds/all_select_cpg.rds",
    "relevant_rds/male_select_cpg.rds",
    "relevant_rds/female_select_cpg.rds",
    "relevant_rds/clock_cpgs.rds"
  )
  missing <- required_clustering_inputs[!file.exists(required_clustering_inputs)]
  if (length(missing)) {
    stop(
      "The quick outputs were saved, but clustering inputs are missing:\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }

  options(agepattern.output_dir = output_dir)

  helper_r <- tempfile(fileext = ".R")
  clustering_r <- tempfile(fileext = ".R")
  knitr::purl(
    "analysis/03_clustering_helpers.Rmd",
    output = helper_r,
    quiet = TRUE
  )
  sys.source(helper_r, envir = .GlobalEnv)

  knitr::purl(
    "analysis/03_clustering.Rmd",
    output = clustering_r,
    quiet = TRUE
  )
  sys.source(clustering_r, envir = .GlobalEnv)

  message("Clustering figures complete.")
} else {
  message(
    "Clustering was skipped. Table 1 and Figure 3 are complete.\n",
    "To create Figures 2 and 5 later, restart R and run:\n",
    "options(agepattern.run_clustering = TRUE)\n",
    "source('run_final_outputs.R')"
  )
}

message("FINAL OUTPUTS FINISHED: ", normalizePath(output_dir))
