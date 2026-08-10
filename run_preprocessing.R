run_preprocessing <- function() {
  source("R/project_paths.R")
  set_project_root()

  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Install the knitr package before running preprocessing: install.packages('knitr')")
  }

  # A reproduction run should stop before expensive split/matrix generation if
  # the raw inputs do not reproduce the manuscript checkpoints.
  old_count_option <- getOption("agepattern.require_original_counts")
  options(agepattern.require_original_counts = TRUE)
  on.exit(options(agepattern.require_original_counts = old_count_option), add = TRUE)

  script_path <- tempfile("agepattern_preprocessing_", fileext = ".R")
  on.exit(unlink(script_path), add = TRUE)

  knitr::purl(
    input = "analysis/01_preprocessing.Rmd",
    output = script_path,
    documentation = 0,
    quiet = TRUE
  )
  source(script_path, echo = FALSE)
}

run_preprocessing()
