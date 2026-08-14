run_preprocessing <- function() {
  source("R/project_paths.R")
  set_project_root()

  # First rebuild dataset-level probe presence directly from the 23 source
  # matrices. This must reproduce the independent 393,628-CpG checkpoint.
  source("rebuild_18of23_sites.R", local = environment())

  # Then rebuild the 23 PP matrices one at a time from those same original
  # sources and verify the 4,176-sample checkpoint. This intentionally stops
  # before split/matrix generation because those files require substantial
  # additional disk space.
  source("rebuild_preprocessed_datasets.R", local = environment())
}

run_preprocessing()
