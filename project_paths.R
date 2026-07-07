find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "README.md")) && dir.exists(file.path(current, "scripts"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not find the AgePatternDNAMethylation project root.")
    current <- parent
  }
}

set_project_root <- function() {
  root <- find_project_root()
  options(agepattern.project_root = root)
  setwd(root)
  invisible(root)
}

project_file <- function(...) {
  root <- getOption("agepattern.project_root")
  if (is.null(root)) root <- set_project_root()
  file.path(root, ...)
}

check_project_inputs <- function() {
  expected <- c("data/Raw", "data/pp_datasets", "data/matrix", "relevant_rds")
  missing <- expected[!dir.exists(project_file(expected))]
  if (length(missing)) stop("Missing project folders:\n- ", paste(missing, collapse = "\n- "))
  invisible(TRUE)
}
