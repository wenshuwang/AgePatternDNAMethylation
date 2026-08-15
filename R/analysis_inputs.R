normalize_site_ids <- function(x, label = "CpG list") {
  ids <- if (is.data.frame(x) || is.matrix(x)) rownames(x) else as.character(x)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) stop(label, " does not contain any usable probe IDs.")
  unique(ids)
}

count_values <- function(x, label) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop(label, " must be a data frame or matrix with a Count column.")
  }
  column <- if ("Count" %in% colnames(x)) "Count" else colnames(x)[1]
  as.numeric(x[, column])
}

analysis_input_mode <- function(mode = Sys.getenv("AGEPATTERN_INPUT_MODE", "auto")) {
  mode <- tolower(trimws(mode))
  if (mode %in% c("original", "og")) mode <- "archived"
  match.arg(mode, c("auto", "local", "archived"))
}

analysis_input_files <- function(archived_dir = project_file("comparison", "kevin_old")) {
  list(
    local = list(
      post_sites = project_file("intermediates", "POST_sites.rds"),
      dataset_presence_sites = project_file(
        "intermediates", "POST_sites_18of23_datasets.rds"
      ),
      sample_coverage_sites = project_file(
        "intermediates", "POST_sites_4176_samples.rds"
      ),
      matrices = c(
        all_mean = project_file("data", "matrix", "allTRUE.csv"),
        all_sd = project_file("data", "matrix", "allFALSE.csv"),
        male_mean = project_file("data", "matrix", "maleTRUE.csv"),
        male_sd = project_file("data", "matrix", "maleFALSE.csv"),
        female_mean = project_file("data", "matrix", "femaleTRUE.csv"),
        female_sd = project_file("data", "matrix", "femaleFALSE.csv")
      )
    ),
    archived = list(
      running_sites = file.path(archived_dir, "running_site_list.rds"),
      sample_sites = file.path(archived_dir, "site_list_by_samples.rds"),
      valid_samples = file.path(archived_dir, "valid_samples.rds"),
      post_sites = file.path(archived_dir, "POST_sites.rds"),
      matrices = c(
        all_mean = file.path(archived_dir, "allTRUE_OG.csv"),
        all_sd = file.path(archived_dir, "allFALSE_OG.csv"),
        male_mean = file.path(archived_dir, "maleTRUE_OG.csv"),
        male_sd = file.path(archived_dir, "maleFALSE_OG.csv"),
        female_mean = file.path(archived_dir, "femaleTRUE_OG.csv"),
        female_sd = file.path(archived_dir, "femaleFALSE_OG.csv")
      )
    )
  )
}

validate_local_inputs <- function(specification) {
  missing <- required_input_paths(specification)[
    !file.exists(required_input_paths(specification))
  ]
  if (length(missing)) {
    stop(
      "Local analysis mode is missing required files:\n- ",
      paste(missing, collapse = "\n- ")
    )
  }

  filtered_sites <- normalize_site_ids(
    readRDS(specification$post_sites),
    "Local POST_sites.rds"
  )
  dataset_presence_sites <- normalize_site_ids(
    readRDS(specification$dataset_presence_sites),
    "Local POST_sites_18of23_datasets.rds"
  )
  sample_coverage_sites <- normalize_site_ids(
    readRDS(specification$sample_coverage_sites),
    "Local POST_sites_4176_samples.rds"
  )

  if (length(dataset_presence_sites) != 393628L) {
    stop(
      "Local dataset-coverage list has ", length(dataset_presence_sites),
      " CpGs; the original raw run requires 393,628."
    )
  }
  if (length(sample_coverage_sites) != 256529L) {
    stop(
      "Local sample-coverage checkpoint has ", length(sample_coverage_sites),
      " CpGs; expected 256,529."
    )
  }
  if (length(filtered_sites) != 393628L ||
      !identical(filtered_sites, dataset_presence_sites)) {
    stop(
      "Local POST_sites.rds must be the 393,628-site manuscript matrix background."
    )
  }
  if (!all(sample_coverage_sites %in% dataset_presence_sites)) {
    stop("The 4,176-sample CpGs are not all present in the manuscript background.")
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required to validate local matrices.")
  }
  for (matrix_path in specification$matrices) {
    matrix_sites <- as.character(
      data.table::fread(matrix_path, select = 1L, showProgress = FALSE)[[1L]]
    )
    if (!identical(matrix_sites, filtered_sites)) {
      stop(
        "Matrix rows do not match local POST_sites.rds: ", matrix_path,
        ". Rebuild the matrices from the validated raw preprocessing outputs."
      )
    }
  }

  list(
    filtered_sites = filtered_sites,
    dataset_presence_sites = dataset_presence_sites,
    observed_counts = c(
      dataset_18_of_23 = length(dataset_presence_sites),
      matrix_sites = length(filtered_sites),
      sample_coverage_sites = length(sample_coverage_sites)
    )
  )
}

required_input_paths <- function(specification) {
  unname(unlist(specification, use.names = FALSE))
}

validate_archived_inputs <- function(specification) {
  missing <- required_input_paths(specification)[
    !file.exists(required_input_paths(specification))
  ]
  if (length(missing)) {
    stop(
      "Archived-original mode is missing required files:\n- ",
      paste(missing, collapse = "\n- ")
    )
  }

  running_sites <- readRDS(specification$running_sites)
  sample_sites <- readRDS(specification$sample_sites)
  valid_samples <- readRDS(specification$valid_samples)
  saved_post_sites <- normalize_site_ids(
    readRDS(specification$post_sites),
    "Archived POST_sites.rds"
  )

  running_counts <- count_values(running_sites, "Archived running_site_list.rds")
  sample_counts <- count_values(sample_sites, "Archived site_list_by_samples.rds")

  dataset_presence_sites <- rownames(running_sites)[
    which(!is.na(running_counts) & running_counts >= 18)
  ]
  analysis_sites <- rownames(sample_sites)[
    which(!is.na(sample_counts) & sample_counts >= 4176)
  ]

  expected <- c(
    running_rows = 486427L,
    split_universe = 485577L,
    dataset_18_of_23 = 393628L,
    analysis_sites = 256529L,
    valid_samples = 4641L
  )
  observed <- c(
    running_rows = nrow(running_sites),
    split_universe = nrow(sample_sites),
    dataset_18_of_23 = length(dataset_presence_sites),
    analysis_sites = length(analysis_sites),
    valid_samples = length(valid_samples)
  )

  if (!identical(as.integer(observed), as.integer(expected))) {
    details <- paste0(names(observed), ": expected ", expected, ", found ", observed)
    stop(
      "Archived checkpoints do not match the original run:\n- ",
      paste(details, collapse = "\n- ")
    )
  }
  if (!setequal(analysis_sites, saved_post_sites)) {
    stop("Archived POST_sites.rds does not match the 4,176-sample filter.")
  }

  list(
    filtered_sites = analysis_sites,
    dataset_presence_sites = dataset_presence_sites,
    observed_counts = observed
  )
}

resolve_analysis_inputs <- function(
    mode = Sys.getenv("AGEPATTERN_INPUT_MODE", "auto"),
    archived_dir = project_file("comparison", "kevin_old")
) {
  mode <- analysis_input_mode(mode)
  files <- analysis_input_files(archived_dir)
  local_ready <- all(file.exists(required_input_paths(files$local)))
  archived_ready <- all(file.exists(required_input_paths(files$archived)))

  if (mode == "auto") {
    # A complete fresh/local run is the default. Archived files are only a
    # fallback when the locally generated matrices do not exist.
    mode <- if (local_ready) "local" else if (archived_ready) "archived" else ""
    if (!nzchar(mode)) {
      stop(
        "No complete filtering input set was found. Provide either:\n",
        "- local matrices under data/matrix plus intermediates/POST_sites.rds, or\n",
        "- original archived files under comparison/kevin_old."
      )
    }
  }

  if (mode == "archived") {
    validation <- validate_archived_inputs(files$archived)
    message(
      "Analysis input mode: archived original run ",
      "(393,628 dataset-coverage CpGs; 256,529 downstream CpGs)."
    )
    return(list(
      mode = mode,
      matrices = files$archived$matrices,
      filtered_sites = validation$filtered_sites,
      dataset_presence_sites = validation$dataset_presence_sites,
      counts = validation$observed_counts
    ))
  }

  validation <- validate_local_inputs(files$local)
  message(
    "Analysis input mode: validated local raw-data run ",
    "(393,628-CpG manuscript matrices; 256,529-site diagnostic checkpoint)."
  )
  list(
    mode = mode,
    matrices = files$local$matrices,
    filtered_sites = validation$filtered_sites,
    dataset_presence_sites = validation$dataset_presence_sites,
    counts = validation$observed_counts
  )
}
