# Rebuild the 23 canonical preprocessed (PP) methylation matrices from the
# original source matrices. This script is resumable: a completed, validated
# GSE*_pp.csv file is skipped when the script is run again.

source("R/project_paths.R")
set_project_root()
source("R/preprocessing_core.R")

required_packages <- c(
  "data.table", "readxl", "GEOquery", "Biobase",
  "sesame", "sesameData", "BiocParallel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Install the missing packages before rebuilding PP files: ",
    paste(missing_packages, collapse = ", ")
  )
}

local_input_dir <- Sys.getenv("AGEPATTERN_ORIGINAL_DATA_DIR", unset = "")
if (!nzchar(local_input_dir)) {
  local_input_dir <- project_file("data", "original_datasets")
}
metadata_dir <- Sys.getenv("AGEPATTERN_METADATA_DIR", unset = "")
if (!nzchar(metadata_dir)) metadata_dir <- project_file("data", "metadata")

pp_dir <- project_file("data", "pp_datasets")
geo_cache_dir <- project_file("intermediates", "geo_cache")
dir.create(pp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(project_file("intermediates"), recursive = TRUE, showWarnings = FALSE)

local_gses <- c(
  "GSE32148", "GSE36054", "GSE40279", "GSE50660", "GSE50759",
  "GSE51057", "GSE53740", "GSE61256", "GSE67705", "GSE73103",
  "GSE80261", "GSE85568", "GSE89253", "GSE90124", "GSE94734",
  "GSE106648", "GSE114134", "GSE124366", "GSE137495", "GSE138279"
)
geo_gses <- c("GSE62924", "GSE51180", "GSE30870")
all_gses <- c(local_gses, geo_gses)

ignored_rebuilt_files <- list.files(
  pp_dir,
  pattern = "_pp_REBUILT[.]csv$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(ignored_rebuilt_files)) {
  warning(
    "Found ", length(ignored_rebuilt_files),
    " old *_pp_REBUILT.csv files. They will not be used, but they occupy disk space."
  )
}

standard_cpg_list <- loadStandardCpGList()
if (length(standard_cpg_list) != 485577L) {
  stop(
    "The HM450 manifest contains ", length(standard_cpg_list),
    " unique probes; expected 485,577."
  )
}

read_cleaned_metadata <- function(gse) {
  candidates <- c(
    file.path(metadata_dir, paste0(gse, "_cleaned_metadata.csv")),
    file.path(metadata_dir, paste0(gse, "_cleaned_metadata.xlsx"))
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) != 1L) {
    stop(
      "Expected one cleaned metadata file for ", gse, " under:\n  ",
      metadata_dir
    )
  }

  path <- candidates[[1L]]
  metadata <- if (grepl("[.]xlsx$", path, ignore.case = TRUE)) {
    as.data.frame(readxl::read_excel(path))
  } else {
    as.data.frame(data.table::fread(path, showProgress = FALSE))
  }
  normalized_names <- tolower(gsub("[^a-z0-9]+", "_", names(metadata)))
  pick_column <- function(candidates, required = FALSE) {
    hit <- match(candidates, normalized_names, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) return(names(metadata)[hit[[1L]]])
    if (required) {
      stop(
        basename(path), " is missing a required column matching: ",
        paste(candidates, collapse = ", ")
      )
    }
    NA_character_
  }
  values_or_na <- function(column) {
    if (is.na(column)) rep(NA_character_, nrow(metadata)) else metadata[[column]]
  }

  sample_column <- pick_column(
    c("sample_id", "status_id", "geo_accession", "gsm", "id"),
    required = TRUE
  )
  age_column <- pick_column("age", required = TRUE)
  gender_column <- pick_column(c("gender", "sex"))
  type_column <- pick_column(c("type", "sample_type"))
  tissue_column <- pick_column(c("tissue", "cell_type", "source_name_ch1"))
  disease_column <- pick_column(c(
    "disease_status", "disease_state", "disease", "diagnosis", "status"
  ))

  metadata <- data.frame(
    sample_id = as.character(metadata[[sample_column]]),
    age = round(suppressWarnings(as.numeric(as.character(metadata[[age_column]])))),
    gender = normalizeGenderLabels(values_or_na(gender_column)),
    type = as.character(values_or_na(type_column)),
    tissue = as.character(values_or_na(tissue_column)),
    disease_status = as.character(values_or_na(disease_column)),
    stringsAsFactors = FALSE
  )
  metadata$dataset <- gse
  metadata <- metadata[
    !is.na(metadata$sample_id) & nzchar(metadata$sample_id), , drop = FALSE
  ]
  metadata
}

full_metadata <- do.call(rbind, lapply(all_gses, read_cleaned_metadata))
rownames(full_metadata) <- NULL
if (anyDuplicated(full_metadata$sample_id)) {
  duplicated_ids <- unique(full_metadata$sample_id[duplicated(full_metadata$sample_id)])
  stop(
    "Cleaned metadata contains duplicated sample IDs: ",
    paste(head(duplicated_ids, 10L), collapse = ", ")
  )
}
valid_samples <- validSampleIds(full_metadata, 80)

if (nrow(full_metadata) != 4946L || length(valid_samples) != 4641L) {
  stop(
    "Cleaned metadata does not reproduce the original sample checkpoints.\n",
    "Expected 4,946 metadata rows and 4,641 valid age-0-to-80 samples; found ",
    nrow(full_metadata), " and ", length(valid_samples), "."
  )
}

saveRDS(full_metadata, project_file("intermediates", "full_metadata.rds"))
saveRDS(valid_samples, project_file("intermediates", "valid_samples.rds"))

find_local_beta <- function(gse) {
  candidates <- list.files(
    local_input_dir,
    pattern = paste0(
      "^o?", gse,
      "_(beta|methylation_data)[.](csv|xlsx|txt|txt[.]gz)$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )
  raw_candidate <- project_file("data", "Raw", gse, "methylation_data.csv")
  candidates <- unique(c(candidates, raw_candidate[file.exists(raw_candidate)]))
  candidates <- candidates[
    file.exists(candidates) & !is.na(file.info(candidates)$size) &
      file.info(candidates)$size > 0
  ]
  if (length(candidates) != 1L) {
    stop(
      "Expected one original beta matrix for ", gse, " under:\n  ",
      local_input_dir, "\nFound ", length(candidates), "."
    )
  }
  candidates[[1L]]
}

read_xlsx_matrix <- function(path) {
  sheets <- readxl::excel_sheets(path)
  data.table::rbindlist(
    lapply(sheets, function(sheet) {
      message("Reading ", basename(path), ", sheet ", sheet)
      data.table::as.data.table(
        readxl::read_excel(path, sheet = sheet, progress = TRUE)
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
}

read_local_beta <- function(path, gse, gse_metadata) {
  extension <- tolower(tools::file_ext(path))
  dataset <- if (extension == "xlsx") {
    read_xlsx_matrix(path)
  } else if (extension %in% c("csv", "txt", "gz")) {
    data.table::fread(path, showProgress = TRUE)
  } else {
    stop("Unsupported source-matrix type: ", path)
  }

  # Two official source tables use generic Sample 1 ... Sample N column names.
  # Restore GSM names only when no GSM names are already present and the
  # cleaned metadata has exactly the same number of samples.
  matrix_samples <- names(dataset)[-1L]
  metadata_samples <- as.character(gse_metadata$sample_id)
  if (gse %in% c("GSE40279", "GSE50660") &&
      !length(intersect(matrix_samples, metadata_samples)) &&
      length(matrix_samples) == length(metadata_samples)) {
    data.table::setnames(dataset, matrix_samples, metadata_samples)
  }
  dataset
}

read_geo_beta <- function(gse) {
  last_error <- NULL
  for (attempt in 1:5) {
    message("GEO download attempt ", attempt, "/5 for ", gse)
    object <- tryCatch(
      GEOquery::getGEO(gse, GSEMatrix = TRUE, destdir = geo_cache_dir),
      error = function(error) {
        last_error <<- error
        NULL
      }
    )
    if (!is.null(object)) {
      values <- Biobase::exprs(object[[1L]])
      dataset <- data.table::as.data.table(values, keep.rownames = "probe_id")
      rm(values, object)
      gc()
      return(dataset)
    }
    if (attempt < 5L) Sys.sleep(attempt * 5L)
  }
  stop("GEO download failed for ", gse, ": ", conditionMessage(last_error))
}

validate_existing_pp <- function(path, expected_samples) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(FALSE)
  }
  header <- names(data.table::fread(path, nrows = 0L, check.names = FALSE))
  if (!setequal(header[-1L], expected_samples)) return(FALSE)
  probe_ids <- data.table::fread(path, select = 1L, showProgress = FALSE)[[1L]]
  identical(as.character(probe_ids), standard_cpg_list)
}

prepare_beta <- function(dataset, gse, expected_samples) {
  probe_ids <- trimws(as.character(dataset[[1L]]))
  dataset[[1L]] <- NULL
  if (anyNA(probe_ids) || any(!nzchar(probe_ids)) || anyDuplicated(probe_ids)) {
    stop(gse, " has missing, blank, or duplicated probe IDs.")
  }

  missing_samples <- setdiff(expected_samples, names(dataset))
  if (length(missing_samples)) {
    stop(
      gse, " is missing ", length(missing_samples), " valid metadata samples. First: ",
      paste(head(missing_samples, 10L), collapse = ", ")
    )
  }
  dataset <- dataset[, expected_samples, with = FALSE]
  beta <- as.matrix(dataset)
  storage.mode(beta) <- "double"
  rownames(beta) <- probe_ids
  rm(dataset, probe_ids)
  gc()

  if (gse %in% c("GSE50759", "GSE89253")) beta <- m_to_b(beta)
  if (gse == "GSE85568") beta <- limit_bounds(beta)
  if (gse == "GSE114134") beta <- liftover(beta)

  beta <- beta[rownames(beta) %in% standard_cpg_list, , drop = FALSE]
  if (anyDuplicated(rownames(beta))) {
    stop(gse, " has duplicated HM450 probe IDs after preprocessing.")
  }
  beta
}

for (index in seq_along(all_gses)) {
  gse <- all_gses[[index]]
  gse_metadata <- full_metadata[full_metadata$dataset == gse, , drop = FALSE]
  expected_samples <- validSampleIds(gse_metadata, 80)
  output_path <- file.path(pp_dir, paste0(gse, "_pp.csv"))

  if (validate_existing_pp(output_path, expected_samples)) {
    message("Already complete (", index, "/23): ", gse)
    next
  }
  if (file.exists(output_path)) {
    stop(
      "An incomplete or incompatible PP file already exists:\n  ", output_path,
      "\nMove or delete only that file, then rerun. Completed datasets are safe."
    )
  }

  message("Rebuilding ", index, "/23: ", gse)
  dataset <- if (gse %in% geo_gses) {
    read_geo_beta(gse)
  } else {
    input_path <- find_local_beta(gse)
    message("Source: ", input_path)
    read_local_beta(input_path, gse, gse_metadata)
  }
  beta <- prepare_beta(dataset, gse, expected_samples)
  rm(dataset)
  gc()

  output <- matrix(
    NA_real_,
    nrow = length(standard_cpg_list),
    ncol = ncol(beta),
    dimnames = list(standard_cpg_list, colnames(beta))
  )
  destination <- match(rownames(beta), standard_cpg_list)
  output[destination, ] <- beta
  rm(beta, destination)
  gc()

  partial_path <- paste0(output_path, ".part")
  if (file.exists(partial_path)) unlink(partial_path)
  data.table::fwrite(
    output,
    partial_path,
    row.names = TRUE,
    showProgress = TRUE
  )
  rm(output)
  gc()
  if (!file.rename(partial_path, output_path)) {
    stop("The PP file was written but could not be finalized: ", partial_path)
  }
  if (!validate_existing_pp(output_path, expected_samples)) {
    stop("Validation failed after writing ", output_path)
  }
  message("Finished and validated ", gse)
}

site_list_by_samples <- rebuildSampleCoverageFromPreprocessed(
  pp_dir = pp_dir,
  valid_samples = valid_samples,
  expected_sites = standard_cpg_list
)

independent_running_path <- project_file(
  "intermediates", "running_site_list_from_original_inputs.rds"
)
if (!file.exists(independent_running_path)) {
  stop(
    "Missing the independent source-presence checkpoint. Run ",
    "rebuild_18of23_sites.R first."
  )
}
running_site_list <- readRDS(independent_running_path)

metadata_tracker <- data.frame(
  matrix(0L, nrow = 101L, ncol = length(all_gses)),
  row.names = paste0("age", 0:100)
)
names(metadata_tracker) <- all_gses
for (gse in all_gses) {
  ages <- full_metadata$age[full_metadata$dataset == gse]
  counts <- table(ages[!is.na(ages) & ages >= 0 & ages <= 100])
  metadata_tracker[paste0("age", names(counts)), gse] <- as.integer(counts)
}

saveRDS(running_site_list, project_file("intermediates", "running_site_list.rds"))
saveRDS(site_list_by_samples, project_file("intermediates", "site_list_by_samples.rds"))
saveRDS(metadata_tracker, project_file("intermediates", "metadata_tracker.rds"))

preprocessing_counts <- reportOriginalPreprocessingCounts(
  running_site_list,
  site_list_by_samples,
  valid_samples,
  require_original_counts = TRUE
)
dataset_presence_sites <- preprocessing_counts$lists$dataset_presence_sites
analysis_sites <- preprocessing_counts$lists$analysis_sites

saveRDS(
  dataset_presence_sites,
  project_file("intermediates", "POST_sites_18of23_datasets.rds")
)
saveRDS(
  analysis_sites,
  project_file("intermediates", "POST_sites_4176_samples.rds")
)
# The manuscript's age-summary matrices and aaCpG filtering use the 393,628
# CpGs present in at least 18 of 23 datasets. Keep the stricter 256,529-site
# sample-coverage checkpoint under its explicit filename for diagnostics.
saveRDS(dataset_presence_sites, project_file("intermediates", "POST_sites.rds"))

male_samples <- intersect(
  full_metadata$sample_id[full_metadata$gender == "male"], valid_samples
)
female_samples <- intersect(
  full_metadata$sample_id[full_metadata$gender == "female"], valid_samples
)
saveRDS(male_samples, project_file("intermediates", "male_samples.rds"))
saveRDS(female_samples, project_file("intermediates", "female_samples.rds"))

cat("\nPP REBUILD COMPLETE AND VALIDATED\n")
cat("- PP files: 23\n")
cat("- Valid samples: ", length(valid_samples), "\n", sep = "")
cat("- CpGs at >=18/23 datasets: ", length(dataset_presence_sites), "\n", sep = "")
cat("- CpGs at >=4,176 samples: ", length(analysis_sites), "\n", sep = "")
cat("Split and matrix files were not generated yet.\n")
