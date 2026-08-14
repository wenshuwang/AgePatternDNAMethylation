# Rebuild the manuscript's 393,628-CpG dataset-presence list from source
# methylation matrices. This step uses probe IDs only; it does not read any
# saved CpG lists or files under comparison/.

source("R/project_paths.R")
set_project_root()

local_input_dir <- Sys.getenv("AGEPATTERN_ORIGINAL_DATA_DIR", unset = "")
if (!nzchar(local_input_dir)) {
  local_input_dir <- project_file("data", "original_datasets")
}

checkpoint_dir <- project_file("intermediates", "presence_checkpoints")
geo_cache_dir <- file.path(checkpoint_dir, "geo_cache")
dir.create(geo_cache_dir, recursive = TRUE, showWarnings = FALSE)

local_gses <- c(
  "GSE32148", "GSE36054", "GSE40279", "GSE50660", "GSE50759",
  "GSE51057", "GSE53740", "GSE61256", "GSE67705", "GSE73103",
  "GSE80261", "GSE85568", "GSE89253", "GSE90124", "GSE94734",
  "GSE106648", "GSE114134", "GSE124366", "GSE137495", "GSE138279"
)
geo_gses <- c("GSE62924", "GSE51180", "GSE30870")
all_gses <- c(local_gses, geo_gses)

checkpoint_path <- function(gse) {
  file.path(checkpoint_dir, paste0(gse, "_original_probe_ids.rds"))
}

clean_probe_ids <- function(ids, source_name) {
  ids <- trimws(as.character(ids))
  keep <- grepl("^(cg[0-9]{8}|ch[.]|rs[0-9]+|ctl_)", ids, ignore.case = TRUE)
  ids <- unique(ids[!is.na(ids) & nzchar(ids) & keep])
  if (!length(ids)) stop("No probe IDs were read from ", source_name)
  ids
}

find_local_input <- function(gse) {
  candidates <- list.files(
    local_input_dir,
    pattern = paste0(
      "^o?", gse,
      "_(beta|methylation_data)[.](csv|xlsx|txt|txt[.]gz)$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Also accept the repository's normal Raw/GSE*/methylation_data.csv layout.
  raw_candidate <- project_file("data", "Raw", gse, "methylation_data.csv")
  candidates <- unique(c(candidates, raw_candidate[file.exists(raw_candidate)]))
  candidates <- candidates[file.exists(candidates)]

  if (!length(candidates)) {
    stop(
      "Missing source methylation matrix for ", gse, ".\n",
      "Place its o", gse, "_beta.csv/xlsx file under:\n  ",
      local_input_dir, "\n",
      "or use data/Raw/", gse, "/methylation_data.csv.\n",
      "Completed checkpoints are safe; rerun this script after adding the file."
    )
  }
  candidates[[1L]]
}

read_first_column <- function(path) {
  if (grepl("[.]xlsx$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Install readxl first: install.packages('readxl')")
    }
    sheets <- readxl::excel_sheets(path)
    return(unlist(lapply(sheets, function(sheet) {
      message("Reading first column of ", basename(path), ", sheet ", sheet)
      values <- readxl::read_excel(
        path,
        sheet = sheet,
        range = readxl::cell_cols(1L),
        col_names = TRUE,
        progress = TRUE
      )
      values[[1L]]
    }), use.names = FALSE))
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table first: install.packages('data.table')")
  }
  data.table::fread(path, select = 1L, showProgress = TRUE)[[1L]]
}

read_lifted_gse114134 <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table first: install.packages('data.table')")
  }
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Install sesame before processing GSE114134.")
  }
  if (grepl("[.]xlsx$", path, ignore.case = TRUE)) {
    stop("GSE114134 must be supplied as CSV or TXT so its full matrix can be lifted.")
  }

  message("Reading and lifting GSE114134 to the HM450 probe coordinates.")
  input <- data.table::fread(path, showProgress = TRUE)
  probe_ids <- as.character(input[[1L]])
  input[[1L]] <- NULL
  beta <- as.matrix(input)
  rownames(beta) <- probe_ids
  rm(input, probe_ids)
  gc()

  lifted <- sesame::mLiftOver(beta, "HM450", impute = FALSE)
  sites <- rownames(lifted)
  rm(beta, lifted)
  gc()
  sites
}

read_geo_sites <- function(gse) {
  if (!requireNamespace("GEOquery", quietly = TRUE) ||
      !requireNamespace("Biobase", quietly = TRUE)) {
    stop("Install GEOquery and Biobase before processing ", gse, ".")
  }

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
      return(rownames(Biobase::exprs(object[[1L]])))
    }
    if (attempt < 5) Sys.sleep(attempt * 5)
  }
  stop("GEO download failed for ", gse, ": ", conditionMessage(last_error))
}

for (gse in all_gses) {
  output <- checkpoint_path(gse)
  if (file.exists(output)) {
    message("Already complete: ", gse)
    next
  }

  message("Processing ", gse)
  sites <- if (gse %in% geo_gses) {
    read_geo_sites(gse)
  } else {
    input_path <- find_local_input(gse)
    if (gse == "GSE114134") {
      read_lifted_gse114134(input_path)
    } else {
      read_first_column(input_path)
    }
  }
  sites <- clean_probe_ids(sites, gse)
  saveRDS(sites, output)
  message("Saved ", length(sites), " probe IDs for ", gse)
  rm(sites)
  gc()
}

# Count how many of the 23 datasets contain each probe. The expected number is
# checked only after the list is calculated; it is never used to select probes.
probe_lists <- lapply(checkpoint_path(all_gses), readRDS)
probe_counts <- table(unlist(probe_lists, use.names = FALSE))
running_site_list <- data.frame(
  Count = as.integer(probe_counts),
  row.names = names(probe_counts)
)
sites_18_of_23 <- rownames(running_site_list)[running_site_list$Count >= 18L]
cpgs_18_of_23 <- sites_18_of_23[
  grepl("^cg[0-9]{8}$", sites_18_of_23, ignore.case = TRUE)
]

dir.create(project_file("intermediates"), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  running_site_list,
  project_file("intermediates", "running_site_list_from_original_inputs.rds")
)
saveRDS(
  cpgs_18_of_23,
  project_file("intermediates", "POST_sites_18of23_datasets.rds")
)

cat("\nDataset-presence rebuild complete\n")
cat("- Datasets:", length(probe_lists), "\n")
cat("- Unique input probe IDs:", nrow(running_site_list), "\n")
cat("- CpGs present in at least 18/23 datasets:", length(cpgs_18_of_23), "\n")

if (nrow(running_site_list) != 486427L || length(cpgs_18_of_23) != 393628L) {
  stop(
    "The source inputs did not reproduce the manuscript checkpoints. ",
    "Expected 486,427 unique probe IDs and 393,628 CpGs at >=18/23."
  )
}
message("Exact manuscript checkpoint reproduced: 393,628 CpGs.")
