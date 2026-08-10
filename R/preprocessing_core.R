suppressPackageStartupMessages(library(data.table))

readCSV <- function(file_path) {
  data.table::fread(file_path)
}

check_bounds <- function(dataset) {
  values <- as.matrix(dataset)
  values <- values[!is.na(values)]
  in_bounds <- !length(values) || all(values >= 0 & values <= 1)
  print(in_bounds)
  invisible(in_bounds)
}

limit_bounds <- function(dataset) {
  pmin(pmax(dataset, 0), 1)
}

m_to_b <- function(dataset) {
  2^as.matrix(dataset) / (1 + 2^as.matrix(dataset))
}

liftover <- function(dataset) {
  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("The sesame package is required for the GSE114134 liftover step.")
  }
  sesame::mLiftOver(dataset, "HM450", impute = FALSE)
}

validSampleIds <- function(md, age_threshold = 80) {
  if (!all(c("sample_id", "age") %in% colnames(md))) {
    stop("Metadata must contain sample_id and age columns before sample filtering.")
  }
  ages <- suppressWarnings(as.numeric(as.character(md$age)))
  ids <- as.character(md$sample_id)
  keep <- !is.na(ids) & nzchar(ids) & !is.na(ages) & ages >= 0 & ages <= age_threshold
  unique(ids[keep])
}

filterBetaToValidSamples <- function(dataset, md, age_threshold = 80) {
  valid_ids <- validSampleIds(md, age_threshold)
  dataset[, colnames(dataset) %in% valid_ids, drop = FALSE]
}

writePreprocessedDataset <- function(dataset, gse_id) {
  output_path <- project_file("data", "pp_datasets", paste0(gse_id, "_pp.csv"))
  data.table::fwrite(dataset, output_path, row.names = TRUE)
  invisible(output_path)
}

updateValidSamples <- function(valid_samples, md, age_threshold) {
  unique(c(as.character(valid_samples), validSampleIds(md, age_threshold)))
}

normalizeGenderLabels <- function(gender) {
  normalized <- trimws(tolower(as.character(gender)))
  normalized[normalized %in% c("m", "male")] <- "male"
  normalized[normalized %in% c("f", "female")] <- "female"
  normalized[is.na(gender) | !nzchar(normalized)] <- NA_character_
  normalized
}

updateMetadata <- function(updated, md, name) {
  age_range <- 0:100
  rounded_ages <- suppressWarnings(round(as.numeric(as.character(md$age))))
  rounded_ages <- rounded_ages[!is.na(rounded_ages) & rounded_ages >= 0 & rounded_ages <= 100]
  age_counts <- table(rounded_ages)
  dataset_column <- rep(0, length(age_range))
  names(dataset_column) <- paste0("age", age_range)
  for (age in names(age_counts)) {
    dataset_column[paste0("age", as.integer(age))] <- age_counts[[age]]
  }
  updated <- cbind(updated, dataset_column)
  colnames(updated)[ncol(updated)] <- name
  updated
}

loadStandardCpGList <- function(
    manifest_path = project_file("data", "ref", "HM450.hg38.manifest.tsv")
) {
  if (!file.exists(manifest_path)) stop("Missing HM450 manifest: ", manifest_path)
  manifest <- read.delim(
    manifest_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fileEncoding = "UTF-8"
  )
  if (!"probeID" %in% colnames(manifest)) {
    stop("HM450 manifest must contain a probeID column.")
  }
  ids <- as.character(manifest$probeID)
  unique(ids[!is.na(ids) & nzchar(ids)])
}

checkRawPreprocessingInputs <- function() {
  gse_ids <- c(
    "32148", "36054", "40279", "50660", "50759", "51057", "53740",
    "61256", "67705", "73103", "80261", "85568", "89253", "90124",
    "94734", "106648", "114134", "124366", "137495", "138279",
    "62924", "51180", "30870"
  )
  required <- unlist(lapply(gse_ids, function(id) {
    project_file(
      "data", "Raw", paste0("GSE", id),
      c("methylation_data.csv", "metadata.csv")
    )
  }))
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Fresh raw preprocessing is missing files:\n- ", paste(missing, collapse = "\n- "))
  }
  invisible(TRUE)
}

derivePreprocessingSiteLists <- function(running_site_list, site_list_by_samples) {
  running_counts <- as.numeric(running_site_list$Count)
  sample_counts <- as.numeric(site_list_by_samples$Count)
  list(
    dataset_presence_sites = rownames(running_site_list)[
      which(!is.na(running_counts) & running_counts >= 18)
    ],
    analysis_sites = rownames(site_list_by_samples)[
      which(!is.na(sample_counts) & sample_counts >= 4176)
    ]
  )
}

reportOriginalPreprocessingCounts <- function(
    running_site_list,
    site_list_by_samples,
    valid_samples,
    require_original_counts = getOption("agepattern.require_original_counts", FALSE)
) {
  lists <- derivePreprocessingSiteLists(running_site_list, site_list_by_samples)
  observed <- c(
    running_rows = nrow(running_site_list),
    split_universe = nrow(site_list_by_samples),
    dataset_18_of_23 = length(lists$dataset_presence_sites),
    analysis_sites = length(lists$analysis_sites),
    valid_samples = length(unique(valid_samples))
  )
  expected <- c(
    running_rows = 486427L,
    split_universe = 485577L,
    dataset_18_of_23 = 393628L,
    analysis_sites = 256529L,
    valid_samples = 4641L
  )
  print(observed)
  if (!identical(as.integer(observed), as.integer(expected))) {
    details <- paste0(names(observed), ": expected ", expected, ", found ", observed)
    message <- paste(
      "Fresh preprocessing does not match the original checkpoint counts:",
      paste0("- ", details, collapse = "\n"),
      sep = "\n"
    )
    if (isTRUE(require_original_counts)) stop(message) else warning(message)
  } else {
    message("Fresh preprocessing matches all original checkpoint counts.")
  }
  invisible(list(lists = lists, observed = observed, expected = expected))
}

makeSplits <- function(cpg_list, full_metadata, valid_samples, sites_per_list, output_folder) {
  file_list <- list.files(
    project_file("data", "pp_datasets"),
    pattern = "_pp\\.csv$",
    full.names = TRUE
  )
  if (length(file_list) != 23L) stop("Expected 23 preprocessed dataset files.")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  split_lists <- split(cpg_list, ceiling(seq_along(cpg_list) / sites_per_list))
  output_paths <- file.path(output_folder, paste0("split_", seq_along(split_lists), ".csv"))

  for (i in seq_along(split_lists)) {
    template_sites <- c("age", split_lists[[i]])
    combined <- data.frame(
      matrix(NA, nrow = length(template_sites), ncol = 0),
      row.names = template_sites
    )
    for (dataset_file in file_list) {
      message("Split ", i, "/", length(split_lists), ": ", basename(dataset_file))
      dataset <- as.data.frame(data.table::fread(dataset_file))
      rownames(dataset) <- dataset[, 1]
      dataset <- dataset[, -1, drop = FALSE]
      dataset <- dataset[rownames(dataset) %in% split_lists[[i]], , drop = FALSE]
      dataset <- dataset[, colnames(dataset) %in% valid_samples, drop = FALSE]
      metadata <- full_metadata[
        match(colnames(dataset), full_metadata$sample_id),
        c("sample_id", "age"),
        drop = FALSE
      ]
      if (anyNA(metadata$sample_id)) stop("A split sample is missing from full_metadata.")
      age_row <- matrix(
        as.numeric(metadata$age),
        nrow = 1,
        dimnames = list("age", colnames(dataset))
      )
      dataset <- rbind(age_row, dataset)
      dataset <- dataset[template_sites, , drop = FALSE]
      combined <- cbind(combined, dataset)
      rm(dataset, metadata, age_row)
      gc()
    }
    data.table::fwrite(
      combined,
      output_paths[[i]],
      row.names = TRUE
    )
    rm(combined)
    gc()
  }
  saveRDS(output_paths, project_file("intermediates", "preprocessing_split_files.rds"))
  invisible(output_paths)
}

createEmptyMatrix <- function(row_number, cpg_list, start_age, end_age) {
  output <- matrix(NA_real_, nrow = row_number, ncol = end_age - start_age + 1)
  rownames(output) <- cpg_list
  colnames(output) <- paste0("age_", start_age:end_age)
  output
}

generateMatrix <- function(
    input_folder,
    cpg_list,
    full_metadata,
    gender,
    average,
    output_folder,
    split_files = NULL
) {
  if (is.null(split_files)) {
    split_files <- list.files(
      input_folder,
      pattern = "^split_[0-9]+\\.csv$",
      full.names = TRUE
    )
    split_numbers <- as.integer(sub("^split_([0-9]+)\\.csv$", "\\1", basename(split_files)))
    split_files <- split_files[order(split_numbers)]
  }
  if (!length(split_files)) stop("No split CSV files found under ", input_folder)
  missing_split_files <- split_files[!file.exists(split_files)]
  if (length(missing_split_files)) {
    stop("Missing preprocessing split files:\n- ", paste(missing_split_files, collapse = "\n- "))
  }
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  valid_samples <- readRDS(project_file("intermediates", "valid_samples.rds"))
  output <- createEmptyMatrix(length(cpg_list), cpg_list, 0, 80)
  normalized_gender <- normalizeGenderLabels(full_metadata$gender)

  for (i in seq_along(split_files)) {
    dataset <- as.data.frame(data.table::fread(split_files[[i]]))
    rownames(dataset) <- dataset[, 1]
    dataset <- dataset[, -1, drop = FALSE]
    if (gender %in% c("male", "female")) {
      gender_samples <- full_metadata$sample_id[normalized_gender == gender]
      dataset <- dataset[, colnames(dataset) %in% intersect(gender_samples, valid_samples), drop = FALSE]
    }
    dataset["age", ] <- round(as.numeric(dataset["age", ]))
    dataset <- dataset[rownames(dataset) %in% c("age", cpg_list), , drop = FALSE]
    identifiers <- setdiff(rownames(dataset), "age")
    for (age in 0:80) {
      desired_columns <- which(dataset["age", ] == age)
      if (!length(desired_columns)) next
      values <- dataset[identifiers, desired_columns, drop = FALSE]
      numeric_values <- as.matrix(values)
      storage.mode(numeric_values) <- "numeric"
      statistic <- if (isTRUE(average)) {
        rowMeans(numeric_values, na.rm = TRUE)
      } else {
        apply(numeric_values, 1, stats::sd, na.rm = TRUE)
      }
      output[identifiers, paste0("age_", age)] <- statistic
    }
    message("Finished matrix split ", i, " of ", length(split_files),
            " for ", gender, ", average=", average)
  }
  data.table::fwrite(
    output,
    file.path(output_folder, paste0(gender, average, ".csv")),
    row.names = TRUE
  )
}
