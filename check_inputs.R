source("project_paths.R")
set_project_root()

gse_ids <- c("32148", "36054", "40279", "50660", "50759", "51057", "53740", "61256", "67705", "73103", "80261", "85568", "89253", "90124", "94734", "106648", "114134", "124366", "137495", "138279", "62924", "51180", "30870")
required <- c(
  unlist(lapply(gse_ids, function(id) file.path("data", "Raw", paste0("GSE", id), c("methylation_data.csv", "metadata.csv")))),
  file.path("data", "matrix", c("allTRUE.csv", "allFALSE.csv", "maleTRUE.csv", "maleFALSE.csv", "femaleTRUE.csv", "femaleFALSE.csv")),
  file.path("relevant_rds", c("all_results.rds", "male_results.rds", "female_results.rds"))
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  cat("Missing inputs:\n", paste0("- ", missing, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("All raw datasets, matrix files, and supplied filtering results are available.\n")
