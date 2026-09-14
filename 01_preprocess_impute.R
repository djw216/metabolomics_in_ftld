#!/usr/bin/env Rscript

# Metabolomics study:  imputation, scaling and RUV-III correction
#
# 
#
# Expected input files (place in data/raw, or set METABOLOMICS_INPUT_DIR):
#   non_imputed_2020.csv
#   non_imputed_2017.csv
#   non_imputed_2015.csv
#   subjects.csv
#   BNI_PEAK_CONC_METAB_LOOKUP.txt
#
# Output:
#   data/processed/metabolomics_preprocessed.rds
#

required_packages <- c("dplyr", "e1071", "missForest", "ruv")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the following packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(e1071)
  library(missForest)
  library(ruv)
})

set.seed(1234)

# ---- User-configurable paths and parameters ---------------------------------
project_dir <- normalizePath(
  Sys.getenv("METABOLOMICS_PROJECT_DIR", unset = "."),
  mustWork = FALSE
)
input_dir <- Sys.getenv(
  "METABOLOMICS_INPUT_DIR",
  unset = file.path(project_dir)
)
output_dir <- Sys.getenv(
  "METABOLOMICS_OUTPUT_DIR",
  unset = file.path(project_dir, "data", "processed")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

missingness_threshold <- 0.30
skewness_threshold <- 2
ruv_k <- 5L

# These compounds were excluded a priori in the original analysis. Document
# any changes to this list in the public analysis protocol.
additional_chemical_ids_to_remove <- c(
  100001604, 100020837, 501, 100001002, 1342, 100006360, 100006361
)

read_metabolon_csv <- function(path) {
  x <- read.csv(path, header = TRUE, fileEncoding = "UTF-8-BOM",
                check.names = TRUE)
  names(x)[1] <- "PARENT_SAMPLE_NAME"
  x
}

median_scale <- function(x) {
  denominator <- stats::median(x, na.rm = TRUE)
  if (!is.finite(denominator) || denominator == 0) {
    stop("A metabolite has a zero or non-finite median and cannot be scaled.")
  }
  x / denominator
}

# ---- Import and harmonise assay batches ------------------------------------
batch_2020 <- read_metabolon_csv(file.path(input_dir, "non_imputed_2020.csv"))
batch_2017 <- read_metabolon_csv(file.path(input_dir, "non_imputed_2017.csv"))
batch_2015 <- read_metabolon_csv(file.path(input_dir, "non_imputed_2015.csv"))

subjects <- read.csv(
  file.path(input_dir, "subjects.csv"),
  header = TRUE, fileEncoding = "UTF-8-BOM", check.names = FALSE
)
subjects$age <- as.numeric(subjects$age)

lookup <- read.delim(
  file.path(input_dir, "BNI_PEAK_CONC_METAB_LOOKUP.txt"),
  header = TRUE, check.names = FALSE
)
lookup$CHEMICAL_ID <- paste0("X", lookup$CHEM_ID)

required_subject_columns <- c("sample_id", "source", "group", "Batch", "age",
                              "Charmed.ID")
missing_subject_columns <- setdiff(required_subject_columns, names(subjects))
if (length(missing_subject_columns) > 0L) {
  stop("subjects.csv is missing: ", paste(missing_subject_columns, collapse = ", "))
}

# Retain only variables measured in all three batches. 
shared_columns <- Reduce(intersect, list(
  names(batch_2015), names(batch_2017), names(batch_2020)
))

shared_columns <- names(batch_2015)[names(batch_2015) %in% shared_columns]
if (!"PARENT_SAMPLE_NAME" %in% shared_columns) {
  stop("PARENT_SAMPLE_NAME is not present in every assay file.")
}

batch_2020 <- batch_2020[, shared_columns, drop = FALSE]
batch_2017 <- batch_2017[, shared_columns, drop = FALSE]
batch_2015 <- batch_2015[, shared_columns, drop = FALSE]

# Exclude 2015 samples without accompanying metadata.
batch_2015 <- batch_2015[
  batch_2015$PARENT_SAMPLE_NAME %in% subjects$sample_id, , drop = FALSE
]

# Remove missingness-filtering step.
raw_all_for_filtering <- bind_rows(batch_2015, batch_2017, batch_2020)
metadata_for_filtering <- subjects %>%
  select(sample_id, source, group, Batch, age)
combined_for_filtering <- raw_all_for_filtering %>%
  left_join(metadata_for_filtering, by = c("PARENT_SAMPLE_NAME" = "sample_id"))

# Metabolites missing in >30% of control samples are excluded. Missingness is
# assessed in controls rather than the whole cohort to avoid diagnosis-driven
# selection of features.
control_data <- combined_for_filtering %>% filter(group == "control")
candidate_metabolites <- setdiff(shared_columns, "PARENT_SAMPLE_NAME")
control_missingness <- vapply(
  control_data[, candidate_metabolites, drop = FALSE],
  function(x) mean(is.na(x)), numeric(1)
)
missingness_exclusions <- names(control_missingness)[
  control_missingness > missingness_threshold
]

predefined_exclusions <- lookup$CHEMICAL_ID[
  lookup$CHEM_ID %in% additional_chemical_ids_to_remove
]
metabolite_names <- setdiff(
  candidate_metabolites,
  union(missingness_exclusions, predefined_exclusions)
)

if (length(metabolite_names) == 0L) stop("No metabolites remain after filtering.")
if (!all(vapply(raw_all_for_filtering[, metabolite_names, drop = FALSE], is.numeric,
                logical(1)))) {
  stop("All retained metabolite columns must be numeric.")
}


skew_input <- bind_rows(batch_2020, batch_2015, batch_2017) %>%
  select(PARENT_SAMPLE_NAME, all_of(metabolite_names))

# Match the original e1071::skewness() defaults and selection rule exactly.
skewness <- vapply(
  skew_input[, metabolite_names, drop = FALSE],
  function(x) e1071::skewness(x, na.rm = TRUE), numeric(1)
)
log_metabolites <- names(skewness)[
  is.finite(skewness) & abs(skewness) > skewness_threshold
]

if (length(log_metabolites) > 0L) {
  non_positive <- vapply(
    skew_input[, log_metabolites, drop = FALSE],
    function(x) any(x <= 0, na.rm = TRUE), logical(1)
  )
  if (any(non_positive)) {
    stop(
      "Log-selected metabolites contain non-positive values: ",
      paste(names(non_positive)[non_positive], collapse = ", "),
      ". Specify and document an appropriate offset before continuing."
    )
  }
}

transform_batch <- function(x) {
  x <- x[, c("PARENT_SAMPLE_NAME", metabolite_names), drop = FALSE]
  x[log_metabolites] <- lapply(x[log_metabolites], log)
  x
}

transformed_2015 <- transform_batch(batch_2015)
transformed_2017 <- transform_batch(batch_2017)
transformed_2020 <- transform_batch(batch_2020)
raw_filtered <- bind_rows(transformed_2015, transformed_2017, transformed_2020)

# Retain untransformed values separately for downstream fold-change estimates.
raw_all <- bind_rows(batch_2020, batch_2015, batch_2017) %>%
  select(PARENT_SAMPLE_NAME, all_of(metabolite_names))

# Random-forest imputation is performed jointly across assay batches. No
# clinical variables are supplied to the imputer.
imputation <- missForest::missForest(
  raw_filtered[, metabolite_names, drop = FALSE],
  verbose = TRUE
)
imputed_metabolites <- as.data.frame(imputation$ximp, check.names = FALSE)
imputed <- bind_cols(
  data.frame(PARENT_SAMPLE_NAME = raw_filtered$PARENT_SAMPLE_NAME),
  imputed_metabolites
)

analysis_metadata <- subjects %>%
  select(sample_id, source, group, Batch, age, Charmed.ID)
data_with_metadata <- imputed %>%
  left_join(analysis_metadata, by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
  mutate(id = if_else(is.na(Charmed.ID), PARENT_SAMPLE_NAME,
                      as.character(Charmed.ID)))

# Median scaling is done after imputation and before RUV-III correction.
scaled_matrix <- data_with_metadata[, metabolite_names, drop = FALSE] %>%
  mutate(across(everything(), median_scale)) %>%
  as.matrix()

replicate_id <- paste(data_with_metadata$id, data_with_metadata$age, sep = "_")
replicate_counts <- table(replicate_id)
if (!any(replicate_counts > 1L)) {
  stop("RUV-III requires at least one technical replicate set.")
}

replicate_design <- ruv::replicate.matrix(replicate_id)
ruv_fit <- ruv::RUVIII(Y = scaled_matrix, M = replicate_design, k = ruv_k)

corrected <- as.data.frame(ruv_fit, check.names = FALSE)
#names(corrected) <- metabolite_names
corrected <- bind_cols(
data.frame(
PARENT_SAMPLE_NAME = data_with_metadata$PARENT_SAMPLE_NAME,
replicate_id = replicate_id
),
corrected
)

#  Retain the first record from each technical replicate set for analysis.
corrected_unique <- corrected %>%
  filter(!duplicated(replicate_id)) %>%
  select(-replicate_id)

processing_parameters <- list(
  random_seed = 1234,
  missingness_threshold = missingness_threshold,
  skewness_threshold = skewness_threshold,
  ruv_k = ruv_k,
  predefined_chemical_id_exclusions = additional_chemical_ids_to_remove,
  log_transformed_metabolites = log_metabolites,
  missForest_oob_error = imputation$OOBerror,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string
)

output <- list(
  corrected_metabolites = corrected_unique,
  raw_metabolites = raw_all,
  subjects = subjects,
  lookup = lookup,
  metabolite_names = metabolite_names,
  exclusions = list(
    missingness = missingness_exclusions,
    predefined = predefined_exclusions
  ),
  parameters = processing_parameters
)

output_file <- file.path(output_dir, "metabolomics_preprocessed.rds")
saveRDS(output, output_file)
message("Saved preprocessing output to: ", output_file)
