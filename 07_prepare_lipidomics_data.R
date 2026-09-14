#!/usr/bin/env Rscript

# Lipidomics study: import, transform and impute the analysis datasets
#


required_packages <- c("dplyr", "e1071", "missForest", "psych")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the following packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(dplyr))
options(stringsAsFactors = FALSE)

# ---- Paths -----------------------------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("LIPIDOMICS_PROJECT_DIR", unset = "."), mustWork = FALSE
)
input_dir <- Sys.getenv(
  "LIPIDOMICS_INPUT_DIR", unset = file.path(project_dir, "data", "raw")
)
processed_dir <- Sys.getenv(
  "LIPIDOMICS_OUTPUT_DIR", unset = file.path(project_dir, "data", "processed")
)
results_dir <- file.path(project_dir, "results", "lipidomics_preparation")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  subjects = "subjects.csv",
  lipid_classes = "lipid_class_all.csv",
  lipid_compositions = "lipid_all_composition.csv",
  lipid_species = "lipid_species_compositions.csv"
)
input_paths <- file.path(input_dir, input_files)
missing_inputs <- input_files[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing input files in ", input_dir, ": ",
       paste(missing_inputs, collapse = ", "))
}

subjects <- read.csv(input_paths[1], fileEncoding = "UTF-8-BOM")
lipids_all <- read.csv(
  input_paths[2], fileEncoding = "UTF-8-BOM"
)
lipid_compositions <- read.csv(
  input_paths[3], fileEncoding = "UTF-8-BOM"
)
lipid_species <- read.csv(
  input_paths[4], fileEncoding = "UTF-8-BOM"
)

subjects$age <- as.numeric(subjects$age)

# left_join() retains the row order of the lipidomics input, as in the source.
lipids_all <- left_join(
  lipids_all, subjects, by = c("PARENT_SAMPLE_NAME" = "sample_id")
)
lipid_compositions <- left_join(
  lipid_compositions, subjects, by = c("PARENT_SAMPLE_NAME" = "sample_id")
)
lipid_species <- left_join(
  lipid_species, subjects, by = c("PARENT_SAMPLE_NAME" = "sample_id")
)

# The source defines the 14 lipid-class variables by their input positions.
lipid_class_names <- names(lipids_all)[2:15]
composition_names <- names(lipid_compositions)[2:15]
lipids_all$total <- rowSums(lipids_all[, lipid_class_names, drop = FALSE])



log_highly_skewed <- function(data, variables, threshold = 2) {
  skews <- vapply(data[, variables, drop = FALSE], function(x) {
    if (is.numeric(x)) e1071::skewness(x, na.rm = TRUE) else NA_real_
  }, numeric(1))
  transformed <- names(skews)[!is.na(skews) & abs(skews) > threshold]
  for (variable in transformed) data[[variable]] <- log(data[[variable]])
  list(data = data, transformed = transformed, skewness = skews)
}

class_transform <- log_highly_skewed(lipids_all, lipid_class_names)
lipids_all <- class_transform$data


composition_transform <- log_highly_skewed(
  lipid_compositions, composition_names
)
lipid_compositions <- composition_transform$data

## removes variables with >30% missingness. Positional bounds are retained so
# that the same input produces the same set of variables.

healthy <- subset(lipid_species, group == "control")
healthy_data <- healthy[, 2:1002, drop = FALSE]
columns_to_remove <- names(healthy_data)[
  colSums(is.na(healthy_data)) / nrow(healthy_data) > 0.30
]
lipid_species <- lipid_species[
  , !(names(lipid_species) %in% columns_to_remove), drop = FALSE
]


lipid_species_names <- names(lipid_species)[2:875]
species_transform <- log_highly_skewed(
  lipid_species, lipid_species_names
)
lipid_species <- species_transform$data

# missForest for inputation
lipid_species_cut <- lipid_species[, 2:875, drop = FALSE]
set.seed(1234)
imputation <- missForest::missForest(lipid_species_cut)
lipid_species_imputed <- cbind(
  lipid_species[, -c(2:875), drop = FALSE], imputation$ximp
)

# The PCA is fitted to the complete, transformed lipid-class dataset before
# either comparison cohort is selected.
newpca <- psych::principal(lipids_all[, 2:15], nfactors = 3)
scores <- newpca$scores
lipids_all <- cbind(lipids_all, scores)
pca_score_names <- colnames(scores)

output <- list(
  lipids_all = lipids_all,
  lipid_compositions = lipid_compositions,
  lipid_species_imputed = lipid_species_imputed,
  lipid_class_names = lipid_class_names,
  composition_names = composition_names,
  lipid_species_names = lipid_species_names,
  pca_score_names = pca_score_names,
  pca = newpca,
  imputation = list(
    OOBerror = imputation$OOBerror,
    NRMSE = imputation$NRMSE,
    PFC = imputation$PFC
  ),
  preprocessing = list(
    species_removed_for_control_missingness = columns_to_remove,
    class_log_transformed = class_transform$transformed,
    composition_log_transformed = composition_transform$transformed,
    species_log_transformed = species_transform$transformed,
    skewness_threshold = 2,
    control_missingness_threshold = 0.30,
    imputation_seed = 1234
  )
)
saveRDS(output, file.path(processed_dir, "lipidomics_preprocessed.rds"))

