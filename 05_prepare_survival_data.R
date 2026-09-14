#!/usr/bin/env Rscript

# Metabolomics study: prepare survival outcomes and counting-process intervals
#
# Prerequisites:
#   01_preprocess_impute.R
#   03_wgcna_network.R
#


required_packages <- c("dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the following packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(dplyr))

# ---- Paths and censoring date ----------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("METABOLOMICS_PROJECT_DIR", unset = "."), mustWork = FALSE
)
processed_dir <- Sys.getenv(
  "METABOLOMICS_OUTPUT_DIR",
  unset = file.path(project_dir, "data", "processed")
)
results_dir <- file.path(project_dir, "results", "survival_preparation")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# This reproduces the original study censoring date. Override using an ISO date,
# for example SURVIVAL_CENSOR_DATE=2026-09-01.
censor_date <- as.Date(
  Sys.getenv("SURVIVAL_CENSOR_DATE", unset = "2025-09-01")
)
if (is.na(censor_date)) stop("SURVIVAL_CENSOR_DATE must use YYYY-MM-DD format.")

preprocessed <- readRDS(
  file.path(processed_dir, "metabolomics_preprocessed.rds")
)
network <- readRDS(file.path(processed_dir, "wgcna_network.rds"))

metabolite_names <- preprocessed$metabolite_names
module_names <- setdiff(names(network$eigengenes), "PARENT_SAMPLE_NAME")

z_scale_frame <- function(data) {
  as.data.frame(lapply(data, function(x) as.numeric(scale(x))),
                check.names = FALSE)
}

# Standardised metabolites and eigengenes across
# the complete available cohort before selecting PSP or other patient subsets.
scaled_metabolites <- z_scale_frame(
  preprocessed$corrected_metabolites[, metabolite_names, drop = FALSE]
)
scaled_metabolites <- bind_cols(
  preprocessed$corrected_metabolites %>% select(PARENT_SAMPLE_NAME),
  scaled_metabolites
)

scaled_modules <- z_scale_frame(
  network$eigengenes[, module_names, drop = FALSE]
)
scaled_modules <- bind_cols(
  network$eigengenes %>% select(PARENT_SAMPLE_NAME),
  scaled_modules
)

required_subject_columns <- c(
  "sample_id", "source", "group_short", "status", "age", "sex", "Batch",
  "primary_date", "age_primary_date", "age_at_death"
)
absent <- setdiff(required_subject_columns, names(preprocessed$subjects))
if (length(absent) > 0L) {
  stop("subjects.csv is missing survival variables: ",
       paste(absent, collapse = ", "))
}

sample_data <- scaled_metabolites %>%
  left_join(scaled_modules, by = "PARENT_SAMPLE_NAME") %>%
  left_join(preprocessed$subjects,
            by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
  mutate(
    id = case_when(
      source == "genfi" & !is.na(subject_id) ~ as.character(subject_id),
      !is.na(Charmed.ID) ~ as.character(Charmed.ID),
      TRUE ~ PARENT_SAMPLE_NAME
    ),
    age = as.numeric(age),
    age_primary_date = as.numeric(age_primary_date),
    age_at_death = as.numeric(age_at_death),
    sex = case_when(
      sex %in% c("m", "M", "Male", "male") ~ "m",
      sex %in% c("f", "F", "Female", "female") ~ "f",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("f", "m")),
    batch = factor(if_else(source == "genfi", "genfi", as.character(Batch))),
    primary_date_parsed = as.Date(primary_date, format = "%d/%m/%Y"),
    age_at_censor = age_primary_date +
      as.numeric(censor_date - primary_date_parsed) / 365.25,
    final_followup_age = if_else(status == "Deceased",
                                 age_at_death, age_at_censor),
    followup_from_sample = final_followup_age - age
  ) %>%
  filter(status %in% c("Alive", "Deceased"))

# Fail rather than silently analyse impossible dates or survival intervals.
invalid_followup <- sample_data %>%
  filter(!is.finite(age) | !is.finite(final_followup_age) |
           followup_from_sample < 0)
if (nrow(invalid_followup) > 0L) {
  write.csv(invalid_followup,
            file.path(results_dir, "invalid_survival_records.csv"),
            row.names = FALSE)
  stop("Invalid survival records found; see invalid_survival_records.csv.")
}

# Counting-process format: each metabolomics measurement applies from its age
# until the next measurement, death, or administrative censoring. Only the last
# interval can contain the death event.
interval_data <- sample_data %>%
  arrange(id, age) %>%
  group_by(id) %>%
  mutate(
    age_at_entry = first(age),
    next_sample_age = lead(age),
    interval_end_age = if_else(!is.na(next_sample_age),
                               next_sample_age, final_followup_age),
    time1 = age - age_at_entry,
    time2 = interval_end_age - age_at_entry,
    event = as.integer(is.na(next_sample_age) & status == "Deceased")
  ) %>%
  ungroup()

invalid_intervals <- interval_data %>%
  filter(!is.finite(time1) | !is.finite(time2) | time1 < 0 | time2 <= time1 |
           interval_end_age > final_followup_age)
if (nrow(invalid_intervals) > 0L) {
  write.csv(invalid_intervals,
            file.path(results_dir, "invalid_counting_process_intervals.csv"),
            row.names = FALSE)
  stop("Invalid start-stop intervals found; see the preparation results folder.")
}

# Baseline data support conventional participant-level Kaplan–Meier curves and
# one-year prediction. This avoids treating repeated intervals as independent
# participants in those analyses.
baseline_data <- interval_data %>%
  group_by(id) %>%
  slice_min(age, n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    survival_from_entry = final_followup_age - age_at_entry,
    baseline_event = as.integer(status == "Deceased")
  )

cohort_summary <- baseline_data %>%
  summarise(
    n_participants = n_distinct(id),
    n_samples = nrow(interval_data),
    n_deaths = sum(baseline_event),
    median_followup_years = median(survival_from_entry, na.rm = TRUE),
    censor_date = as.character(censor_date)
  )
write.csv(cohort_summary, file.path(results_dir, "cohort_summary.csv"),
          row.names = FALSE)

output <- list(
  interval_data = interval_data,
  baseline_data = baseline_data,
  sample_data = sample_data,
  metabolite_names = metabolite_names,
  module_names = module_names,
  lookup = preprocessed$lookup,
  parameters = list(
    censor_date = as.character(censor_date),
    metabolite_scaling = "z score across complete corrected cohort",
    module_scaling = "z score across complete cohort",
    time_scale = "years since first metabolomics sample",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
)
saveRDS(output, file.path(processed_dir, "survival_prepared.rds"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
