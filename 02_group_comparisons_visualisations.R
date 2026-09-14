#!/usr/bin/env Rscript

# Metabolomics study: prespecified group comparisons and visualisations
#
# Run `01_preprocess_impute.R` first. 

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "lmerTest", "GSVA"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the following packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(lmerTest)
  library(GSVA)
})

set.seed(123)

# ---- Paths ------------------------------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("METABOLOMICS_PROJECT_DIR", unset = "."),
  mustWork = FALSE
)
processed_dir <- Sys.getenv(
  "METABOLOMICS_OUTPUT_DIR",
  unset = file.path(project_dir, "data", "processed")
)
results_dir <- Sys.getenv(
  "METABOLOMICS_RESULTS_DIR",
  unset = file.path(project_dir, "results")
)
table_dir <- file.path(results_dir, "tables")
figure_dir <- file.path(results_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input <- readRDS(file.path(processed_dir, "metabolomics_preprocessed.rds"))
metabolite_names <- input$metabolite_names
lookup <- input$lookup

group_data <- input$corrected_metabolites %>%
  left_join(input$subjects, by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
  mutate(
    id = if_else(source == "genfi", as.character(subject_id),
                 as.character(Charmed.ID)),
    id = if_else(is.na(id) | id == "", PARENT_SAMPLE_NAME, id),
    sex = recode(as.character(sex), Male = "m", Female = "f"),
    batch = if_else(source == "genfi", "genfi", as.character(Batch))
  )

# ---- Reusable modelling functions ------------------------------------------
first_observation <- function(data) {
  # Visit is preferred when present; age is the fallback ordering variable.
  ordering_variable <- if ("Visit" %in% names(data)) "Visit" else "age"
  data %>%
    arrange(id, .data[[ordering_variable]]) %>%
    group_by(id) %>%
    slice_head(n = 1L) %>%
    ungroup()
}

fit_one_outcome <- function(data, outcome, group_variable, reference,
                            covariates, random_intercept = NULL,
                            scale_outcome = TRUE) {
  model_data <- data
  model_data[[group_variable]] <- stats::relevel(
    droplevels(factor(model_data[[group_variable]])), ref = reference
  )
  if (nlevels(model_data[[group_variable]]) != 2L) {
    stop(group_variable, " must contain exactly two levels after filtering.")
  }
  model_data$.outcome <- model_data[[outcome]]
  if (scale_outcome) model_data$.outcome <- as.numeric(scale(model_data$.outcome))

  fixed_terms <- c(group_variable, covariates)
  rhs <- paste(fixed_terms, collapse = " + ")
  if (!is.null(random_intercept)) {
    rhs <- paste0(rhs, " + (1 | ", random_intercept, ")")
  }
  formula <- stats::as.formula(paste(".outcome ~", rhs))

  fit <- if (is.null(random_intercept)) {
    stats::lm(formula, data = model_data)
  } else {
    lmerTest::lmer(formula, data = model_data)
  }

  coefficients <- summary(fit)$coefficients
  term_row <- grep(paste0("^", group_variable), rownames(coefficients))
  if (length(term_row) != 1L) {
    stop("Could not identify one coefficient for ", group_variable,
         " in the model for ", outcome, ".")
  }

  p_column <- grep("^Pr\\(", colnames(coefficients), value = TRUE)
  c(
    estimate = coefficients[term_row, "Estimate"],
    t_value = coefficients[term_row, "t value"],
    p_value = coefficients[term_row, p_column]
  )
}

fit_metabolites <- function(data, group_variable, reference, covariates,
                            random_intercept = NULL, scale_outcome = TRUE) {
  result_matrix <- vapply(
    metabolite_names,
    function(metabolite) fit_one_outcome(
      data = data,
      outcome = metabolite,
      group_variable = group_variable,
      reference = reference,
      covariates = covariates,
      random_intercept = random_intercept,
      scale_outcome = scale_outcome
    ),
    numeric(3)
  )

  tibble(
    variable = metabolite_names,
    estimate = result_matrix["estimate", ],
    t_value = result_matrix["t_value", ],
    p_value = result_matrix["p_value", ],
    p_fdr = p.adjust(p_value, method = "fdr")
  ) %>%
    left_join(lookup, by = c("variable" = "CHEMICAL_ID"))
}

fit_pathways <- function(data, metabolite_results, group_variable, reference,
                         covariates, random_intercept = NULL) {
  pathway_sets <- split(
    metabolite_results$variable,
    metabolite_results$SUB_PATHWAY
  )
  pathway_sets <- pathway_sets[!is.na(names(pathway_sets)) & names(pathway_sets) != ""]

  expression_matrix <- t(scale(data[, metabolite_names, drop = FALSE]))
  colnames(expression_matrix) <- data$PARENT_SAMPLE_NAME
  rownames(expression_matrix) <- metabolite_names

  parameter <- GSVA::ssgseaParam(
    exprData = expression_matrix,
    geneSets = pathway_sets,
    minSize = 3,
    maxSize = 500
  )
  score_matrix <- GSVA::gsva(parameter, verbose = FALSE)
  scores <- as.data.frame(t(score_matrix), check.names = FALSE)
  pathway_names <- names(scores)
  scores$PARENT_SAMPLE_NAME <- rownames(scores)
  model_data <- data %>% left_join(scores, by = "PARENT_SAMPLE_NAME")

  result_matrix <- vapply(
    pathway_names,
    function(pathway) fit_one_outcome(
      data = model_data,
      outcome = pathway,
      group_variable = group_variable,
      reference = reference,
      covariates = covariates,
      random_intercept = random_intercept,
      scale_outcome = TRUE
    ),
    numeric(3)
  )

  tibble(
    pathway = pathway_names,
    estimate = result_matrix["estimate", ],
    t_value = result_matrix["t_value", ],
    p_value = result_matrix["p_value", ],
    p_fdr = p.adjust(p_value, method = "fdr")
  )
}

calculate_fold_changes <- function(analysis_data, group_variable, reference,
                                   comparison_level, excluded_batches = NULL,
                                   exclude_genfi = FALSE) {
  raw <- input$raw_metabolites %>%
    filter(PARENT_SAMPLE_NAME %in% analysis_data$PARENT_SAMPLE_NAME) %>%
    left_join(input$subjects, by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
    mutate(
      id = if_else(source == "genfi", as.character(subject_id),
                   as.character(Charmed.ID)),
      id = if_else(is.na(id) | id == "", PARENT_SAMPLE_NAME, id),
      analysis_batch = if_else(source == "genfi", "genfi", as.character(Batch))
    ) %>%
    group_by(id) %>%
    slice_max(age, n = 1L, with_ties = FALSE) %>%
    ungroup()

  group_map <- analysis_data %>%
    distinct(PARENT_SAMPLE_NAME, .keep_all = TRUE) %>%
    select(PARENT_SAMPLE_NAME, all_of(group_variable))
  # The subject table may already contain a column with this name (for example,
  # group_short). Remove it so the analysis-defined version is not suffixed
  # `.x`/`.y` by the join.
  if (group_variable %in% names(raw)) {
    raw <- raw %>% select(-all_of(group_variable))
  }
  raw <- raw %>% left_join(group_map, by = "PARENT_SAMPLE_NAME")
  if (length(excluded_batches) > 0L) {
    raw <- raw %>% filter(!Batch %in% excluded_batches)
  }
  if (exclude_genfi) raw <- raw %>% filter(source != "genfi")

  long_summary <- raw %>%
    filter(.data[[group_variable]] %in% c(reference, comparison_level)) %>%
    group_by(Batch, .data[[group_variable]]) %>%
    summarise(across(all_of(metabolite_names), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")

  batches <- unique(long_summary$Batch)
  batch_ratios <- lapply(batches, function(current_batch) {
    x <- long_summary %>% filter(Batch == current_batch)
    if (!all(c(reference, comparison_level) %in% x[[group_variable]])) return(NULL)
    denominator <- x[x[[group_variable]] == reference, metabolite_names,
                     drop = FALSE]
    numerator <- x[x[[group_variable]] == comparison_level, metabolite_names,
                   drop = FALSE]
    as.numeric(numerator[1, ] / denominator[1, ])
  })
  batch_ratios <- batch_ratios[!vapply(batch_ratios, is.null, logical(1))]
  if (length(batch_ratios) == 0L) {
    warning("No batch contained both contrast groups; fold changes are NA.")
    return(tibble(variable = metabolite_names, fold_change = NA_real_,
                  log2_fold_change = NA_real_))
  }

  ratio_matrix <- do.call(rbind, batch_ratios)
  fold_change <- colMeans(ratio_matrix, na.rm = TRUE)
  tibble(
    variable = metabolite_names,
    fold_change = fold_change,
    log2_fold_change = log2(fold_change)
  )
}

save_volcano <- function(results, title, filename) {
  plot_data <- results %>% mutate(neg_log10_fdr = -log10(p_fdr))
  p <- ggplot(plot_data,
              aes(x = log2_fold_change, y = neg_log10_fdr,
                  colour = SUPER_PATHWAY)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text_repel(
      data = plot_data %>% filter(p_fdr < 0.01),
      aes(label = CHEMICAL_NAME), size = 3, show.legend = FALSE,
      max.overlaps = Inf
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted",
               colour = "grey50") +
    labs(title = title, x = "log2 fold change",
         y = "-log10 FDR-adjusted p value", colour = "Super-pathway") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")
  ggsave(file.path(figure_dir, filename), p, width = 8, height = 6, dpi = 300)
}

run_comparison <- function(name, data, group_variable, reference, covariates,
                           random_intercept = NULL, fold_change = NULL) {
  data[[group_variable]] <- droplevels(factor(data[[group_variable]]))
  message("Running ", name, ": ", paste(levels(data[[group_variable]]),
                                         collapse = " versus "))
  metabolite_results <- fit_metabolites(
    data, group_variable, reference, covariates, random_intercept
  )
  if (!is.null(fold_change)) {
    fc <- calculate_fold_changes(
      data, group_variable, reference,
      comparison_level = fold_change$comparison_level,
      excluded_batches = fold_change$excluded_batches,
      exclude_genfi = fold_change$exclude_genfi
    )
    metabolite_results <- metabolite_results %>% left_join(fc, by = "variable")
  }
  pathway_results <- fit_pathways(
    data, metabolite_results, group_variable, reference, covariates,
    random_intercept
  )
  write.csv(metabolite_results,
            file.path(table_dir, paste0(name, "_metabolites.csv")),
            row.names = FALSE)
  write.csv(pathway_results,
            file.path(table_dir, paste0(name, "_pathways.csv")),
            row.names = FALSE)
  list(metabolites = metabolite_results, pathways = pathway_results)
}

# ---- Diagnostic groups used in the comparisons -----------------------------
# `all_patient_labels` defines the symptomatic FTLD cohort. 
all_patient_labels <- c("PSP", "bvftd", "GEN_2", "cbs", "nfvPPA", "ftld_mimic")

standard_fold_change <- list(
  excluded_batches = character(),
  exclude_genfi = FALSE
)

# ---- Comparison 1: all symptomatic patients versus controls ----------------
all_patients_data <- group_data %>%
  filter(group_short %in% c("control", "GEN_0", all_patient_labels)) %>%
  mutate(comparison_group = if_else(group_short == "control",
                                    "control", "patient"))

all_patients <- run_comparison(
  name = "all_patients_vs_controls",
  data = all_patients_data,
  group_variable = "group",
  reference = "control",
  covariates = c("age", "sex", "batch"),
  random_intercept = "id",
  fold_change = c(
    list(comparison_level = "patient"),
    standard_fold_change
  )
)
save_volcano(
  all_patients$metabolites,
  title = "All symptomatic FTLD patients versus controls",
  filename = "all_patients_vs_controls_volcano.png"
)

# ---- Comparison 2: PSP versus controls -------------------------------------
psp_control_data <- group_data %>%
  filter(group_short %in% c("control", "GEN_0", "PSP")) %>%
  mutate(comparison_group = if_else(group_short == "PSP",
                                    "PSP", "control"))

psp_control <- run_comparison(
  name = "psp_vs_controls",
  data = psp_control_data,
  group_variable = "comparison_group",
  reference = "control",
  covariates = c("age", "sex", "batch"),
  random_intercept = "id",
  fold_change = c(
    list(comparison_level = "PSP"),
    standard_fold_change
  )
)
save_volcano(
  psp_control$metabolites,
  title = "PSP versus controls",
  filename = "psp_vs_controls_volcano.png"
)

# ---- Comparison 3: bvFTD (sporadic bvFTD and GEN_2) versus controls --------
bvftd_control_data <- group_data %>%
  filter(group_short %in% c("control", "bvftd", "GEN_2", "GEN_0")) %>%
  mutate(comparison_group = if_else(group_short == "control",
                                    "control", "bvftd"))

standard_fold_change <- list(
  excluded_batches = c("b_fifteen"),
  exclude_genfi = FALSE
)

bvftd_control <- run_comparison(
  name = "bvftd_vs_controls",
  data = bvftd_control_data,
  group_variable = "group",
  reference = "control",
  covariates = c("age", "sex", "batch"),
  random_intercept = "id",
  fold_change = c(
    list(comparison_level = "patient"),
    standard_fold_change
  )
)
save_volcano(
  bvftd_control$metabolites,
  title = "bvFTD (including GEN_2) versus controls",
  filename = "bvftd_vs_controls_volcano.png"
)

# ---- Comparison 4: PSP versus bvFTD ----------------------------------------
# GEN_2 participants are combined with the sporadic bvFTD group.
psp_bvftd_data <- group_data %>%
  filter(group_short %in% c("PSP", "bvftd", "GEN_2")) %>%
  mutate(comparison_group = if_else(group_short == "PSP", "PSP", "bvftd"))

psp_bvftd <- run_comparison(
  name = "psp_vs_bvftd",
  data = psp_bvftd_data,
  group_variable = "comparison_group",
  reference = "bvftd",
  covariates = c("age", "sex", "batch"),
  random_intercept = "id",
  fold_change = c(
    list(comparison_level = "PSP"),
    standard_fold_change
  )
)
save_volcano(
  psp_bvftd$metabolites,
  title = "PSP versus bvFTD (including GEN_2)",
  filename = "psp_vs_bvftd_volcano.png"
)

# ---- Comparison 5: presymptomatic carriers versus non-carriers -------------
# Only the earliest observation per participant is retained, avoiding repeated
# observations in this cross-sectional presymptomatic comparison. GEN_0 is the
# reference, so positive estimates indicate higher values in GEN_1.
presymptomatic_data <- group_data %>%
  filter(group_short %in% c("GEN_0", "GEN_1")) %>%
  mutate(comparison_group = if_else(group_short == "GEN_1",
                                    "GEN_1", "GEN_0"))

presymptomatic_data %>% group_by(subject_id) %>% slice_max(Visit) -> presymptomatic_data

presymptomatic <- run_comparison(
  name = "presymptomatic_vs_noncarriers",
  data = presymptomatic_data,
  group_variable = "comparison_group",
  reference = "GEN_0",
  covariates = c("age", "sex"),
  fold_change = c(
    list(comparison_level = "GEN_1"),
    standard_fold_change
  )
)
save_volcano(
  presymptomatic$metabolites,
  title = "Presymptomatic mutation carriers versus non-carriers",
  filename = "presymptomatic_vs_noncarriers_volcano.png"
)

# ---- Combined pathway visualisation ----------------------------------------
all_pathways <- bind_rows(
  all_patients$pathways %>% mutate(comparison = "All patients vs controls"),
  psp_control$pathways %>% mutate(comparison = "PSP vs controls"),
  bvftd_control$pathways %>% mutate(comparison = "bvFTD vs controls"),
  psp_bvftd$pathways %>% mutate(comparison = "PSP vs bvFTD"),
  presymptomatic$pathways %>%
    mutate(comparison = "Presymptomatic vs non-carriers")
) %>%
  mutate(significance = case_when(
    p_fdr < 0.001 ~ "***",
    p_fdr < 0.01 ~ "**",
    p_fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

pathway_plot <- ggplot(
  all_pathways,
  aes(x = t_value, y = reorder(pathway, t_value), colour = t_value)
) +
  geom_point(size = 2.5) +
  geom_text(aes(label = significance), nudge_x = 0.35, colour = "black",
            size = 4) +
  facet_wrap(~comparison, scales = "free_y", ncol = 2) +
  labs(x = "t statistic", y = "Pathway") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(figure_dir, "all_comparisons_pathway_scores.png"),
       pathway_plot, width = 12, height = 10, dpi = 300)
write.csv(all_pathways, file.path(table_dir, "all_pathway_results.csv"),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
