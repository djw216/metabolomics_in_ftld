#!/usr/bin/env Rscript

# Metabolomics study: WGCNA module associations and visualisations
#
# Prerequisites:
#   01_preprocess_impute.R
#   03_wgcna_network.R
#


required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "lmerTest", "emmeans",
  "e1071", "psych", "MASS", "car", "ggsignif"
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
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(lmerTest)
  library(emmeans)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# ---- Paths ------------------------------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("METABOLOMICS_PROJECT_DIR", unset = "."), mustWork = FALSE
)
input_dir <- Sys.getenv(
  "METABOLOMICS_INPUT_DIR",
  unset = file.path(project_dir, "data", "raw")
)
processed_dir <- Sys.getenv(
  "METABOLOMICS_OUTPUT_DIR",
  unset = file.path(project_dir, "data", "processed")
)
results_dir <- file.path(project_dir, "results", "wgcna_associations")
table_dir <- file.path(results_dir, "tables")
figure_dir <- file.path(results_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

preprocessed <- readRDS(
  file.path(processed_dir, "metabolomics_preprocessed.rds")
)
network <- readRDS(file.path(processed_dir, "wgcna_network.rds"))

module_names <- setdiff(names(network$eigengenes), "PARENT_SAMPLE_NAME")
if (length(module_names) == 0L) stop("No module eigengenes were found.")

group_data <- preprocessed$corrected_metabolites %>%
  left_join(network$eigengenes, by = "PARENT_SAMPLE_NAME") %>%
  left_join(preprocessed$subjects,
            by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
  mutate(
    id = if_else(source == "genfi", as.character(subject_id),
                 as.character(Charmed.ID)),
    id = if_else(is.na(id) | id == "", PARENT_SAMPLE_NAME, id),
    Charmed.ID = as.character(Charmed.ID),
    sex = recode(as.character(sex), Male = "m", Female = "f"),
    batch = if_else(source == "genfi", "genfi", as.character(Batch))
  )

# ---- General modelling helpers ---------------------------------------------
first_observation <- function(data) {
  ordering_variable <- if ("Visit" %in% names(data)) "Visit" else "age"
  data %>%
    arrange(id, .data[[ordering_variable]]) %>%
    group_by(id) %>%
    slice_head(n = 1L) %>%
    ungroup()
}

last_observation <- function(data) {
  ordering_variable <- if ("Visit" %in% names(data)) "Visit" else "age"
  data %>%
    arrange(id, .data[[ordering_variable]]) %>%
    group_by(id) %>%
    slice_tail(n = 1L) %>%
    ungroup()
}

coefficient_row <- function(fit, pattern) {
  coefs <- summary(fit)$coefficients
  row <- grep(pattern, rownames(coefs))
  if (length(row) != 1L) return(c(estimate = NA, t_value = NA, p_value = NA))
  p_col <- grep("^Pr\\(", colnames(coefs), value = TRUE)
  c(estimate = coefs[row, "Estimate"],
    t_value = coefs[row, "t value"],
    p_value = coefs[row, p_col])
}

fit_module_predictor <- function(data, outcome, covariates,
                                 random_intercept = NULL,
                                 predictor_scaled = TRUE,
                                 outcome_scaled = TRUE) {
  rows <- lapply(module_names, function(module) {
    model_data <- data
    model_data$.outcome <- model_data[[outcome]]
    model_data$.module <- model_data[[module]]
    if (outcome_scaled) model_data$.outcome <- as.numeric(scale(model_data$.outcome))
    if (predictor_scaled) model_data$.module <- as.numeric(scale(model_data$.module))
    rhs <- paste(c(".module", covariates), collapse = " + ")
    if (!is.null(random_intercept)) {
      rhs <- paste0(rhs, " + (1 | ", random_intercept, ")")
    }
    formula <- as.formula(paste(".outcome ~", rhs))
    values <- tryCatch({
      fit <- if (is.null(random_intercept)) lm(formula, model_data) else
        lmerTest::lmer(formula, model_data, REML = FALSE)
      coefficient_row(fit, "^\\.module$")
    }, error = function(e) c(estimate = NA, t_value = NA, p_value = NA))
    data.frame(module = module, t(values), check.names = FALSE)
  })
  bind_rows(rows) %>%
    mutate(across(c(estimate, t_value, p_value), as.numeric),
           p_fdr = p.adjust(p_value, method = "fdr"))
}

# Fits each module eigengene as the outcome and one continuous variable as the
# predictor. This matches the disease-stage model in the original notebook.
fit_predictor_across_modules <- function(data, predictor, covariates,
                                         random_intercept = NULL) {
  rows <- lapply(module_names, function(module) {
    model_data <- data
    model_data$.module <- as.numeric(scale(model_data[[module]]))
    model_data$.predictor <- as.numeric(scale(model_data[[predictor]]))
    rhs <- paste(c(".predictor", covariates), collapse = " + ")
    if (!is.null(random_intercept)) {
      rhs <- paste0(rhs, " + (1 | ", random_intercept, ")")
    }
    formula <- as.formula(paste(".module ~", rhs))
    values <- tryCatch({
      fit <- if (is.null(random_intercept)) lm(formula, model_data) else
        lmerTest::lmer(formula, model_data, REML = FALSE)
      coefficient_row(fit, "^\\.predictor$")
    }, error = function(e) c(estimate = NA, t_value = NA, p_value = NA))
    data.frame(module = module, t(values), check.names = FALSE)
  })
  bind_rows(rows) %>%
    mutate(across(c(estimate, t_value, p_value), as.numeric),
           p_fdr = p.adjust(p_value, method = "fdr"))
}

fit_metabolite_predictors <- function(data, outcome, covariates) {
  metabolite_names <- network$metabolite_names
  rows <- lapply(metabolite_names, function(metabolite) {
    model_data <- data
    model_data$.outcome <- model_data[[outcome]]
    model_data$.metabolite <- model_data[[metabolite]]
    formula <- as.formula(paste(
      ".outcome ~ .metabolite +", paste(covariates, collapse = " + ")
    ))
    values <- tryCatch({
      fit <- lm(formula, model_data)
      coefficient_row(fit, "^\\.metabolite$")
    }, error = function(e) c(estimate = NA, t_value = NA, p_value = NA))
    data.frame(CHEMICAL_ID = metabolite, t(values), check.names = FALSE)
  })
  bind_rows(rows) %>%
    mutate(across(c(estimate, t_value, p_value), as.numeric),
           p_fdr = p.adjust(p_value, method = "fdr")) %>%
    left_join(preprocessed$lookup, by = "CHEMICAL_ID")
}

save_metabolite_association_plot <- function(results, title, filename) {
  plot_data <- results %>% mutate(neg_log10_fdr = -log10(p_fdr))
  p <- ggplot(plot_data, aes(x = t_value, y = neg_log10_fdr,
                             colour = SUPER_PATHWAY)) +
    geom_point(alpha = 0.8) +
    geom_text_repel(
      data = plot_data %>% filter(p_fdr < 0.05),
      aes(label = CHEMICAL_NAME), size = 3, show.legend = FALSE,
      max.overlaps = Inf
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted",
               colour = "grey50") +
    labs(title = title, x = "t statistic",
         y = "-log10 FDR-adjusted p value", colour = "Super-pathway") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")
  ggsave(file.path(figure_dir, filename), p, width = 8, height = 6, dpi = 300)
}

save_diagnostic_group_plot <- function(data, module, omnibus_fdr,
                                       posthoc_results, filename) {
  display_labels <- c(
    control = "Control", bvftd = "bvFTD", cbs = "CBS",
    nfvPPA = "nfvPPA", PSP = "PSP"
  )
  significant <- posthoc_results %>%
    filter(.data$module == .env$module, p.value < 0.05) %>%
    mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        TRUE ~ "*"
      )
    )

  p <- ggplot(data, aes(x = group_short2, y = .data[[module]],
                        colour = group_short2)) +
    geom_boxplot(outlier.shape = NA, colour = "grey30") +
    geom_jitter(width = 0.12, alpha = 0.75) +
    scale_x_discrete(labels = display_labels) +
    labs(
      title = paste("Diagnostic-group differences:", module),
      subtitle = paste0("Omnibus FDR-adjusted P = ",
                        format.pval(omnibus_fdr, digits = 3)),
      x = NULL, y = module
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")

  if (nrow(significant) > 0L) {
    comparisons <- strsplit(as.character(significant$contrast), " - ",
                            fixed = TRUE)
    y_range <- range(data[[module]], na.rm = TRUE)
    y_span <- diff(y_range)
    if (!is.finite(y_span) || y_span == 0) y_span <- 1
    y_positions <- y_range[2] + y_span * (0.08 +
      0.08 * (seq_len(nrow(significant)) - 1))
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons,
      annotations = significant$stars,
      y_position = y_positions,
      tip_length = 0.02,
      textsize = 4,
      colour = "black"
    ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15 +
        0.08 * nrow(significant))))
  }
  ggsave(file.path(figure_dir, filename), p, width = 8, height = 6, dpi = 300)
}

save_strongest_module_plot <- function(data, results, x, x_label, title,
                                       filename) {
  valid <- results %>% filter(is.finite(p_value)) %>% arrange(p_value)
  if (nrow(valid) == 0L) return(invisible(NULL))
  module <- valid$module[1]
  p <- ggplot(data, aes(x = .data[[module]], y = .data[[x]])) +
    geom_point(aes(colour = group_short), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, colour = "black") +
    labs(title = title, x = module, y = x_label,
         subtitle = paste0("FDR-adjusted P = ",
                           format.pval(valid$p_fdr[1], digits = 3))) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir, filename), p, width = 7, height = 5, dpi = 300)
}

# ---- Omnibus diagnostic-group comparison -----------------------------------
diagnostic_group_data <- group_data %>%
  mutate(
    group_short2 = case_when(
      group_short == "GEN_2" ~ "bvftd",
      group_short == "GEN_0" ~ "control",
      tolower(group_short) == "psp" ~ "PSP",
      tolower(group_short) == "nfvppa" ~ "nfvPPA",
      tolower(group_short) == "cbs" ~ "cbs",
      TRUE ~ as.character(group_short)
    )
  ) %>%
  filter(group_short2 %in% c("control", "bvftd", "cbs", "nfvPPA", "PSP")) %>%
  mutate(group_short2 = factor(
    group_short2,
    levels = c("control", "bvftd", "cbs", "nfvPPA", "PSP")
  ))

diagnostic_group_counts <- diagnostic_group_data %>%
  count(group_short2, name = "n_samples", .drop = FALSE)
write.csv(diagnostic_group_counts,
          file.path(table_dir, "diagnostic_group_sample_counts.csv"),
          row.names = FALSE)
missing_diagnostic_groups <- diagnostic_group_counts$group_short2[
  diagnostic_group_counts$n_samples == 0
]
if (length(missing_diagnostic_groups) > 0L) {
  stop("The omnibus model is missing required diagnostic groups: ",
       paste(missing_diagnostic_groups, collapse = ", "))
}

diagnostic_fits <- lapply(module_names, function(module) {
  model_data <- diagnostic_group_data
  model_data$.module <- model_data[[module]]
  tryCatch(
    lmerTest::lmer(
      .module ~ group_short2 + age + sex + batch + (1 | id),
      data = model_data, REML = TRUE
    ),
    error = function(e) NULL
  )
})
names(diagnostic_fits) <- module_names

diagnostic_omnibus <- lapply(module_names, function(module) {
  fit <- diagnostic_fits[[module]]
  if (is.null(fit)) {
    return(data.frame(module = module, numerator_df = NA_real_,
                      denominator_df = NA_real_, f_value = NA_real_,
                      p_value = NA_real_))
  }
  test <- anova(fit)
  data.frame(
    module = module,
    numerator_df = test["group_short2", "NumDF"],
    denominator_df = test["group_short2", "DenDF"],
    f_value = test["group_short2", "F value"],
    p_value = test["group_short2", "Pr(>F)"]
  )
}) %>%
  bind_rows() %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))

diagnostic_posthoc <- lapply(module_names, function(module) {
  fit <- diagnostic_fits[[module]]
  if (is.null(fit)) return(NULL)
  as.data.frame(
    pairs(emmeans::emmeans(fit, ~group_short2), adjust = "tukey")
  ) %>%
    mutate(module = module, .before = 1)
}) %>%
  bind_rows() %>%
  left_join(
    diagnostic_omnibus %>% select(module, omnibus_p_value = p_value,
                                  omnibus_p_fdr = p_fdr),
    by = "module"
  )

write.csv(
  diagnostic_omnibus,
  file.path(table_dir, "diagnostic_group_module_omnibus.csv"),
  row.names = FALSE
)
write.csv(
  diagnostic_posthoc,
  file.path(table_dir, "diagnostic_group_module_posthoc_emmeans.csv"),
  row.names = FALSE
)



# ---- Pathology-group comparison --------------------------------------------
pathology_data <- group_data %>%
  filter(path_combined %in% c("tau", "tdp", "ftld_mimic"),
         group_short != "GEN_0") %>%
  mutate(path_combined = factor(path_combined,
                                levels = c("ftld_mimic", "tau", "tdp")))

pathology_omnibus <- lapply(module_names, function(module) {
  pathology_data$.module <- pathology_data[[module]]
  fit <- tryCatch(
    lmerTest::lmer(.module ~ path_combined + age + sex + batch + (1 | id),
                   pathology_data, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(module = module, p_value = NA_real_))
  }
  aov_table <- anova(fit)
  data.frame(module = module,
             p_value = aov_table["path_combined", "Pr(>F)"])
}) %>% bind_rows() %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))
write.csv(pathology_omnibus,
          file.path(table_dir, "pathology_group_omnibus.csv"), row.names = FALSE)

pathology_posthoc <- lapply(module_names, function(module) {
  pathology_data$.module <- pathology_data[[module]]
  fit <- tryCatch(
    lmerTest::lmer(.module ~ path_combined + age + sex + batch + (1 | id),
                   pathology_data, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  as.data.frame(pairs(emmeans::emmeans(fit, ~path_combined), adjust = "tukey")) %>%
    mutate(module = module)
}) %>% bind_rows() %>%
  mutate(p_fdr_across_all_module_contrasts = p.adjust(p.value, method = "fdr"))
write.csv(pathology_posthoc,
          file.path(table_dir, "pathology_group_posthoc.csv"), row.names = FALSE)

