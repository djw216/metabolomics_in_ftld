#!/usr/bin/env Rscript

# Metabolomics study: survival models, one-year prediction and visualisations
#
# Prerequisite: run 05_prepare_survival_data.R.


required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "survival", "survminer",
  "GSVA", "glmnet", "pROC"
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
  library(survival)
  library(survminer)
  library(GSVA)
  library(glmnet)
  library(pROC)
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
results_dir <- file.path(project_dir, "results", "survival_models")
table_dir <- file.path(results_dir, "tables")
figure_dir <- file.path(results_dir, "figures")
model_dir <- file.path(results_dir, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

survival_data <- readRDS(file.path(processed_dir, "survival_prepared.rds"))
interval_data <- survival_data$interval_data
baseline_data <- survival_data$baseline_data
metabolite_names <- survival_data$metabolite_names
module_names <- c("MEblack", "MEgreen", "MEmagenta", 
                             "MEyellow", "MEred", "MEturquoise")
lookup <- survival_data$lookup

# ---- Cox-model helpers ------------------------------------------------------
fit_one_cox_predictor <- function(data, predictor,
                                  covariates = c("age_at_entry", "sex")) {
  model_data <- data
  model_data$.predictor <- model_data[[predictor]]
  rhs <- paste(c(".predictor", covariates, "cluster(id)"), collapse = " + ")
  formula <- as.formula(paste(
    "Surv(time1, time2, event) ~", rhs
  ))
  tryCatch({
    fit <- coxph(formula, data = model_data, ties = "efron", x = TRUE)
    beta <- unname(coef(fit)[".predictor"])
    standard_error <- sqrt(vcov(fit)[".predictor", ".predictor"])
    z_value <- beta / standard_error
    p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
    data.frame(
      variable = predictor,
      log_hazard_ratio = beta,
      hazard_ratio = exp(beta),
      lower_95 = exp(beta - 1.96 * standard_error),
      upper_95 = exp(beta + 1.96 * standard_error),
      robust_standard_error = standard_error,
      z_value = z_value,
      p_value = p_value
    )
  }, error = function(e) {
    data.frame(
      variable = predictor,
      log_hazard_ratio = NA_real_, hazard_ratio = NA_real_,
      lower_95 = NA_real_, upper_95 = NA_real_,
      robust_standard_error = NA_real_, z_value = NA_real_,
      p_value = NA_real_
    )
  })
}

fit_cox_screen <- function(data, variables,
                           covariates = c("age_at_entry", "sex")) {
  bind_rows(lapply(variables, function(variable) {
    fit_one_cox_predictor(data, variable, covariates)
  })) %>%
    mutate(p_fdr = p.adjust(p_value, method = "fdr"))
}

save_cox_volcano <- function(results, title, filename) {
  plot_data <- results %>% mutate(neg_log10_fdr = -log10(p_fdr))
  p <- ggplot(plot_data,
              aes(x = log_hazard_ratio, y = neg_log10_fdr,
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
    labs(title = title, x = "Log hazard ratio per SD",
         y = "-log10 FDR-adjusted p value", colour = "Super-pathway") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")
  ggsave(file.path(figure_dir, filename), p, width = 8, height = 6, dpi = 300)
}

# ---- PSP metabolite-level Cox models ---------------------------------------
##This can be modified to obtain results for all patients
psp_intervals <- interval_data %>% filter(group_short == "PSP")
if (n_distinct(psp_intervals$id) < 10L || sum(psp_intervals$event) < 5L) {
  stop("Too few PSP participants or deaths for the planned survival models.")
}

metabolite_cox <- fit_cox_screen(psp_intervals, metabolite_names) %>%
  left_join(lookup, by = c("variable" = "CHEMICAL_ID"))
write.csv(metabolite_cox,
          file.path(table_dir, "psp_metabolite_cox_models.csv"),
          row.names = FALSE)
save_cox_volcano(
  metabolite_cox, "PSP metabolite associations with survival",
  "psp_metabolite_survival.png"
)

# ---- PSP super-pathway scores and Cox models -------------------------------
super_pathways <- split(
  metabolite_cox$variable,
  metabolite_cox$SUPER_PATHWAY
)
super_pathways <- super_pathways[
  !is.na(names(super_pathways)) & names(super_pathways) != ""
]
psp_intervals2 <- psp_intervals
psp_intervals2[,2:476] <- scale(psp_intervals[,2:476])
psp_expression <- t(as.matrix(psp_intervals2[, metabolite_names, drop = FALSE]))
colnames(psp_expression) <- psp_intervals$PARENT_SAMPLE_NAME
rownames(psp_expression) <- metabolite_names
pathway_parameter <- GSVA::ssgseaParam(
  exprData = psp_expression,
  geneSets = super_pathways,
  minSize = 3,
  maxSize = 500
)
pathway_scores <- as.data.frame(
  t(GSVA::gsva(pathway_parameter, verbose = FALSE)),
  check.names = FALSE
)
pathway_names <- names(pathway_scores)
pathway_scores$PARENT_SAMPLE_NAME <- rownames(pathway_scores)
psp_pathway_data <- psp_intervals %>%
  left_join(pathway_scores, by = "PARENT_SAMPLE_NAME")
psp_pathway_data[pathway_names] <- lapply(
  psp_pathway_data[pathway_names], function(x) as.numeric(scale(x))
)
pathway_cox <- fit_cox_screen(psp_pathway_data, pathway_names)
write.csv(pathway_cox,
          file.path(table_dir, "psp_super_pathway_cox_models.csv"),
          row.names = FALSE)

# ---- PSP module Cox models --------------------------------------------------
module_cox <- fit_cox_screen(psp_intervals, module_names)
write.csv(module_cox, file.path(table_dir, "psp_module_cox_models.csv"),
          row.names = FALSE)



# ---- One-year PSP survival prediction --------------------------------------
psprs_file <- file.path(input_dir, "psprs_scores.csv")

psp_scores <- read.csv(psprs_file, fileEncoding = "UTF-8-BOM") %>%
  select(sample_id, psp_42_tot, age_diff) %>%
  filter(abs(age_diff) < 1.00001)

prediction_data <- survival_data$sample_data %>%
  filter(group_short == "PSP") %>%
  left_join(psp_scores, by = c("PARENT_SAMPLE_NAME" = "sample_id")) %>%
  filter(!is.na(psp_42_tot)) %>%
  group_by(id) %>%
  slice_min(abs(age_diff), n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    event_1_year = case_when(
      status == "Deceased" & followup_from_sample < 1 ~ 1L,
      followup_from_sample >= 1 ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(event_1_year), !is.na(age), !is.na(sex))

ifelse(prediction_data$sex=="f", 0, 1) -> prediction_data$sex

feature_sets <- list(
  Baseline = c("age", "sex"),
  Baseline_PSPRS = c("age", "sex", "psp_42_tot"),
  Baseline_Metabolites = c("age", "sex", metabolite_names),
  Baseline_PSPRS_Metabolites =
    c("age", "sex", "psp_42_tot", metabolite_names)
)

df1 <- prediction_data
y <- df1$event_1_year

set.seed(123)

R <- 5

model_results <- list()
summary_all <- data.frame()

for (grp_name in names(feature_sets)) {
  
  cat("\nRunning model:", grp_name, "\n")
  
  cols <- feature_sets[[grp_name]]
  X <- as.matrix(df1[, cols])
  y <- df1$event_1_year
  
  # Retain participant and repetition identifiers
  results <- data.frame(
    participant_index = integer(),
    repetition = integer(),
    prob = numeric(),
    true = integer()
  )
  
  for (r in seq_len(R)) {
    
    cat(" Repetition:", r, "\n")
    set.seed(123 + r - 1)
    
    for (i in seq_len(nrow(X))) {
      
      test_idx <- i
      train_idx <- setdiff(seq_len(nrow(X)), i)
      
      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      
      X_test <- X[test_idx, , drop = FALSE]
      y_test <- y[test_idx]
      
      cv_fit <- cv.glmnet(
        x = X_train,
        y = y_train,
        family = "binomial",
        alpha = 1,
        nfolds = 5
      )
      
      model <- glmnet(
        x = X_train,
        y = y_train,
        family = "binomial",
        alpha = 1,
        lambda = cv_fit$lambda.min
      )
      
      prob_test <- as.numeric(
        predict(
          model,
          newx = X_test,
          type = "response"
        )
      )
      
      results <- rbind(
        results,
        data.frame(
          participant_index = i,
          repetition = r,
          prob = prob_test,
          true = y_test
        )
      )
    }
  }
  
  # One prediction per participant
  participant_results <- results %>%
    group_by(participant_index) %>%
    summarise(
      true = first(true),
      prob = mean(prob),
      prediction_sd = sd(prob),
      .groups = "drop"
    )
  
  # Calculate performance from N participants, not N × R predictions
  roc_obj <- pROC::roc(
    response = participant_results$true,
    predictor = participant_results$prob,
    levels = c(0, 1),
    direction = "<"
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  summary_all <- rbind(
    summary_all,
    data.frame(
      Model = grp_name,
      N = nrow(participant_results),
      AUC = auc_val
    )
  )
  
  model_results[[grp_name]] <- list(
    participant_results = participant_results,
    repeated_predictions = results,
    roc = roc_obj,
    auc = auc_val
  )
}

