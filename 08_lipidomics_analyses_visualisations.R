#!/usr/bin/env Rscript

# Lipidomics study: prespecified comparisons and visualisations
#
# Analyses:
#   1. all lipid species;
#   2. the first component of the lipid-class PCA;
#   3. lipid-class composition by PERMANOVA and component-wise models.
#
# Each analysis compares:
#   - all patients versus controls, excluding GEN_1;
#   - presymptomatic mutation carriers (GEN_1) versus non-carriers (GEN_0).
#
# Prerequisite: run 07_prepare_lipidomics_data.R.

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggtext", "compositions", "vegan"
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
})
options(stringsAsFactors = FALSE)

# ---- Paths -----------------------------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("LIPIDOMICS_PROJECT_DIR", unset = "."), mustWork = FALSE
)
processed_dir <- Sys.getenv(
  "LIPIDOMICS_OUTPUT_DIR", unset = file.path(project_dir, "data", "processed")
)
results_dir <- file.path(project_dir, "results", "lipidomics_analyses")
table_dir <- file.path(results_dir, "tables")
figure_dir <- file.path(results_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(file.path(processed_dir, "lipidomics_preprocessed.rds"))
lipid_species_imputed <- prepared$lipid_species_imputed
lipids_all <- prepared$lipids_all
lipid_compositions <- prepared$lipid_compositions
composition_names <- prepared$composition_names

required_metadata <- c(
  "group", "group_short", "age", "sex", "source", "subject_id",
  "BlindedCode", "Visit"
)
for (object_name in c("lipid_species_imputed", "lipids_all",
                      "lipid_compositions")) {
  absent <- setdiff(required_metadata, names(get(object_name)))
  if (length(absent) > 0L) {
    stop(object_name, " is missing metadata: ", paste(absent, collapse = ", "))
  }
}

extract_first_group_term <- function(fit) {
  coefficient_table <- summary(fit)$coefficients
  group_row <- grep("group", rownames(coefficient_table))[1]
  if (is.na(group_row)) stop("No group coefficient was estimable.")
  data.frame(
    estimate = coefficient_table[group_row, "Estimate"],
    tval = coefficient_table[group_row, "t value"],
    pval = coefficient_table[group_row, "Pr(>|t|)"]
  )
}

fit_species_models <- function(data, positions, covariates) {
  rows <- lapply(positions, function(position) {
    model_data <- data
    model_data$.outcome <- model_data[[position]]
    fit <- lm(
      reformulate(c("group", covariates), response = ".outcome"),
      data = model_data
    )
    cbind(variable = names(model_data)[position], extract_first_group_term(fit))
  })
  bind_rows(rows) %>%
    mutate(pval_fdr = p.adjust(pval, method = "fdr")) %>%
    select(variable, pval, tval, estimate, pval_fdr)
}

write_model_summary <- function(fit, filename) {
  coefficient_table <- as.data.frame(summary(fit)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  coefficient_table <- coefficient_table %>%
    select(term, everything())
  write.csv(coefficient_table, file.path(table_dir, filename), row.names = FALSE)
}

# ---- 1. Lipid-species comparisons -----------------------------------------
# 

gen_species_variables <- names(lipid_species_imputed)[57:930]
patient_species_variables <- names(lipid_species_imputed)[57:930]


species_gen <- lipid_species_imputed %>%
  filter(group_short %in% c("GEN_0", "GEN_1"))
species_gen_results <- fit_species_models(
  species_gen, 57:930, c("age", "sex")
)
write.csv(species_gen_results,
          file.path(table_dir, "lipid_species_GEN1_vs_GEN0.csv"),
          row.names = FALSE)

species_all <- lipid_species_imputed %>%
  filter(group_short != "GEN_1") %>%
  group_by(subject_id) %>%
  slice_max(Visit) %>%
  ungroup()
species_all_results <- fit_species_models(
  species_all, 57:930, c("age", "sex", "source")
)
write.csv(species_all_results,
          file.path(table_dir, "lipid_species_all_patients_vs_controls.csv"),
          row.names = FALSE)

# ---- 2. Lipid-class PCA ----------------------------------------------------
# The PCA itself was fitted in script 07 to the full transformed cohort, as in
# the source.
if (!"RC1" %in% names(lipids_all)) stop("The PCA did not produce an RC1 score.")

pca_gen <- lipids_all %>%
  filter(group_short %in% c("GEN_0", "GEN_1")) %>%
  group_by(BlindedCode) %>%
  slice_max(Visit) %>%
  ungroup()
pca_gen_fit <- lm(RC1 ~ group + age + sex, data = pca_gen)
write_model_summary(pca_gen_fit, "lipid_class_PCA_RC1_GEN1_vs_GEN0.csv")

# Plotting 
pca_gen_plot <- ggplot(pca_gen, aes(x = group, y = RC1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = group)) +
  theme_minimal(base_size = 18) +
  labs(y = "Lipid Class Component1") +
  theme(
    panel.grid = element_blank(), legend.title = element_blank(),
    legend.position = "none", axis.title.x = element_blank()
  ) +
  annotate(
    "richtext", label = "*t*=-2.1 FDR *P*=0.036",
    x = 1.5, y = 3.2, size = 6
  ) +
  scale_x_discrete(
    labels = c("control" = "Non-Carriers", "patient" = "Presymptomatic")
  )
ggsave(file.path(figure_dir, "lipid_class_PCA_RC1_GEN1_vs_GEN0.png"),
       pca_gen_plot, width = 7, height = 6, dpi = 300)

# All patients
pca_all <- lipids_all %>%
  filter(group_short != "GEN_1") %>%
  group_by(subject_id) %>%
  slice_max(Visit) %>%
  ungroup()
pca_all_fit <- lm(RC1 ~ group + age + sex + source, data = pca_all)
write_model_summary(
  pca_all_fit, "lipid_class_PCA_RC1_all_patients_vs_controls.csv"
)
pca_all_plot <- ggplot(pca_all, aes(x = group, y = RC1, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_minimal(base_size = 18) +
  labs(y = "Lipid Class Component1", x = NULL) +
  theme(panel.grid = element_blank(), legend.position = "none")
ggsave(
  file.path(figure_dir, "lipid_class_PCA_RC1_all_patients_vs_controls.png"),
  pca_all_plot, width = 7, height = 6, dpi = 300
)

# ---- 3. Lipid compositions: PERMANOVA and component-wise models -----------
run_composition_analysis <- function(data, group_variable,
                                     permanova_covariates, lm_covariates,
                                     group_coefficient, seed, file_stub,
                                     plot_title, plot_x) {
  message("Running ", file_stub, ": CLR transformation")
  clr_data <- compositions::clr(data[, composition_names, drop = FALSE])
  comp_matrix <- clr_data
  
  set.seed(seed)
  message("Running ", file_stub, ": PERMANOVA")
  if (identical(group_variable, "group_short")) {
    permanova <- vegan::adonis2(
      comp_matrix ~ group_short + sex + age,
      data = data, method = "euclidean", by = "margin",
      permutations = 5000
    )
  } else if (identical(group_variable, "group")) {
    permanova <- vegan::adonis2(
      comp_matrix ~ group + sex + age + source,
      data = data, method = "euclidean", by = "margin",
      permutations = 5000
    )
  } else {
    stop("Unsupported composition grouping variable: ", group_variable)
  }
  writeLines(
    capture.output(print(permanova)),
    file.path(table_dir, paste0(file_stub, "_PERMANOVA.txt"))
  )
  write.csv(as.data.frame(permanova),
            file.path(table_dir, paste0(file_stub, "_PERMANOVA.csv")))
  
  message("Running ", file_stub, ": component-wise models")
  component_rows <- lapply(seq_len(ncol(clr_data)), function(j) {
    model_data <- data
    model_data$.clr_component <- clr_data[, j]
    if (identical(group_variable, "group_short")) {
      fit <- lm(.clr_component ~ group_short + age + sex, data = model_data)
    } else {
      fit <- lm(.clr_component ~ group + age + sex + source,
                data = model_data)
    }
    coefficient_table <- summary(fit)$coefficients
    if (!group_coefficient %in% rownames(coefficient_table)) {
      stop("Expected coefficient '", group_coefficient,
           "' was not found for ", composition_names[j], ".")
    }
    data.frame(
      variable = composition_names[j],
      pval = coefficient_table[group_coefficient, "Pr(>|t|)"],
      tval = coefficient_table[group_coefficient, "t value"]
    )
  })
  component_results <- bind_rows(component_rows) %>%
    mutate(pval_fdr = p.adjust(pval, method = "BH"))
  write.csv(component_results,
            file.path(table_dir, paste0(file_stub, "_components.csv")),
            row.names = FALSE)
  
  message("Running ", file_stub, ": component plot")
  plot_data <- data.frame(
    comparison_group = data[[group_variable]],
    as.data.frame(clr_data, check.names = FALSE),
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = all_of(composition_names), names_to = "variable",
      values_to = "value"
    ) %>%
    left_join(component_results %>% select(variable, pval_fdr),
              by = "variable") %>%
    mutate(sig = case_when(
      pval_fdr < 0.001 ~ "***",
      pval_fdr < 0.01 ~ "**",
      pval_fdr < 0.05 ~ "*",
      TRUE ~ ""
    ))
  star_positions <- plot_data %>%
    group_by(variable) %>%
    summarise(ypos = max(value) + 0.1, sig = first(sig), .groups = "drop")
  
  composition_plot <- ggplot(
    plot_data,
    aes(x = comparison_group, y = value, color = comparison_group)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.2) +
    facet_wrap(~variable, scales = "free_y", ncol = 7) +
    geom_text(
      data = star_positions,
      aes(x = 1.5, y = ypos, label = sig),
      size = 7, inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(), legend.position = "none"
    ) +
    labs(
      title = plot_title,
      y = "Centre Log-Ratio Transformed Proportion of Total Lipids",
      x = plot_x
    )
  ggsave(file.path(figure_dir, paste0(file_stub, "_components.png")),
         composition_plot, width = 14, height = 8, dpi = 300)
  
  invisible(list(permanova = permanova, components = component_results))
}

composition_gen <- lipid_compositions %>%
  filter(group_short %in% c("GEN_0", "GEN_1")) %>%
  group_by(BlindedCode) %>%
  slice_max(Visit) %>%
  ungroup()
run_composition_analysis(
  data = composition_gen,
  group_variable = "group_short",
  permanova_covariates = c("sex", "age"),
  lm_covariates = c("age", "sex"),
  group_coefficient = "group_shortGEN_1",
  seed = 1,
  file_stub = "lipid_composition_GEN1_vs_GEN0",
  plot_title = "Presymptomatic Mutation Carriers v Non-Carriers",
  plot_x = "Group"
)

composition_all <- lipid_compositions %>%
  filter(group_short != "GEN_1") %>%
  mutate(
    group_short = ifelse(group_short == "GEN_2", "bvftd", group_short),
    sex = ifelse(sex == "f", 1, 0)
  ) %>%
  group_by(subject_id) %>%
  slice_max(Visit) %>%
  ungroup()
run_composition_analysis(
  data = composition_all,
  group_variable = "group",
  permanova_covariates = c("sex", "age", "source"),
  lm_covariates = c("age", "sex", "source"),
  group_coefficient = "grouppatient",
  seed = 123,
  file_stub = "lipid_composition_all_patients_vs_controls",
  plot_title = "All Patients v Controls",
  plot_x = "Group"
)


