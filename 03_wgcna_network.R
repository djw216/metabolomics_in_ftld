#!/usr/bin/env Rscript

# Metabolomics study: WGCNA network construction and module annotation
#
# Prerequisite: run 01_preprocess_impute.R. 
#
# Outputs:
#   data/processed/wgcna_network.rds
#   results/wgcna_network/tables/
#   results/wgcna_network/figures/

required_packages <- c("dplyr", "ggplot2", "WGCNA", "scales")
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
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# ---- Paths and network parameters ------------------------------------------
project_dir <- normalizePath(
  Sys.getenv("METABOLOMICS_PROJECT_DIR", unset = "."), mustWork = FALSE
)
processed_dir <- Sys.getenv(
  "METABOLOMICS_OUTPUT_DIR",
  unset = file.path(project_dir, "data", "processed")
)
results_dir <- file.path(project_dir, "results", "wgcna_network")
table_dir <- file.path(results_dir, "tables")
figure_dir <- file.path(results_dir, "figures")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

candidate_powers <- 1:20
soft_power <- as.integer(Sys.getenv("WGCNA_SOFT_POWER", unset = "3"))
minimum_module_size <- as.integer(
  Sys.getenv("WGCNA_MIN_MODULE_SIZE", unset = "10")
)
deep_split <- as.integer(Sys.getenv("WGCNA_DEEP_SPLIT", unset = "1"))
hub_threshold <- 0.75

# ---- Load corrected metabolite data ----------------------------------------
preprocessed <- readRDS(
  file.path(processed_dir, "metabolomics_preprocessed.rds")
)
metabolite_names <- preprocessed$metabolite_names

missing_metabolites <- setdiff(
  metabolite_names, names(preprocessed$corrected_metabolites)
)
if (length(missing_metabolites) > 0L) {
  stop("Corrected data are missing metabolites: ",
       paste(missing_metabolites, collapse = ", "))
}

sample_metadata <- preprocessed$corrected_metabolites %>%
  select(PARENT_SAMPLE_NAME) %>%
  left_join(preprocessed$subjects,
            by = c("PARENT_SAMPLE_NAME" = "sample_id"))

# WGCNA expects observations in rows and metabolites in columns.
dat_expr <- as.data.frame(
  preprocessed$corrected_metabolites[, metabolite_names, drop = FALSE],
  check.names = FALSE
)
rownames(dat_expr) <- preprocessed$corrected_metabolites$PARENT_SAMPLE_NAME

quality <- WGCNA::goodSamplesGenes(dat_expr, verbose = 3)
if (!quality$allOK) {
  excluded_samples <- rownames(dat_expr)[!quality$goodSamples]
  excluded_metabolites <- names(dat_expr)[!quality$goodGenes]
  writeLines(excluded_samples,
             file.path(table_dir, "excluded_samples.txt"))
  writeLines(excluded_metabolites,
             file.path(table_dir, "excluded_metabolites.txt"))
  dat_expr <- dat_expr[quality$goodSamples, quality$goodGenes, drop = FALSE]
  sample_metadata <- sample_metadata[quality$goodSamples, , drop = FALSE]
  metabolite_names <- names(dat_expr)
}

# ---- Select and document the soft-thresholding power -----------------------
soft_threshold <- WGCNA::pickSoftThreshold(
  dat_expr,
  powerVector = candidate_powers,
  networkType = "unsigned",
  verbose = 5
)
fit_indices <- as.data.frame(soft_threshold$fitIndices)
write.csv(fit_indices,
          file.path(table_dir, "soft_threshold_diagnostics.csv"),
          row.names = FALSE)

png(file.path(figure_dir, "soft_threshold_diagnostics.png"),
    width = 1800, height = 850, res = 160)
par(mfrow = c(1, 2))
with(fit_indices, {
  plot(Power, -sign(slope) * SFT.R.sq,
       xlab = "Soft-thresholding power",
       ylab = "Signed scale-free topology fit (R²)", type = "n")
  text(Power, -sign(slope) * SFT.R.sq, labels = Power, col = "red")
  abline(h = 0.80, col = "red", lty = 2)
  plot(Power, mean.k., xlab = "Soft-thresholding power",
       ylab = "Mean connectivity", type = "n")
  text(Power, mean.k., labels = Power, col = "red")
})
dev.off()


adjacency_matrix <- WGCNA::adjacency(
  dat_expr, power = soft_power, type = "unsigned"
)
tom <- WGCNA::TOMsimilarity(adjacency_matrix, TOMType = "unsigned")
diss_tom <- 1 - tom
metabolite_tree <- hclust(as.dist(diss_tom), method = "average")

dynamic_modules <- cutreeDynamic(
  dendro = metabolite_tree,
  distM = diss_tom,
  deepSplit = deep_split,
  pamRespectsDendro = FALSE,
  minClusterSize = minimum_module_size,
  respectSmallClusters = TRUE,
  method = "hybrid"
)
module_colors <- WGCNA::labels2colors(dynamic_modules)

png(file.path(figure_dir, "metabolite_dendrogram_modules.png"),
    width = 1800, height = 1000, res = 160)
WGCNA::plotDendroAndColors(
  metabolite_tree, module_colors, "Dynamic tree cut",
  dendroLabels = FALSE, hang = 0.03,
  main = "Metabolite clustering and WGCNA modules"
)
dev.off()

# ---- Eigengenes, membership and annotation ---------------------------------
eigengenes <- WGCNA::moduleEigengenes(
  dat_expr, colors = module_colors
)$eigengenes
eigengenes <- WGCNA::orderMEs(eigengenes)

module_membership <- stats::cor(
  dat_expr, eigengenes, use = "pairwise.complete.obs"
)
membership_long <- as.data.frame(as.table(module_membership))
names(membership_long) <- c("CHEMICAL_ID", "module_eigengene", "kME")

module_table <- tibble(
  CHEMICAL_ID = names(dat_expr),
  module_color = module_colors
) %>%
  left_join(preprocessed$lookup, by = "CHEMICAL_ID") %>%
  left_join(membership_long, by = "CHEMICAL_ID") %>%
  mutate(
    own_module_eigengene = paste0("ME", module_color),
    is_own_module = module_eigengene == own_module_eigengene
  )

own_membership <- module_table %>%
  filter(is_own_module) %>%
  select(-is_own_module, -own_module_eigengene) %>%
  arrange(module_color, desc(abs(kME)))

module_summary <- own_membership %>%
  group_by(module_color) %>%
  summarise(
    n_metabolites = n(),
    n_hubs_abs_kME_gt_0_75 = sum(abs(kME) > hub_threshold, na.rm = TRUE),
    top_sub_pathway = {
      tab <- sort(table(SUB_PATHWAY), decreasing = TRUE)
      if (length(tab)) names(tab)[1] else NA_character_
    },
    .groups = "drop"
  )

hub_metabolites <- own_membership %>%
  filter(abs(kME) > hub_threshold) %>%
  arrange(module_color, desc(abs(kME)))

write.csv(own_membership, file.path(table_dir, "module_membership.csv"),
          row.names = FALSE)
write.csv(module_summary, file.path(table_dir, "module_summary.csv"),
          row.names = FALSE)
write.csv(hub_metabolites, file.path(table_dir, "hub_metabolites.csv"),
          row.names = FALSE)

# Module composition is displayed as proportions.
composition <- own_membership %>%
  filter(!is.na(SUPER_PATHWAY), SUPER_PATHWAY != "") %>%
  count(module_color, SUPER_PATHWAY, name = "n") %>%
  group_by(module_color) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
write.csv(composition, file.path(table_dir, "module_composition.csv"),
          row.names = FALSE)


eigengene_table <- bind_cols(
  data.frame(PARENT_SAMPLE_NAME = rownames(dat_expr)),
  as.data.frame(eigengenes, check.names = FALSE)
)

output <- list(
  eigengenes = eigengene_table,
  module_colors = setNames(module_colors, names(dat_expr)),
  module_membership = own_membership,
  module_summary = module_summary,
  sample_metadata = sample_metadata,
  metabolite_names = names(dat_expr),
  network_parameters = list(
    network_type = "unsigned",
    soft_power = soft_power,
    minimum_module_size = minimum_module_size,
    deep_split = deep_split,
    hub_threshold = hub_threshold,
    candidate_powers = candidate_powers
  ),
  soft_threshold_fit = fit_indices,
  metabolite_tree = metabolite_tree,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

saveRDS(output, file.path(processed_dir, "wgcna_network.rds"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))

