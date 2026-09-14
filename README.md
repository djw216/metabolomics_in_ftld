# Metabolomics and Lipidomics in FTLD
Repository containing R files processing and visualisation for "Predictive early and late metabolomic changes in frontotemporal lobar degeneration syndromes"

All processing and analysis performed in R version 4.4.1

Metabolomic analysis:
Random Forest imputation, log transformation, missingness, RUV-III batch correction is performed in 01_preprocess_impute

Metabolome-wide group comparisons and single sample metabolite set enrichment analysis are performed in 02_group_comparisons_visualisations, and weighted correlation network analysis in 03_wgcna_network and 04_wgcna_associations_visualisations.

Survival analysis is in 05_prepare_survival_data and 06_survival_models_visualisations.

Lipidomic analysis:
Random Forest imputation, log transformation, missingness is set out in 07_prepare_lipidomics_data, with group comparisons for lipid species, lipid classes, and lipid compositions for all patients and in presymptomatic mutation carriers in 08_lipidomics_analyses_visualisations
