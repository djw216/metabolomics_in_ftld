# Metabolomics and Lipidomics in FTLD
Repository containing R files processing and visualisation for "Predictive early and late metabolomic changes in frontotemporal lobar degeneration syndromes"

All processing and analysis performed in R version 4.4.1

**Metabolomic analysis**
01_preprocess_impute.R: data preparation, missingness filtering, log transformation, random forest imputation and RUV-III batch correction.

02_group_comparisons_visualisations.R: metabolome-wide group comparisons, single-sample metabolite set enrichment analysis and associated visualisations.

03_wgcna_network.R: construction of weighted metabolite correlation networks.

04_wgcna_associations_visualisations.R: analysis and visualisation of metabolite-module associations.

05_prepare_survival_data.R: preparation of clinical outcomes and metabolomic data for survival analyses.

06_survival_models_visualisations.R: survival modelling, prediction analyses and associated visualisations.


**Lipidomic analysis**
07_prepare_lipidomics_data.R: data preparation, missingness filtering, log transformation and random forest imputation.

08_lipidomics_analyses_visualisations.R: comparisons of lipid species, lipid-class principal components and lipid compositions between all patients and controls, and between presymptomatic mutation carriers and non-carriers.
