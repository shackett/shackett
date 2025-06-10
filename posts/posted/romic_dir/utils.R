prep_shiny_gasch2k <- function () {
  
  gasch_2000 <- readr::read_tsv(
    file = "http://www.shackett.org/files/gasch2000.txt",
    col_types = readr::cols() # to accept default column types
  )
  gasch_matrix <- gasch_2000 %>%
    select(-UID, -NAME, -GWEIGHT) %>%
    as.matrix()
  rownames(gasch_matrix) <- gasch_2000$UID 
 
  goslim_mappings <- readr::read_tsv(
    "https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab",
    col_names = c("ORF", "common", "SGD", "category", "geneset", "GO", "class")
  ) %>%
    select(-GO) %>%
    group_by(ORF, category) %>%
    slice(1) %>%
    tidyr::spread(category, geneset) %>%
    select(
      ORF, common, SGD, class,
      cellular_compartment = C,
      molecular_function = F,
      biological_process = P
    ) %>%
    ungroup()
  
  feature_metadata <- gasch_2000 %>%
    select(UID) %>%
    left_join(goslim_mappings, by = c("UID" = "ORF"))
  
  experiment_labels <- tibble::tibble(sample = colnames(gasch_matrix)) %>%
    mutate(experiment = case_when(
      str_detect(sample, "hs\\-1") ~ "Heat Shock (A) (duration)",
      str_detect(sample, "hs\\-2") ~ "Heat Shock (B) (duration)",
      str_detect(sample, "^37C to 25C") ~ "Cold Shock (duration)",
      str_detect(sample, "^heat shock") ~ "Heat Shock (severity)",
      str_detect(sample, "^29C to 33C") ~ "29C to 33C (duration)",
      str_detect(sample, "^29C \\+1M sorbitol to 33C \\+ 1M sorbitol") ~ "29C + Sorbitol to 33C + Sorbitol (duration)",
      str_detect(sample, "^29C \\+1M sorbitol to 33C \\+ \\*NO sorbitol") ~ "29C + Sorbitol to 33C (duration)",
      str_detect(sample, "^constant 0.32 mM H2O2") ~ "Hydrogen peroxide (duration)",
      str_detect(sample, "^1 ?mM Menadione") ~ "Menadione (duration)",
      str_detect(sample, "^2.5mM DTT") ~ "DTT (A) (duration)",
      str_detect(sample, "^dtt") ~ "DTT (B) (duration)",
      str_detect(sample, "diamide") ~ "Diamide (duration)",
      str_detect(sample, "^1M sorbitol") ~ "Sorbitol (duration)",
      str_detect(sample, "^Hypo-osmotic shock") ~ "Hypo-Osmotic Shock (duration)",
      str_detect(sample, "^aa starv") ~ "Amino Acid Starvation (duration)",
      str_detect(sample, "^Nitrogen Depletion") ~ "Nitrogen Depletion (duration)",
      str_detect(sample, "^[Dd]iauxic [Ss]hift") ~ "Diauxic Shift (duration)",
      str_detect(sample, "ypd-2") ~ "YPD (duration)",
      str_detect(sample, "ypd-1") ~ "YPD stationary phase (duration)",
      str_detect(sample, "overexpression") ~ "TF Overexpression",
      str_detect(sample, "car-1") ~ "Carbon Sources (A)",
      str_detect(sample, "car-2") ~ "Carbon Sources (B)",
      str_detect(sample, "ct-1") ~ "Temperature Gradient",
      str_detect(sample, "ct-2") ~ "Temperature Gradient, Steady State"
    )) %>%
    group_by(experiment) %>%
    mutate(experiment_order = 1:n()) %>%
    ungroup() %>%
    mutate(
      experiment_order = ifelse(is.na(experiment), NA, experiment_order),
      experiment = ifelse(is.na(experiment), "Other", experiment)
    )
  
  tall_gasch <- gasch_2000 %>%
    select(-NAME, -GWEIGHT) %>%
    tidyr::gather("sample", "expression", -UID) %>%
    dplyr::filter(!is.na(expression))
  
  triple_omic <- create_triple_omic(
    measurement_df = tall_gasch,
    feature_df = feature_metadata,
    sample_df = experiment_labels,
    feature_pk = "UID",
    sample_pk = "sample"
  )
  
  imputed_triple <- triple_omic %>%
    # overwrite existing expression so that we don't have the
    # raw expression changes which contains lots of missing values
    impute_missing_values(impute_var_name = "expression")
  
  triple_omic_gene_subset <- imputed_triple %>%
    filter_tomic(
      filter_type = "category",
      filter_table = "features",
      filter_variable = "cellular_compartment",
      filter_value = c(
        "mitochondrion"
      ))
  
  heatshock_triple_omic <- imputed_triple %>%
    filter_tomic(
      filter_type = "category",
      filter_table = "samples",
      filter_variable = "experiment",
      filter_value = c(
        "Heat Shock (A) (duration)",
        "Heat Shock (B) (duration)",
        "Heat Shock (severity)",
        "Temperature Gradient"
      ))
  
  saveRDS(triple_omic_gene_subset, file = "/home/sean/projects/shackett_blog/posts/staging/romic_dir/pcs/gasch2K.Rds")
  saveRDS(heatshock_triple_omic, file = "/home/sean/projects/shackett_blog/posts/staging/romic_dir/heatmap/heatshock_gasch2K.Rds")
  
}