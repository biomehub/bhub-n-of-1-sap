
run_power_simulation <- function(
  n_sim = 10,
  n_threads = 10,
  overall_seed = 1,
  sim_settings = data.frame(n_patients = 80),
  overwrite = FALSE
) {
  set.seed(overall_seed)
  sim_seeds <- sample(n_sim * 10, size = n_sim)
  results <- vector("list", nrow(sim_settings))
  for (i in seq_len(nrow(sim_settings))) {
    .n_patients <- sim_settings[i, "n_patients"]
    
    plan(multisession, workers = n_threads)
    simulations <- future_map_dfr(
      1:n_sim, ~ {
        run_seed <- sim_seeds[.x]
        filename <- paste0("simulation-run-", i, "-", .x, "-seed-", run_seed, ".rds")
        output_dir <- file.path("output", "simulation_runs", paste0("setting-", i), paste0("run-", .x))
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        output_rds_file <- file.path(output_dir, filename)
        if (file.exists(output_rds_file) & isFALSE(overwrite)) {
          .trial <- readRDS(output_rds_file)
        } else {
          .trial <- sample_trial_dataset(n_patients = .n_patients, .seed = run_seed) %>% 
            run_primary_analysis()
            saveRDS(.trial, output_rds_file)
        }

        data.frame(
          setting_id = i,
          run_id = .x,
          run_seed = run_seed,
          n_completed = .trial$n_completed,
          n_dropped = .trial$n_dropped,
          patient_pval = .trial$analysis$patient_pval,
          microbiome_pval = .trial$analysis$microbiome_pval,
          patient_padj = .trial$analysis$patient_padj,
          microbiome_padj = .trial$analysis$microbiome_padj,
          microbiome_interaction_effect = lme4::fixef(.trial$analysis$fit)[["m:dietB"]],
          sd_random_slope = attr(lme4::VarCorr(.trial$analysis$fit)$patient, "stddev")[["dietB"]]
        )
      },
      .options = furrr::furrr_options(seed = overall_seed),
      .progress = TRUE
    )
    plan(sequential)
    closeAllConnections()
    
    write_tsv(
      simulations,
      here(
        str_glue(
          "output/simulation_runs/power-sim-setting-{i}-seed-{overall_seed}.tsv"
        )
      )
    )

    results[[i]] <- simulations %>% 
      summarise(
        setting_id = unique(setting_id),
        n_simulation_runs = n(),
        n_completed_avg = mean(n_completed),
        n_dropped_avg = mean(n_dropped),
        dropout_rate_avg = mean(n_dropped/(n_completed + n_dropped)),
        power_patient = mean(patient_pval < 0.05),
        power_microbiome = mean(microbiome_pval < 0.05),
        power_patient_padj = mean(patient_padj < 0.05),
        power_microbiome_padj = mean(microbiome_padj < 0.05),
        microbiome_interaction_effect = mean(microbiome_interaction_effect),
        avg_sd_random_slope = mean(sd_random_slope)
      )
  }

  write_tsv(
      bind_rows(results),
      here(
        str_glue(
          "output/simulation-results.tsv"
        )
      )
    )

}