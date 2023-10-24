library(tidyverse)
library(here)
library(lme4)
library(furrr)
library(patchwork)

source(here("R/functions.R"))

THREADS <- 10
n_sim <- 600
overall_seed <- 29102023
set.seed(overall_seed)
sim_seeds <- sample(n_sim * 10, size = n_sim)

settings2 <- expand.grid(
  n_pct = seq(50, 80, 5),
  n_cycles = 5
)
settings2 <- data.frame(n_pct = 80, n_cycles = 5)
results <- vector("list", nrow(settings2))
result805 <- NULL
for (i in seq_len(nrow(settings2))) {
  .n_pct <- settings2[i, "n_pct"]
  .n_cyc <- settings2[i, "n_cycles"]
  x <- paste0(
    "\nSetting ", i,
    " n patients = ", .n_pct,
    " n cycles = ", .n_cyc,
    "\n"
  )
  cat(x)
  
  plan(multisession, workers = THREADS)
  simulations <- future_map_dfr(
    1:n_sim, ~ {
      .trial <- simulate_n_of_1_trial(
        n_patients = .n_pct,
        n_cycles = .n_cyc,
        sd_patient = 12.75,
        sd_random_slope = 4,
        sd_residual = 6,
        intercept = 60,
        overall_treatment_effect = 0,
        ite_explainable_fraction = 0.2,
        outcome_microbiome_association = 0,
        .seed = sim_seeds[.x],
        plot = FALSE,
        raw_data = FALSE,
        keep_fit = TRUE,
        bootstrap_samples = 1,
        bootstrap_cpus = 1
      )
      s <- summary(.trial$fit)$coefficients
      vcov <- lme4::VarCorr(.trial$fit)$patient
      vcov2 <- lme4::VarCorr(.trial$fit2)$patient
      data.frame(
        setting_id = i,
        n_pct = .n_pct,
        n_cycles = .n_cyc,
        pval_ite = .trial$ite_pval,
        pval_ite2 = .trial$ite_pval2,
        pval_hte = .trial$hte_pval,
        sd_slope = attr(vcov, "stddev")[[2]],
        sd_slope2 = attr(vcov2, "stddev")[[2]],
        overall_effect = s[2, 1],
        interaction = s[4, 1]
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
        "output/tmp/power-sim-setting-{i}.tsv"
      )
    )
  )
  if (.n_pct == 80 && .n_cyc == 5){
    result805 <- simulations
  }
  results[[i]] <- simulations %>% 
    summarise(
      setting_id = unique(setting_id),
      n_pct = unique(n_pct),
      n_cycles = unique(n_cycles),
      power_ite = mean(pval_ite < 0.05),
      power_ite2 = mean(pval_ite2 < 0.05),
      power_hte = mean(pval_hte < 0.05),
      avg_sd_slope = mean(sd_slope),
      avg_overall_effect = mean(overall_effect),
      avg_interaction = mean(interaction)
    )
}

# sigle trial
result805

# results across trials

results <- bind_rows(results)

write_tsv(
  results,
  here("output/simulation-results2.tsv")
)
.p <- results %>%
  select(setting_id, n_pct, n_cycles, contains("power")) %>% 
  pivot_longer(contains('power')) %>% 
  mutate(
    name = factor(
      str_extract(name, "ite\\d*|hte") %>% 
        str_to_upper() %>% 
        str_replace("2", "-2")
    ),
    n_cycles = ordered(n_cycles),
    n_pct = ordered(n_pct)
  ) %>% 
  ggplot(aes(n_cycles, n_pct, fill = value*100)) +
  geom_tile() +
  geom_text(aes(label = paste0(round(value*100), "%")),
            color = 'black') +
  facet_wrap(~name, scales = "free") +
  scale_fill_viridis_b() +
  theme_bw(base_size = 14) +
  labs(
    x = "Cycles",
    y = "Patients",
    fill = "Power (%)"
  )

ggsave(
  here("output/power-plot2.png"),
  .p, width = 9, height = 4, dpi = 600
)
