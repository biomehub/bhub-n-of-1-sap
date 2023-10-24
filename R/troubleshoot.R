library(tidyverse)
library(here)
library(lme4)
library(furrr)
library(patchwork)

# TODO: simulate ite_predictor and pct-treat rneff from bivariate normal
# TODO: treat patient as factor

complement <- function(y, rho, x = NULL) {
  if (is.null(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

get_n_of_1_dataset <- function(
    n_patients,
    n_cycles,
    sd_patient,
    sd_random_slope,
    sd_residual,
    intercept,
    overall_treatment_effect,
    .seed,
    ite_explainable_fraction = 0.5
) {
  withr::with_seed(
    seed = .seed,
    code = {
      expand.grid(
        patient = 1:n_patients,
        treatment = LETTERS[1:n_treats],
        cycle = 1:n_cycles
      ) %>% 
        dplyr::arrange(patient, cycle, treatment) %>% 
        dplyr::group_by(patient) %>% 
        dplyr::mutate(
          patient_microbiome_score = rnorm(
            1, sd = sd_random_slope*sqrt(ite_explainable_fraction)  # sqrt(var * expl_frac)
          ),
          unexplainable_ite_variability = rnorm(
            1, sd = sd_random_slope*sqrt(1 - ite_explainable_fraction) # sqrt(var * (1 - expl_frac))
          ),
          patient_slope_rneff = ifelse(
            treatment == "A", 0,
            patient_microbiome_score + unexplainable_ite_variability
          ),
          patient_rneff = rnorm(1, sd = sd_patient), # TODO: correlation intercept-slope
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
          intercept = intercept,
          treatment_effect = ifelse(treatment == "A", 0, overall_treatment_effect),
          residual_error = rnorm(nrow(.), sd = sd_residual),
          true_individual_effect = treatment_effect + patient_slope_rneff,
          y = intercept + treatment_effect + patient_rneff + patient_slope_rneff + residual_error
        )
    }
  )
}

simulate_n_of_1_trial <- function(
    n_patients,
    n_cycles,
    sd_patient,
    sd_random_slope,
    sd_residual,
    intercept,
    overall_treatment_effect,
    .seed,
    ite_explainable_fraction = 0.5,
    plot = FALSE,
    raw_data = TRUE,
    keep_fit = FALSE,
    bootstrap_samples = 100,
    bootstrap_cpus = 1,
    n_treats = 2
) {
  stopifnot("Can only simulate 2 treatments" = n_treats == 2L)
  # simulate base data
  df <- withr::with_seed(
    seed = .seed,
    code = {
      expand.grid(
        patient = 1:n_patients,
        treatment = LETTERS[1:n_treats],
        cycle = 1:n_cycles
      ) %>% 
        dplyr::arrange(patient, cycle, treatment) %>% 
        dplyr::group_by(patient) %>% 
        dplyr::mutate(
          patient_microbiome_score = rnorm(
            1, sd = sd_random_slope*sqrt(ite_explainable_fraction)  # sqrt(var * expl_frac)
          ),
          unexplainable_ite_variability = rnorm(
            1, sd = sd_random_slope*sqrt(1 - ite_explainable_fraction) # sqrt(var * (1 - expl_frac))
          ),
          patient_slope_rneff = ifelse(
            treatment == "A", 0,
            patient_microbiome_score + unexplainable_ite_variability
          ),
          patient_rneff = rnorm(1, sd = sd_patient), # TODO: correlation intercept-slope
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
          intercept = intercept,
          treatment_effect = ifelse(treatment == "A", 0, overall_treatment_effect),
          residual_error = rnorm(nrow(.), sd = sd_residual),
          true_individual_effect = treatment_effect + patient_slope_rneff,
          y = intercept + treatment_effect + patient_rneff + patient_slope_rneff + residual_error
        )
    }
  )
  
  fit_full <- suppressMessages({
    lme4::lmer(
      y ~ treatment + treatment:patient_microbiome_score + (1 + treatment|patient),
      data = df,
      REML = TRUE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  fit_full_mle <- refitML(fit_full)
  
  fit_reduced_ite <- suppressMessages({
    lme4::lmer(
      y ~ treatment + treatment:patient_microbiome_score + (1|patient),
      data = df,
      REML = FALSE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  fit_reduced_microbiome <- suppressMessages({
    lme4::lmer(
      y ~ treatment + (1 + treatment|patient),
      data = df,
      REML = FALSE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  
  ite_pval <- anova(fit_full_mle, fit_reduced_ite)[["Pr(>Chisq)"]][2]
  hte_pval <- anova(fit_full_mle, fit_reduced_microbiome)[["Pr(>Chisq)"]][2]
  
  get_preds <- function(myfit) {
    .pred_data <- df %>% 
      select(treatment, patient, patient_microbiome_score) %>% 
      unique()
    .pred_data$y <- lme4:::predict.merMod(
      myfit,
      newdata = .pred_data,
      re.form = NULL
    )
    .pred_data %>% 
      tidyr::pivot_wider(
        names_from = treatment,
        values_from = y
      ) %>% 
      dplyr::arrange(patient) %>% 
      dplyr::mutate(d = B - A) %>% 
      dplyr::pull(d)
  }
  
  bootstrap_predictions <- lme4::bootMer(
    fit_full,
    nsim = bootstrap_samples, 
    FUN = get_preds, 
    ncpus = bootstrap_cpus
  )
  
  df_individual_effects <- tibble(
    patient = 1:n_patients,
    estimate = bootstrap_predictions$t0,
    se = apply(bootstrap_predictions$t, 2, sd),
    lower = estimate - qnorm(.975)*se,
    upper = estimate + qnorm(.975)*se
  )
  
  df_summaries <- df %>% 
    select(patient, treatment, cycle, y) %>% 
    pivot_wider(names_from = treatment, values_from = y) %>% 
    mutate(d = B - A) %>% 
    group_by(patient) %>% 
    summarise(
      avg_delta = mean(d),
      se_delta = sd(d)/sqrt(n_cycles),
      lower_delta = avg_delta - qnorm(0.975)*se_delta,
      upper_delta = avg_delta + qnorm(0.975)*se_delta,
    )
  
  df_individual_effects <- inner_join(
    df_individual_effects,
    df %>% 
      filter(treatment == "B") %>% 
      select(contains("patient"), contains("ite"), contains("true")) %>% 
      unique(),
    by = "patient"
  ) %>% 
    inner_join(
      df_summaries,
      by = "patient"
    )
  
  output <- list(
    ite_pval = ite_pval,
    hte_pval = hte_pval
  )
  
  if (isTRUE(keep_fit)) {
    output[["fit"]] <- fit_full
  }
  
  if (isTRUE(raw_data)) {
    output[["individual_effects"]] = df_individual_effects
    output[["df"]] <- df
    output[["summary"]] <- df_summaries
  }
  
  if (isTRUE(plot)) {
    theme_set(theme_bw(base_size = 14))
    plot_indiv_effects <- df_individual_effects %>% 
      ggplot(
        aes(
          fct_reorder(factor(patient), estimate, max), 
          estimate, 
          ymin = lower, ymax = upper
        )
      ) +
      geom_hline(yintercept = 0, lty = 2, alpha = 0.3, lwd = 1) +
      geom_pointrange() +
      geom_point(
        aes(y = true_individual_effect),
        pch = 18, size = 2, color = "red"
      ) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_blank(),
            plot.caption.position = "plot",
            plot.caption = element_text(hjust = 0)) +
      labs(
        x = "Individual patients",
        y = "Estimated individual diet effect",
        subtitle = "Average difference in posprandial blood glucose for each patient",
        caption = "Red diamonds are the true individual diet effects"
      )
    plot_indiv_obs_values <- df %>% 
      mutate(
        x = fct_reorder(
          factor(patient),
          y,
          mean
        )
      ) %>% 
      ggplot(aes(x, y)) +
      geom_boxplot(aes(color = paste0("Diet ", treatment)), alpha = 0) +
      scale_color_manual(
        values = c(
          "Diet A" = "#1B9E77", "Diet B" = "#D95F02"
        )
      ) +
      scale_y_continuous(breaks = scales::pretty_breaks(5)) +
      theme(
        axis.text.x = element_blank(),
        legend.position = "bottom"
      ) +
      labs(
        subtitle = "Observed posprandial blood glucose for each patient",
        x = "Individual patients",
        y = "Blood glucose (mg/dL)",
        color = NULL
      )
    plot_indiv_obs_deltas <- df %>% 
      select(patient, treatment, cycle, y) %>% 
      pivot_wider(names_from = treatment, values_from = y) %>% 
      mutate(d = B - A, patient = factor(patient)) %>% 
      ggplot(aes(fct_reorder(factor(patient), d, median), d)) +
      geom_hline(yintercept = 0, lty = 2, alpha = 0.3, lwd = 1) +
      geom_boxplot() +
      geom_jitter(width = .2) +
      theme(
        axis.text.x = element_blank()
      ) +
      labs(
        subtitle = "Observed differences in posprandial blood glucose for each patient",
        x = "Individual patients",
        y = "Delta blood glucose (mg/dL)\n(B-A)"
      )
    
    plot_outcome_distribution <- df %>% 
      mutate(treatment = paste0("Diet ", treatment)) %>% 
      ggplot(aes(y)) +
      geom_histogram(bins = 30, fill = "steelblue") +
      facet_wrap(~treatment) +
      labs(x = "Postprandial blood glucose (mg/dL)",
           y = "Frequency")
    
    plot_observed_differences_by_cycle <- df %>% 
      select(patient, treatment, cycle, y) %>% 
      pivot_wider(names_from = treatment, values_from = y) %>% 
      mutate(d = B - A, patient = factor(patient)) %>% 
      select(-A, -B) %>% 
      pivot_wider(names_from = cycle, values_from = d, names_prefix = "Cycle ") %>% 
      column_to_rownames("patient") %>% 
      GGally::ggpairs(
        lower = list(
          continuous = GGally::wrap(
            GGally::ggally_smooth, color = "gray40"
          ) 
        ),
        title = "Observed differences in posprandial blood glucose for each patient by trial cycle",
        xlab = "Delta blood glucose (B-A)", 
        ylab = "Delta blood glucose (B-A)",
        progress = FALSE
      )
    
    plot_observed_values_by_cycle <- df %>% 
      select(patient, treatment, cycle, y) %>% 
      pivot_wider(names_from = treatment, values_from = y) %>% 
      mutate(d = B - A,
             cycle = paste0("Cycle", cycle),
             x = fct_reorder(factor(patient), d, median)) %>% 
      ggplot(aes(x = x, xend = x, ymin = A, ymax = B)) +
      geom_errorbar(
        color = "gray40",
        position = position_dodge(width=.9), 
        show.legend = F, aes(group =cycle)
      ) +
      geom_point(aes(y = B, group = cycle, fill = "Diet B"),
                 pch = 21,
                 position = position_dodge(width=.9)) +
      geom_point(aes(y = A, group = cycle, fill = "Diet A"),
                 pch = 21,
                 position = position_dodge(width=.9)) +
      theme(
        axis.text.x = element_blank(),
        legend.position = 'top'
      ) +
      scale_fill_manual(
        values = list("Diet B" = "red", "Diet A" = "steelblue")
      ) +
      facet_wrap(~cycle) +
      labs(
        y = "Postprandial blood glucose (mg/dL)",
        x = "Individual patients",
        fill = NULL
      )
    
    plot_shrinkage <- df_individual_effects %>% 
      ggplot(
        aes(
          x = avg_delta,
          xmin = lower_delta,
          xmax = upper_delta,
          y = estimate,
          ymin = lower,
          ymax = upper
        )
      ) +
      geom_abline(lty = 2, alpha = 0.3) +
      geom_pointrange() +
      geom_errorbarh() +
      geom_point(
        aes(y = true_individual_effect),
        color = "red"
      ) +
      labs(
        x = "Summary-based estimates",
        y = "Shrunk estimates"
      )
    
    output[["plots"]] <- list(
      individual_effects = plot_indiv_effects,
      observed_differences = plot_indiv_obs_deltas,
      observed_values = plot_indiv_obs_values,
      observed_values_by_treatment = plot_outcome_distribution,
      observed_differences_by_cycle = plot_observed_differences_by_cycle,
      observed_values_by_cycle = plot_observed_values_by_cycle,
      shrinkage = plot_shrinkage
    )
  }
  return(output)
}

# TODO: check why more cycles
# is giving higher SE's
# check if higher n patients gives lower SE's too
# TODO: check accuracy of individual tret eff estimates

trial <- simulate_n_of_1_trial(
  n_patients = 50,
  n_cycles = 5,
  sd_patient = 8,
  sd_random_slope = 4,
  sd_residual = 5.85,
  intercept = 125,
  overall_treatment_effect = -5,
  ite_explainable_fraction = 0.25,
  .seed = 9323,
  plot = TRUE,
  raw_data = TRUE,
  keep_fit = TRUE,
  bootstrap_samples = 50,
  bootstrap_cpus = 1
)

cool_plot <- (
  trial$plots$individual_effects+coord_cartesian(ylim = c(-15, 15))
)/(
  trial$plots$observed_differences+coord_cartesian(ylim = c(-15, 15))
)

trial$plots$individual_effects$data$true_individual_effect %>% 
  range()


# par setting 1

settings <- expand.grid(
  n_pct = seq(45, 100, 5),
  n_cycles = 4:7
)

settings$power_ite <- NA
settings$power_hte <- NA

for (i in seq_len(nrow(settings))) {
  .n_pct <- settings[i, "n_pct"]
  .n_cyc <- settings[i, "n_cycles"]
  x <- paste0(
    "Setting ", i,
    " n patients = ", .n_pct,
    " n cycles = ", .n_cyc,
    "\n"
  )
  cat(x)
  
  n_sim <- 500
  plan(multisession, workers = 7)
  simulations <- future_map(
    1:n_sim, ~ {
      simulate_n_of_1_trial(
        n_patients = .n_pct,
        n_cycles = .n_cyc,
        sd_patient = sqrt(30),
        sd_random_slope = sqrt(6.5),
        sd_residual = sqrt(30),
        intercept = 130,
        overall_treatment_effect = -5,
        ite_explainable_fraction = 0.25,
        .seed = .x,
        plot = FALSE,
        raw_data = FALSE,
        keep_fit = FALSE,
        bootstrap_samples = 1,
        bootstrap_cpus = 1
      )
    },
    .options = furrr::furrr_options(seed = 12345),
    .progress = TRUE
  )
  plan(sequential)
  closeAllConnections()
  
  settings[i, "power_ite"] <- mean(
    map_dbl(simulations, ~ .x$ite_pval) < 0.05
  )
  settings[i, "power_hte"] <- mean(
    map_dbl(simulations, ~ .x$hte_pval) < 0.05
  )
}



settings %>% 
  pivot_longer(contains('power')) %>% 
  mutate(
    name = factor(
      str_extract(name, "ite|hte") %>% 
        str_to_upper()
    ),
    n_cycles = ordered(n_cycles),
    n_pct = ordered(n_pct)
  ) %>% 
  ggplot(aes(n_cycles, n_pct, fill = value)) +
  geom_tile() +
  geom_text(aes(label = paste0(round(value*100), "%"))) +
  facet_wrap(~name) +
  scale_fill_viridis_b()


# par setting 2

"
delta individual effects of ~10
between 10% least-responder
and 10% top-responder
"

qnorm(c(.1, .9), sd=4)
qnorm(c(.1, .9), sd=5.85)
qnorm(c(.1, .9), sd=8)

settings2 <- expand.grid(
  n_pct = seq(45, 100, 5),
  n_cycles = 4:7
)

settings2$power_ite <- NA
settings2$power_hte <- NA

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
  
  n_sim <- 2000
  plan(multisession, workers = 7)
  simulations <- future_map(
    1:n_sim, ~ {
      simulate_n_of_1_trial(
        n_patients = .n_pct,
        n_cycles = .n_cyc,
        sd_patient = 8,
        sd_random_slope = 4,
        sd_residual = 5.85,
        intercept = 125,
        overall_treatment_effect = -5,
        ite_explainable_fraction = 0.25,
        .seed = .x + 10,
        plot = FALSE,
        raw_data = FALSE,
        keep_fit = FALSE,
        bootstrap_samples = 1,
        bootstrap_cpus = 1
      )
    },
    .options = furrr::furrr_options(seed = 2),
    .progress = TRUE
  )
  plan(sequential)
  closeAllConnections()
  
  settings2[i, "power_ite"] <- mean(
    map_dbl(simulations, ~ .x$ite_pval) < 0.05
  )
  settings2[i, "power_hte"] <- mean(
    map_dbl(simulations, ~ .x$hte_pval) < 0.05
  )
}



settings2 %>% 
  pivot_longer(contains('power')) %>% 
  mutate(
    name = factor(
      str_extract(name, "ite|hte") %>% 
        str_to_upper()
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
  "output/power.png",
  width = 9, height = 4, dpi = 600
)

trial <- simulate_n_of_1_trial(
  n_patients = 100,
  n_cycles = 5,
  sd_patient = 8,
  sd_random_slope = 4,
  sd_residual = 5.85,
  intercept = 125,
  overall_treatment_effect = -5,
  ite_explainable_fraction = 0.25,
  .seed = 342342,
  plot = T,
  raw_data = T,
  keep_fit = T,
  bootstrap_samples = 300,
  bootstrap_cpus = 1
)

hist(rnorm(3e6, -5, sd = 4), probability = T, breaks = 50,
     main = "Distrbuição do efeito individual sobre GPP\npara uma dieta com efeito médio de -5 mg/dL",
     xlab = "Delta Glicemia pós-prandial (mg/dL)")

qnorm(c(.1, .9), sd=4)
qnorm(c(.05, .95), sd=4)
qnorm(c(.05/2, 1 - .05/2),-5, sd=4)