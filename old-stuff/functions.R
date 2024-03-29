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
    outcome_microbiome_association = 0.5,
    prediabetes_prevalence = 0,
    prediabetes_intercept_multiplier = 1,
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
          patient_rneff = rnorm(1, sd = sd_patient)
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
          intercept = ifelse(
            patient > n_patients*prediabetes_prevalence,
            intercept,
            intercept * prediabetes_intercept_multiplier 
          ),
          treatment_effect = ifelse(treatment == "A", 0, overall_treatment_effect),
          residual_error = rnorm(nrow(.), sd = sd_residual),
          true_individual_effect = treatment_effect + patient_slope_rneff,
          microbiome_effect = outcome_microbiome_association * patient_microbiome_score,
          y = intercept + treatment_effect + microbiome_effect + patient_rneff + patient_slope_rneff + residual_error
        ) %>% 
        dplyr::select(-microbiome_effect)
    }
  )
  
  fit_full <- suppressMessages({
    lme4::lmer(
      y ~ treatment*patient_microbiome_score + (1 + treatment|patient),
      data = df,
      REML = TRUE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  fit_full_mle <- refitML(fit_full)
  
  fit_reduced_ite <- suppressMessages({
    lme4::lmer(
      y ~ treatment*patient_microbiome_score + (1|patient),
      data = df,
      REML = FALSE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  fit_reduced_microbiome <- suppressMessages({
    lme4::lmer(
      y ~ treatment + patient_microbiome_score + (1 + treatment|patient),
      data = df,
      REML = FALSE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  
  ite_pval <- anova(fit_full_mle, fit_reduced_ite)[["Pr(>Chisq)"]][2]
  hte_pval <- anova(fit_full_mle, fit_reduced_microbiome)[["Pr(>Chisq)"]][2]
  
  # testing without microbiome score
  fit_full2 <- suppressMessages({
    lme4::lmer(
      y ~ treatment + (1 + treatment|patient),
      data = df,
      REML = TRUE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  fit_reduced2 <- suppressMessages({
    lme4::lmer(
      y ~ treatment + (1|patient),
      data = df,
      REML = TRUE,
      control = lmerControl(optCtrl = list(maxfn = 200))
    )
  })
  ite_pval2 <- anova(fit_full2, fit_reduced2, refit = FALSE)[["Pr(>Chisq)"]][2]
  
  
  get_preds <- function(myfit) {
    .pred_data <- expand.grid(
      patient = unique(df$patient),
      treatment = unique(df$treatment),
      patient_microbiome_score = mean(df$patient_microbiome_score)
    )
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
    ite_pval2 = ite_pval2,
    hte_pval = hte_pval
  )
  
  if (isTRUE(keep_fit)) {
    output[["fit"]] <- fit_full
    output[["fit2"]] <- fit_full2
  }
  
  if (isTRUE(raw_data)) {
    output[["individual_effects"]] = df_individual_effects
    output[["df"]] <- df
    output[["summary"]] <- df_summaries
  }
  
  if (isTRUE(plot)) {
    theme_set(theme_bw(base_size = 14))
    plot_indiv_effects <- df_individual_effects %>% 
      mutate(
        # because we predicted at avg(microbiome score)
        true_individual_effect = true_individual_effect - patient_microbiome_score + mean(patient_microbiome_score)
      ) %>% 
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
      geom_abline(lty = 2, alpha = 0.5) +
      geom_pointrange() +
      geom_errorbarh() +
      geom_point(
        aes(y = true_individual_effect),
        color = "red"
      ) +
      labs(
        x = "Summary-based estimates",
        y = "Shrunk estimates"
      ) +
      coord_cartesian(ylim = )
    
    plot_cond_effects <- marginaleffects::plot_cme(
      fit_full, 
      "treatment", 
      "patient_microbiome_score",
      re.form = ~ 0
    )
    
    output[["plots"]] <- list(
      individual_effects = plot_indiv_effects,
      conditional_effects = plot_cond_effects,
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
