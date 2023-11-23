sample_trial_dataset <- function(
    n_patients = 80,
    n_periods = 10,
    a0 = 200,
    b0 = 0,
    t1 = 5,
    t2 = 5,
    sd_epsilon = 15,           # within-patient observations vary   +/- 30  (avg. patient: 170 -- 230)
    sd1 = sqrt(40^2 - t1^2),   # individual averages vary           +/- 80  (patients' avgs: 120 -- 280)
    sd2 = sqrt(10^2 - t2^2),   # individual effects  vary           +/- 20  (MCID = 40)
    g_j = rep(0, n_periods),
    possible_sequences = c("ABBABAABBA", "BAABABBAAB"),
    dropout_rate = 0.2,
    .seed = 123
) {
  withr::with_seed(
    seed = .seed,
    code = {

      # generate metadata with random treatment sequence
      df <- expand.grid(
        patient = factor(
          paste0("pt-", 1:n_patients),
          levels = paste0("pt-", 1:n_patients)
        ),
        period = 1:n_periods
      ) %>%
        dplyr::arrange(patient, period) %>%
        dplyr::group_by(patient) %>%
        dplyr::mutate(
          sequence = factor(sample(possible_sequences, 1)),
          cycle = dplyr::case_when(
            period <= 2 ~ 1,
            period <= 4 ~ 2,
            period <= 6 ~ 3,
            period <= 8 ~ 4,
            period <= 10 ~ 5,
            TRUE ~ NA
          ),
          cycle = factor(cycle)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          diet = purrr::map2_chr(period, sequence, \(j, s) {
            stringr::str_split(s, "")[[1]][j]
          }),
          period = factor(period),
          m = NA,
          u1 = NA,
          u2 = NA,
          a = NA,
          b = NA,
          x = NA,
          epsilon = NA,
          iauc = NA
        )

      # sample random variables
      M <- as.vector(scale(rnorm(n = n_patients)))
      u1 <- rnorm(n = n_patients, sd = sd1)
      u2 <- rnorm(n = n_patients, sd = sd2)
      epsilon <- matrix(
        rnorm(n = n_patients * n_periods, sd = sd_epsilon),
        nrow = n_patients, ncol = n_periods
      )
      a <- a0 + t1 * M + u1
      b <- b0 + t2 * M + u2
      ## compute Y given everything else
      X <- Y <- matrix(nrow = n_patients, ncol = n_periods)
      for (i in seq_len(n_patients)) {
        patient_i_ix <- df$patient == paste0("pt-", i)
        df[patient_i_ix, "m"] <- M[i]
        df[patient_i_ix, "u1"] <- u1[i]
        df[patient_i_ix, "u2"] <- u2[i]
        df[patient_i_ix, "a"] <- a[i]
        df[patient_i_ix, "b"] <- b[i]
        for (j in seq_len(n_periods)) {
          ij_ix <- df$patient == paste0("pt-", i) & df$period == j
          stopifnot(sum(ij_ix) == 1)
          x_ij <- df$diet[ij_ix] == "B"
          df[ij_ix, "x"] <- x_ij
          df[ij_ix, "epsilon"] <- epsilon[i, j]
          df[ij_ix, "iauc"] <- a[i] + b[i] * x_ij + g_j[j] + epsilon[i, j]
        }
      }

      df_dropped <- NULL
      n_dropped <- 0
      n_completed <- n_patients
      if (dropout_rate > 0) {
        n_dropped <- round(n_patients * dropout_rate)
        n_completed <- n_patients - n_dropped
        dropped_patients <- base::sample(seq_len(n_patients), n_dropped, replace = FALSE)
        drop_ix <- df$patient %in% paste0("pt-", dropped_patients)
        df_dropped <- df[drop_ix, ]
        df <- df[!drop_ix, ]
        stopifnot(dplyr::n_distinct(df$patient) == n_completed)
        stopifnot(dplyr::n_distinct(df_dropped$patient) == n_dropped)
      }
    }
  )

  

  output <- list(
    df = df,
    n_completed = n_completed,
    df_dropped = df_dropped,
    n_dropped = n_dropped,
    n_patients = n_patients,
    n_periods = n_periods,
    sd_epsilon = sd_epsilon,
    sd1 = sd1,
    sd2 = sd2,
    a0 = a0,
    b0 = b0,
    t1 = t1,
    t2 = t2,
    g_j = g_j,
    possible_sequences = possible_sequences,
    .seed = .seed
  )
  return(output)
}


run_primary_analysis <- function(
  simulation_output,
  summarise_random_effects = TRUE,
  return_fit_object = TRUE
) {
  # patient-by-treatment interaction
  .fit_full1 <- lmerTest::lmer(
    iauc ~ 1 + diet + period + (1 + diet | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lme4::lmerControl(optCtrl = list(maxfn = 500))
  )
  .fit_null1 <- lmerTest::lmer(
    iauc ~ 1 + diet + period + (1 | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lme4::lmerControl(optCtrl = list(maxfn = 500))
  )
  patient_pval <- anova(.fit_full1, .fit_null1)[["Pr(>Chisq)"]][2]
  # microbiome-by-treatment interaction
  .fit_full2 <- lmerTest::lmer(
    iauc ~ 1 + m + diet + period + m:diet + (1 + diet | patient),
    data = simulation_output$df,
    REML = TRUE,  # used for estimation and figures
    control = lme4::lmerControl(optCtrl = list(maxfn = 200))
  )
  .fit_null2 <- lmerTest::lmer(
    iauc ~ 1 + m + diet + period + (1 + diet | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lme4::lmerControl(optCtrl = list(maxfn = 200))
  )
  microbiome_pval <- anova(lme4::refitML(.fit_full2), .fit_null2)[["Pr(>Chisq)"]][2]

  # parameter estimates from full model
  fixed_effects <- lme4::fixef(.fit_full2)
  reneff_sds <- attr(lme4::VarCorr(.fit_full2)$patient, "stddev")
  
  # adjust p values for multiple comparisons
  padj <- p.adjust(
    c("patient" = patient_pval, "microbiome" = microbiome_pval),
    method = "holm"
  )
  analysis_summary <- list(
    patient_pval = patient_pval,
    patient_padj = padj[["patient"]],
    microbiome_pval = microbiome_pval,
    microbiome_padj = padj[["microbiome"]],
    a0 = fixed_effects[["(Intercept)"]],
    b0 = fixed_effects[["dietB"]],
    t1 = fixed_effects[["m"]],
    t2 = fixed_effects[["m:dietB"]],
    sd1 = reneff_sds[["(Intercept)"]],
    sd2 = reneff_sds[["dietB"]],
    sd_epsilon = sigma(.fit_full2)
  )

  if (isTRUE(summarise_random_effects)) {
    analysis_summary[["random_effects"]] <- tibble::as_tibble(
      lme4::ranef(.fit_full2)$patient,
      rownames = "patient"
    ) %>%
    dplyr::rename(a = `(Intercept)`, b = dietB)
  }

  if (isTRUE(return_fit_object)) {
    analysis_summary[["fit"]] <- .fit_full2
  }

  analysis_summary[["fit_convergence"]] <- get_convergence_status(.fit_full2)
  
  simulation_output[["analysis"]] <- analysis_summary
  return(simulation_output)
}

estimate_individual_effects <- function(
    simulation_output,
    m_type = c("average", "observed"),
    bootstrap_samples=200
  ) {

  m_type <- match.arg(m_type)
  if (m_type == "average") {
    m_values = mean(simulation_output$df$m)
  } else {
    stop("m_type='observed' not implemented.")
    m_values = simulation_output$df %>%
      dplyr::arrange(patient) %>% 
      dplyr::select(patient, m) %>% 
      unique() %>% 
      pull(m)
  }

  # get predicted individual effects 
  .get_preds <- function(myfit) {
    .pred_data <- expand.grid(
      patient = unique(simulation_output$df$patient),
      diet = unique(simulation_output$df$diet),
      period = factor(1, levels = levels(simulation_output$df$period)),
      m = mean(simulation_output$df$m)
    )
    .pred_data$y <- lme4:::predict.merMod(
      myfit,
      newdata = .pred_data,
      re.form = NULL
    )
    .pred_data %>% 
      tidyr::pivot_wider(
        names_from = diet,
        values_from = y
      ) %>% 
      dplyr::arrange(patient) %>% 
      dplyr::mutate(d = B - A) %>% 
      dplyr::pull(d)
  }
  bootstrap_predictions <- lme4::bootMer(
    simulation_output$analysis$fit,
    nsim = bootstrap_samples, 
    FUN = .get_preds, 
    ncpus = 1
  )
  
  #  summary statistics
  df_summaries <- simulation_output$df %>% 
    select(patient, diet, cycle, iauc) %>% 
    pivot_wider(names_from = diet, values_from = iauc) %>% 
    mutate(d = B - A) %>% 
    group_by(patient) %>% 
    summarise(
      avg_delta = mean(d),
      se_delta = sd(d)/sqrt(simulation_output$n_periods/2),
      lower_delta = avg_delta - qnorm(0.975)*se_delta,
      upper_delta = avg_delta + qnorm(0.975)*se_delta,
    )
  # individual effects from model predictions
  df_individual_effects <- tibble::tibble(
    patient = unique(simulation_output$df$patient),
    estimate = bootstrap_predictions$t0,
    se = apply(bootstrap_predictions$t, 2, sd),
    lower = estimate - qnorm(.975) * se,
    upper = estimate + qnorm(.975) * se
  )  %>% 
    dplyr::inner_join(
      simulation_output$df %>% 
        select(patient, sequence, m, u1, u2, a, b) %>% 
        unique(),
      by = "patient"
    ) %>% 
    dplyr::inner_join(
      df_summaries,
      by = "patient"
    )

  simulation_output[["individual_effects"]] <- df_individual_effects

  return(simulation_output)
}


plot_simulation <- function(
    simulation_output
) {
  df_individual_effects <- simulation_output$individual_effects %>% 
    mutate(
      # because we predicted at avg(microbiome score)
      true_individual_effect = b - m*simulation_output$t2 + mean(m)*simulation_output$t2
    )
  df <- simulation_output$df
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 14))
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
    geom_pointrange(
      aes(
        x = as.numeric(fct_reorder(factor(patient), estimate, max))-.5,
        y = avg_delta,
        ymin = lower_delta,
        ymax = upper_delta
        ),
      color = "#28B463",
    ) +
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
        iauc,
        mean
      )
    ) %>% 
    ggplot(aes(x, iauc)) +
    geom_boxplot(aes(color = paste0("Diet ", diet)), alpha = 0) +
    scale_color_manual(
      values = c(
        "Diet A" = "#1B9E77", "Diet B" = "#D95F02"
      )
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    theme(
      axis.text.x = element_blank(),
      legend.position = c(.1, .85)
    ) +
    labs(
      subtitle = "Observed iAUC values for each patient",
      x = "Individual patients",
      y = "iAUC (mmol*min/L)",
      color = NULL
    )
  plot_indiv_obs_deltas <- df %>% 
    select(patient, diet, cycle, iauc) %>% 
    pivot_wider(names_from = diet, values_from = iauc) %>% 
    mutate(d = B - A, patient = patient) %>% 
    ggplot(aes(fct_reorder(patient, d, median), d)) +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.3, lwd = 1) +
    geom_boxplot() +
    geom_jitter(width = .2) +
    theme(
      axis.text.x = element_blank()
    ) +
    labs(
      subtitle = "Observed differences in iAUC",
      x = "Individual patients",
      y = "Delta iAUC (mmol*min/L)\n(B-A)"
    )
  
  plot_outcome_distribution <- df %>% 
    mutate(treatment = paste0("Diet ", diet)) %>% 
    ggplot(aes(iauc)) +
    geom_histogram(bins = 30, fill = "steelblue") +
    facet_wrap(~diet) +
    labs(x = "Postprandial blood glucose (mg/dL)",
          y = "Frequency")
  
  plot_observed_differences_by_cycle <- df %>% 
    select(patient, diet, cycle, iauc) %>% 
    pivot_wider(names_from = diet, values_from = iauc) %>% 
    mutate(d = B - A, patient = patient) %>% 
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
      xlab = "Delta iAUC (B-A)", 
      ylab = "Delta iAUC (B-A)",
      progress = FALSE
    )
  
  plot_observed_values_by_cycle <- df %>% 
    select(patient, diet, cycle, iauc) %>% 
    pivot_wider(names_from = diet, values_from = iauc) %>% 
    mutate(d = B - A,
            cycle = paste0("Cycle", cycle),
            x = fct_reorder(patient, d, median)) %>% 
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
  

  .lim <- 50
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
      aes(y = b),
      color = "red"
    ) +
    labs(
      x = "Summary-based estimates",
      y = "Shrunk estimates"
    ) +
    coord_cartesian(ylim = c(-.lim, .lim), xlim = c(-.lim, .lim))
  
  plot_cond_effects <- marginaleffects::plot_cme(
    simulation_output$analysis$fit, 
    "diet", 
    "m",
    re.form = ~ 0
  )
  
  simulation_output[["plots"]] <- list(
    individual_effects = plot_indiv_effects,
    conditional_effects = plot_cond_effects,
    observed_differences = plot_indiv_obs_deltas,
    observed_values = plot_indiv_obs_values,
    observed_values_by_treatment = plot_outcome_distribution,
    observed_differences_by_cycle = plot_observed_differences_by_cycle,
    observed_values_by_cycle = plot_observed_values_by_cycle,
    shrinkage = plot_shrinkage
  )

  return(simulation_output)
}
