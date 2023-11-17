get_n_of_1_dataset <- function(
    n_patients = 80,
    n_periods = 10,
    sd_epsilon = 10,
    sd1 = 20,
    sd2 = 5,
    a0 = 200,
    b0 = 0,
    t1 = 2,
    t2 = 2,
    g_j = rep(0, n_periods),
    ite_explainable_fraction = 0.5,
    possible_sequences = c("ABBABAABBA", "BAABABBAAB"),
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
      M <- rnorm(n = n_patients, sd = 1)
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
          df[ij_ix, "iauc"] <- a[i] + b[i] * x_ij + epsilon[i, j]
        }
      }

    }
  )

  output <- list(
    df = df,
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
    ite_explainable_fraction = ite_explainable_fraction,
    possible_sequences = possible_sequences,
    .seed = .seed
  )
  return(output)
}


run_primary_analyses <- function(
  simulation_output,
  summarise_random_effects = TRUE,
  return_fit_object = TRUE
) {
  # patient-by-treatment interaction
  .fit_full1 <- lmerTest::lmer(
    iauc ~ 1 + diet + (1 + diet | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lmerControl(optCtrl = list(maxfn = 200))
  )
  .fit_null1 <- lmerTest::lmer(
    iauc ~ 1 + diet + (1 | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lmerControl(optCtrl = list(maxfn = 200))
  )
  patient_pval <- anova(.fit_full1, .fit_null1)[["Pr(>Chisq)"]][2]
  # microbiome-by-treatment interaction
  .fit_full2 <- lmerTest::lmer(
    iauc ~ 1 + m + diet + m:diet + (1 + diet | patient),
    data = simulation_output$df,
    REML = TRUE,  # used for estimation and figures
    control = lmerControl(optCtrl = list(maxfn = 200))
  )
  .fit_null2 <- lmerTest::lmer(
    iauc ~ 1 + m + diet + (1 + diet | patient),
    data = simulation_output$df,
    REML = FALSE,
    control = lmerControl(optCtrl = list(maxfn = 200))
  )
  microbiome_pval <- anova(lme4::refitML(.fit_full2), .fit_null2)[["Pr(>Chisq)"]][2]

  # parameter estimates from full model
  fixed_effects <- lme4::fixef(.fit_full2)
  reneff_sds <- attr(lme4::VarCorr(.fit_full2)$patient, "stddev")
  analysis_summary <- list(
    patient_pval = patient_pval,
    microbiome_pval = microbiome_pval,
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
    analysis_summary[["fit"]] <- .fit_null2
  }
  
  simulation_output[["analysis"]] <- analysis_summary
  return(simulation_output)
}

estimate_individual_effects <- function(simulation_output, m_type = c("average", "observed")) {
  m_type <- match.arg(m_type)
  if (m_type == "average") {
    m_values = mean(simulation_output$df$m)
  } else {
    m_values = simulation_output$df %>%
      dplyr::arrange(patient) %>% 
      dplyr::select(patient, m) %>% 
      unique() %>% 
      pull(m)
  }
  # get predicted individual effects 
  bootstrap_predictions <- lme4::bootMer(
    simulation_output$analysis$fit,
    nsim = bootstrap_samples, 
    FUN = \(myfit) {
      .pred_data <- expand.grid(
        patient = levels(simulation_output$df$patient),
        diet = unique(simulation_output$df$diet)
      )
      .pred_data$m = m_values
      
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
    }, 
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
  tibble::tibble(
    patient = paste0("pt-", 1:simulation_output$n_patients),
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
