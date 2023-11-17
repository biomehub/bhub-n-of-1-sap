plot_simulation <- function(
    simulation_output
) {
  df_individual_effects <- simulation_output$individual_effects %>% 
    mutate(
      # because we predicted at avg(microbiome score)
      true_individual_effect = b - m*simulation_output$t2
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