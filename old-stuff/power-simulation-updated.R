# metdata settings
.seed <- 1
n_patients <- 80 # 200
n_cycles <- 5 # 50
n_treats <- 2
N <- n_patients * n_treats * n_cycles
# parameters
alpha_0 <- 60
tau_1 <- 0
beta_0 <- 10
tau_2 <- 0
sigma_1 <- 13
sigma_2 <- 4
sigma_12 <- 0
Sigma_u <- matrix(c(sigma_1^2, sigma_12, sigma_12, sigma_2^2), ncol = 2)
sigma_epsilon <- 10

df <- expand.grid(
        patient = 1:n_patients,
        cycle = 1:n_cycles,
        diet = LETTERS[1:n_treats]
      ) %>%
      arrange(patient, cycle, diet) %>%
  # add patient-specific terms (Microbiome covariate and random effects)
  group_by(patient) %>%
  mutate(
    # microbiome covariate for HTE testing
    microbiome = rnorm(1, mean = 0, sd = 1),
    # random intercepts and slopes
    u = list(
      set_names(
        MASS::mvrnorm(1, c(0, 0), Sigma_u),
        c('u1', 'u2')
      )
    )
  ) %>%
  ungroup() %>%
  unnest_wider(u) %>%
  mutate(
    alpha_i = alpha_0 + tau_1 * microbiome + u1,
    beta_i = beta_0 + tau_2 * microbiome + u2,
    epsilon = rnorm(nrow(.), mean = 0, sd = sigma_epsilon),
    y = alpha_i + beta_i * as.numeric(diet == "B") + epsilon
  ) %>%
  dplyr::select(
    patient, cycle, diet, y, microbiome, alpha_i, beta_i, u1, u2, epsilon
  )
df_pct <- df %>%
  dplyr::select(patient, microbiome, alpha_i, beta_i, u1, u2)  %>% unique

f <- lmer(
  y ~ 1 + microbiome + diet + microbiome:diet + (1 + diet | patient), data = df
)
summary(f)

df2 <- withr::with_seed(
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