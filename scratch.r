library(marginaleffects)
pred1 <- predictions(
    simulation_output$analysis$fit,
    newdata = datagrid(
        patient = unique(simulation_output$df$patient),
        diet = unique(simulation_output$df$diet),
        m = 0
    )
)
pred1 %>%
    ggplot(aes(patient, estimate, ymin = conf.low, ymax = conf.high, color = diet)) +
    geom_point(
        data = simulation_output$df,
        aes(x = patient, y = iauc, color = diet), alpha = .3, size = 1,
        inherit.aes = F,
        position = position_dodge(width = 0.15)
    ) +
    geom_pointrange(
        position = position_dodge(width = 0.15)
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

comp1 <- comparisons(
    simulation_output$analysis$fit,
    newdata = datagrid(
        patient = unique(simulation_output$df$patient),
        diet = unique(simulation_output$df$diet),
        m = 0
    ),
    vcov = "Kenward-Roger",
    variables = "diet",
    by = "patient"
)
comp1 %>%
    as_tibble() %>% 
    filter(term=="diet") %>% 
    #select(-diet, -rowid, -predicted) %>% 
    unique() %>% 
    inner_join(
        simulation_output$individual_effects %>% select(-estimate),
        by = "patient"
    ) %>% 
    ggplot(aes(patient, estimate, ymin = conf.low, ymax=conf.high)) +
    geom_pointrange() +
    geom_point(aes(y=b), color = "red")


n_pts <- 100; n_per_patient <- 5; N <- n_pts * n_per_patient
pts_ix <- rep(1:n_pts, each = n_per_patient)
random_interpt <- rnorm(n_pts)[pts_ix]
random_slope <- rnorm(n_pts)[pts_ix]
x <- sample(0:1, size = n_pts, replace = T)[pts_ix]
error <- rnorm(N)
df <- data.frame(
    y = 10 + random_interpt + (5 + random_slope) * x + error,
    treatment = x,
    subject = factor(pts_ix)
)
fit <- lmer(y ~ 1 + treatment + (1 + treatment | subject), data = df)
bootMer(fit, nsim = 200, FUN = \(.x) predict(.x), seed = 1)
comparisons(fit) %>% as_tibble
