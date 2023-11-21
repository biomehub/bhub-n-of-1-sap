
library(targets)
tar_option_set(
  packages = c("tidyverse", "here", "furrr", "lmerTest")
)
tar_source()

overall_seed <- 21112023

list(
  tar_target(
    name = power_simulation,
    command = run_power_simulation(
      n_sim = 10,
      n_threads = 10,
      overall_seed = overall_seed,
      overwrite = FALSE
    )
  )
)