source("jls_functions.R")
load_dependencies()
# jags library must also be installed
# see https://mcmc-jags.sourceforge.io/

targets_filename <- file.path(".", "beetles-targets.csv.gz")

targets <- read_csv(targets_filename)


targets_ts <- targets %>% 
              as_tsibble(index = datetime, key = c(variable,site_id))

latency_date  <- Sys.Date() - months(12) 
forecast_date <- latency_date - months(3)

past <- targets_ts %>% 
        filter(datetime < forecast_date) |>
        pivot_wider(names_from = "variable", values_from = "observation")

future <- targets_ts %>% 
          filter(datetime > forecast_date & datetime < latency_date) |>
          pivot_wider(names_from = "variable", values_from = "observation")


null_abundance <- fit_runjags_null_abundance(past   = past,
                                             future = future)
null_richness  <- fit_runjags_null_richness(past    = past,
                                             future = future)


