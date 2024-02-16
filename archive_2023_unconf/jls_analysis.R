# The JAGS library is also required for running this code
# see https://mcmc-jags.sourceforge.io/ for installation directions.

source("jls_functions.R")
load_dependencies()

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

sites               <- unique(past$site_id)
nsites              <- length(sites)
null_abundance_fits <- named_null_list(sites)
null_richness_fits  <- named_null_list(sites)

# each site with a separate model set

for (i in 1:nsites) {

  past_i   <- past[past$site_id == sites[i], ]
  future_i <- future[future$site_id == sites[i], ]
  message("fitting abundance for site ", sites[i])
  null_abundance_fits[[i]] <- fit_runjags_null_abundance(past   = past_i,
                                                         future = future_i)
  message("fitting richness for site ", sites[i])
  null_richness_fits[[i]]  <- fit_runjags_null_richness(past    = past_i,
                                                        future = future_i)

}

