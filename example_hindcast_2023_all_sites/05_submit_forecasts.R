# Author info ----
# Eric R Sokol (esokol@battelleecology.org)
# Example for EFI RCN NEON Forecasting Challenge
# Populations and communities, Ground Beetles theme
# for a workshop at ESA 2024 meeting




# libraries ----
library(tidyverse)
library(fable)

# helper functions for format for EFI submission
library(neon4cast) # remotes::install_github("eco4cast/neon4cast", dep=TRUE)

mod_dir <- "fc_efi_null/"
mod_files <- list.files(mod_dir)


# submit the file
forecast_file <- paste0(mod_dir,"beetles-2024-02-25-bet_example_mod_naive.csv.gz")
neon4cast::submit(forecast_file = forecast_file)

forecast_file <- paste0(mod_dir,"beetles-2024-02-25-bet_example_mod_null.csv.gz")
neon4cast::submit(forecast_file = forecast_file)

  
fc_mod <- read_csv(forecast_file)
fc_mod %>%
  filter(site_id %in% c("YELL","OSBS","HARV","BARR")) %>%
  as_tsibble(index = datetime,
             key = c(model_id, parameter, site_id)) %>%
  ggplot(aes(datetime, prediction, color = parameter)) +
  geom_line() +
  facet_grid(site_id ~ model_id, scales = "free") +
  theme(legend.position = "none")







mod_dir <- "fc_efi_lm_ensemble/"
mod_files <- list.files(mod_dir)


# submit the file
forecast_file <- paste0(mod_dir,"beetles-2024-02-25-bet_lm_precip_clim_ensemble.csv.gz")
neon4cast::submit(forecast_file = forecast_file)

forecast_file <- paste0(mod_dir,"beetles-2024-02-25-bet_lm_temp_clim_ensemble.csv.gz")
neon4cast::submit(forecast_file = forecast_file)

forecast_file <- paste0(mod_dir,"beetles-2024-02-25-bet_lm_temp_precip_clim_ensemble.csv.gz")
neon4cast::submit(forecast_file = forecast_file)






