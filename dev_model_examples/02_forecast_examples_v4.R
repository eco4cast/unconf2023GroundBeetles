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




# choose site and forecast date ----
my_site = "OSBS"

forecast_date <- "2022-01-01" #fit up through 2021, forecast 2022 data
forecast_enddate <- "2025-01-01"




# Get targets and format data for model fitting ----

# read target data
# # Retrieving the original targets data
# targets <- read_csv("https://data.ecoforecast.org/neon4cast-targets/beetles/beetles-targets.csv.gz")
# write_csv(targets, "beetles-targets.csv")

# read in and filter to my site
targets <-  read_csv("beetles-targets.csv") %>% 
  dplyr::filter(site_id == my_site,
                datetime < "2022-12-31") # excluding provisional data 

# visualize the targets dataset
targets %>% 
  as_tsibble(index = datetime, key = variable) %>%
  autoplot() +
  facet_grid(variable ~ ., scales = "free_y")

# create a tsibble training dataset from targets
targets_train <- targets %>%
  filter(datetime < forecast_date) %>%
  pivot_wider(names_from = variable, values_from = observation) %>%
  as_tsibble(index = datetime)




# Null models ----
# Compute a simple mean/sd model and a random walk model as null models

# specify and fit models
# Using a log(x + 1) transform on the abundance data
mod_fits <- targets_train %>% 
  tsibble::fill_gaps() %>%
  fabletools::model(
    mod_null = fable::MEAN(log1p(abundance)),
    mod_naive = fable::NAIVE(log1p(abundance))) # random walk model, requires gapfill

# make a forecast
fc_null <- mod_fits %>%
  fabletools::forecast(h = "3 years") 

# visualize the forecast
fc_null %>% 
  autoplot(targets_train) +
  facet_grid(.model ~ ., scales = "free_y")

# format for submission to EFI
# for non-normal distributions, efi_format function draws samples to create
# n time series to provide an estimate of uncertainty
fc_null_efi <- fc_null %>% 
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format()

# visualize the EFI-formatted submission
fc_null_efi %>% 
  as_tsibble(index = datetime,
             key = c(model_id, parameter)) %>%
  ggplot(aes(datetime, prediction, color = parameter)) +
  geom_line() +
  facet_grid(model_id ~ .)




# Using climate model outputs as predictors of beetle abundance ----

# Get climate data ----
path_to_clim_data <- "C:/Users/esokol/Box/00_MY_NEON/Forecasting_Beetles/future_climate_data/"

# list files
clim_file_list <- list.files(path_to_clim_data)

# filter to files for target site
clim_files_to_stack <- clim_file_list[grepl(my_site, clim_file_list)]

# loop through to stack the files
clim_long <- data.frame()
for(i_file in clim_files_to_stack){
  clim_long <- 
    bind_rows(
      clim_long,
      read_csv(paste0(path_to_clim_data,i_file)) %>%
        filter(datetime <= forecast_enddate))
}

# make a tsibble object
clim_long_ts <- clim_long %>%
  as_tsibble(index = datetime, 
             key = c(variable, model_id))

# make wide
clim_wide <- clim_long %>%
  select(-unit) %>%
  pivot_wider(names_from = variable, values_from = prediction)

# visualize climate data
clim_long_ts %>%
  ggplot(aes(datetime, prediction, color = model_id)) + 
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  geom_vline(xintercept = lubridate::as_date(forecast_date),
             lty = 2)




# Simple linear model ----

# We will only use one climate model out for now
clim_model_id <- "CMCC_CM2_VHR4"

# subset into past and future datasets, based on forecast_date
clim_past <- clim_wide %>%
  filter(model_id == clim_model_id,
         datetime < forecast_date,
         datetime > "2012-01-01")

clim_future <- clim_wide %>%
  filter(model_id == clim_model_id,
         datetime >= forecast_date,
         datetime <= forecast_enddate)

# combine target and climate data into a training dataset
targets_clim_train <- targets_train %>%
  left_join(clim_past)

# specify and fit model
mod_fit_candidates <- targets_clim_train %>%
  fabletools::model(
    mod_temp = fable::TSLM(log1p(abundance) ~ temperature_2m_mean),
    mod_precip = fable::TSLM(log1p(abundance) ~ precipitation_sum),
    mod_both = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + precipitation_sum),
    mod_temp_season = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + season(period = "1 year")),
    mod_precip_season = fable::TSLM(log1p(abundance) ~ precipitation_sum + season(period = "1 year")),
    mod_season = fable::TSLM(log1p(abundance) ~ season(period = "1 year")))

# look at fit stats
fabletools::report(mod_fit_candidates)

# visualize model fit
fabletools::augment(mod_fit_candidates) %>%
  ggplot(aes(x = datetime)) +
  geom_line(aes(y = abundance, color = "Obs")) +
  geom_line(aes(y = .fitted, color = .model))

# plot observed vs. fitted
fabletools::augment(mod_fit_candidates) %>%
  ggplot(aes(x = abundance, y = .fitted)) +
  geom_point(aes(color = .model)) +
  geom_abline(slope = 1) +
  xlim(0,.35) + ylim(0,.35)

# none of the models do particularly well

# focus on temperature model for now
mod_best_lm <- mod_fit_candidates %>% select(mod_temp)
report(mod_best_lm)

# make a forecast
fc_best_lm <- mod_best_lm %>%
  fabletools::forecast(
    new_data = 
      future_climate_wide %>%
      dplyr::filter(model_id == "CMCC_CM2_VHR4") %>%
      as_tsibble(index = datetime)) 

# visualize the forecast
fc_best_lm %>% 
  autoplot(past) +
  facet_grid(.model ~ .)

# format for submission to EFI
# for non-normal distributions, efi_format function draws samples to create
# n time series to provide an estimate of uncertainty
# https://projects.ecoforecast.org/neon4cast-ci/instructions.html
fc_climate_mods_efi <- fc_best_lm %>% 
  mutate(.model = paste0("bet_lm_temp2m_",.model)) %>%
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format() %>%
  mutate(
    project_id = "neon4cast",
    reference_datetime = forecast_date,
    duration = "P1W")

# visualize the EFI-formatted submission
fc_climate_mods_efi %>% 
  as_tsibble(index = datetime,
             key = c(model_id, parameter)) %>%
  ggplot(aes(datetime, prediction, color = parameter)) +
  geom_line() +
  facet_grid(model_id ~ .) 




# submit a forecast ----

# Before submitting, register your forecast model here https://forms.gle/kg2Vkpho9BoMXSy57
# registered theme_name is "beetles"

#file name format is: theme_name-year-month-day-model_id.csv
theme_name <- "beetles"
file_date <- Sys.Date()
model_id <- paste0("bet_lm_temp2m_",clim_model_id)

forecast_file <- paste0(theme_name,"-",file_date,"-",model_id,".csv.gz")
write_csv(fc_climate_mods_efi, forecast_file)

# submit the file
neon4cast::submit(forecast_file = forecast_file)
