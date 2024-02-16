##############################################################################
# libraries ----
library(tidyverse)
library(fable)

# helper functions for format for EFI submission
library(neon4cast) # remotes::install_github("eco4cast/neon4cast", dep=TRUE)


##############################################################################
# choose site and forecast date ----
my_site = "OSBS"

forecast_date <- "2022-01-01" #fit up through 2021, forecast 2022 data
forecast_enddate <- "2025-01-01"


##############################################################################
# Get targets and format data for model fitting ----

# read target data
# # Retrieving the original targets data
# targets <- read_csv("https://data.ecoforecast.org/neon4cast-targets/beetles/beetles-targets.csv.gz")
# write_csv(targets, "beetles-targets.csv")

# read in and filter to my site
targets <-  read_csv("beetles-targets.csv") %>% 
  dplyr::filter(site_id == my_site) %>%
  dplyr::select(-c(site_id, iso_week))

# combine past targets and climate data for model fitting
past <- targets %>%
  filter(datetime < forecast_date) %>%
  pivot_wider(names_from = variable, values_from = observation) %>%
  left_join(past_climate_wide) %>%
  as_tsibble(index = datetime)




##############################################################################
# set local paths to forecast data
path_to_past_weather_data <- "C:/Users/esokol/Box/00_MY_NEON/Forecasting_Beetles/historical_weather_data/"
path_to_future_clim_data <- "C:/Users/esokol/Box/00_MY_NEON/Forecasting_Beetles/future_climate_data/"



##############################################################################
# Get past climate data ----

# Past climate data
# list files
weather_file_list <- list.files(path_to_past_weather_data)
# filter to files for target site
weather_files_to_stack <- weather_file_list[grepl(my_site, weather_file_list)]
# loop through to stack the files, if necessary
past_climate_long <- data.frame()
for(i_file in weather_files_to_stack){
  past_climate_long <- 
    bind_rows(
      past_climate_long,
      read_csv(paste0(path_to_past_weather_data,
                      i_file)) %>%
        filter(datetime < forecast_date))
}

# Aggregate past climate data to daily to match forecast climate data
past_climate_wide <- past_climate_long %>%
  mutate(ymd = lubridate::as_date(datetime)) %>%
  select(datetime, ymd, model_id,
         variable, prediction) %>%
  pivot_wider(names_from = variable, values_from = prediction) %>%
  group_by(ymd, model_id) %>%
  summarize(temperature_2m_mean = mean(temperature_2m),
            precipitation_sum = sum(precipitation)) %>%
  rename(datetime = ymd)


##############################################################################
# Get future climate data
# list files
clim_file_list <- list.files(path_to_future_clim_data)
# filter to files for target site
clim_files_to_stack <- clim_file_list[grepl(my_site, clim_file_list)]
# loop through to stack the files
future_climate_long <- data.frame()
for(i_file in clim_files_to_stack){
  future_climate_long <- 
    bind_rows(
      future_climate_long,
      read_csv(paste0(path_to_future_clim_data,
                      i_file)) %>%
        filter(datetime >= forecast_date,
               datetime <= forecast_enddate))
}

# filter dataset to desired dates
future_climate_wide <- future_climate_long %>%
  filter(datetime >= forecast_date,
         datetime < forecast_enddate) %>% #forecast through end of 2025
  select(-c(site_id, unit)) %>%
  pivot_wider(names_from = variable, values_from = prediction)


##############################################################################
# visualize climate data

# stack past and future climate datasets for plotting
climate_ts <- 
  past_climate_wide %>%
  pivot_longer(cols = c(temperature_2m_mean, precipitation_sum),
               names_to = "variable",
               values_to = "value") %>%
  mutate(data_type = "past") %>%
  bind_rows(
    future_climate_long %>%
      mutate(data_type = "future") %>%
      rename(value = prediction)) %>%
  select(-unit) %>%
  as_tsibble(index = datetime,
             key = c(variable, model_id))

# plot climate time series
climate_ts %>%
  ggplot(aes(datetime, value, color = model_id)) + 
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y")


##############################################################################
# Null models
# Compute a simple mean/sd model per site... obviously silly given huge seasonal aspect

# specify and fit models
mod_fits <- past %>% 
  tsibble::fill_gaps() %>%
  fabletools::model(
    mod_null = fable::MEAN(log1p(abundance)),
    mod_naive = fable::NAIVE(log1p(abundance)))

# make a forecast
fc_null <- mod_fits %>%
  fabletools::forecast(h = "3 years") 

# visualize the forecast
fc_null %>% 
  autoplot(past) +
  facet_grid(.model ~ .)

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


##############################################################################
# Simple linear model ----

# specify and fit model
mod_fit_candidates <- past %>%
  fabletools::model(
    mod_temp = fable::TSLM(log1p(abundance) ~ temperature_2m_mean),
    mod_precip = fable::TSLM(log1p(abundance) ~ precipitation_sum),
    mod_both = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + precipitation_sum),
    mod_temp_season = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + season(period = "1 year")),
    mod_precip_season = fable::TSLM(log1p(abundance) ~ precipitation_sum + season(period = "1 year")),
    mod_season = fable::TSLM(log1p(abundance) ~ season(period = "1 year")))

# look at fit stats
fabletools::report(mod_fit_candidates)

# visualize model fit on past data
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

# focus on temperature model for now
mod_best_lm <- mod_fit_candidates %>% select(mod_temp)
report(mod_best_lm)


###################
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

###################
# Make "scenarios" based on each climate model forecast for future data
climate_model_scenarios <- fabletools::scenarios(
  CMCC_CM2_VHR4 = future_climate_wide %>% 
    dplyr::filter(model_id == "CMCC_CM2_VHR4") %>%
    as_tsibble(index = datetime),
  EC_Earth3P_HR = future_climate_wide %>% 
    dplyr::filter(model_id == "EC_Earth3P_HR") %>%
    as_tsibble(index = datetime),
  FGOALS_f3_H = future_climate_wide %>% 
    dplyr::filter(model_id == "FGOALS_f3_H") %>%
    as_tsibble(index = datetime),
  HiRAM_SIT_HR = future_climate_wide %>% 
    dplyr::filter(model_id == "HiRAM_SIT_HR") %>%
    as_tsibble(index = datetime),
  MPI_ESM1_2_XR = future_climate_wide %>% 
    dplyr::filter(model_id == "MPI_ESM1_2_XR") %>%
    as_tsibble(index = datetime),
  NICAM16_8S = future_climate_wide %>% 
    dplyr::filter(model_id == "NICAM16_8S") %>%
    as_tsibble(index = datetime))

# Use the mod_best_lm to forecast over each scenario
fc_climate_mods <- mod_best_lm %>%
  fabletools::forecast(new_data = climate_model_scenarios)

# visualize the forecasts based on different climate models
past %>%
  autoplot(abundance) +
  autolayer(fc_climate_mods) +
  facet_grid(.scenario ~ .)


# format for submission to EFI
# for non-normal distributions, efi_format function draws samples to create
# n time series to provide an estimate of uncertainty
fc_climate_mods_efi <- fc_climate_mods %>% 
  mutate(.model = paste0(.model,"_",model_id)) %>%
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format()

# visualize the EFI-formatted submission
fc_climate_mods_efi %>% 
  as_tsibble(index = datetime,
             key = c(model_id, parameter)) %>%
  ggplot(aes(datetime, prediction, color = parameter)) +
  geom_line() +
  facet_grid(model_id ~ .)


##############################################################################
# submit a forecast ----

# registered "team name" is "NEON Beetle Tutorial"
# registered theme is "Beetle Communities"

# Start by writing the forecast to file
file_date <- Sys.Date() #forecast$reference_datetime[1]
forecast_file <- paste0("beetles","-",file_date,"-",model_id,".csv.gz")
write_csv(forecast, forecast_file)

