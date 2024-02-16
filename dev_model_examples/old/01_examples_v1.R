##############################################################################
# choose site ----
my_site = "OSBS"


##############################################################################
# set forecast date ----

forecast_date <- "2022-01-01" #fit up through 2021, forecast 2022 data
forecast_enddate <- "2025-12-31"


##############################################################################
# libraries ----
library(tidyverse)
library(fable)

# # get past weather and climate forecasts
# library(RopenMeteo) # remotes::install_github("FLARE-forecast/RopenMeteo")

# helper functions for format for EFI submission
library(neon4cast) # remotes::install_github("eco4cast/neon4cast", dep=TRUE)


##############################################################################
# Get data ----

# get NEON site metadata file 
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")

# # look at lat long for target site
# neon_site_info %>% filter(field_site_id == my_site) %>% select(field_site_id, field_latitude, field_longitude) %>%
#   as.data.frame()

# get past weather for my_site 
# past_climate <- RopenMeteo::get_historical_weather(
#   latitude = 29.68928,
#   longitude = -81.99343,
#   start_date = "2012-01-01",end_date = forecast_date,
#   site_id = my_site,variables = c("temperature_2m","precipitation"))
# write_csv(past_climate, "clim_dat_PAST_OSBS_test.csv")
past_climate_long <- read_csv("clim_dat_PAST_OSBS_test.csv") 

past_climate_wide <- past_climate_long %>%
  mutate(ymd = lubridate::as_date(datetime)) %>%
  select(datetime, ymd, site_id, 
         variable, prediction) %>%
  pivot_wider(names_from = variable, values_from = prediction) %>%
  group_by(ymd, site_id) %>%
  summarize(temperature_2m_mean = mean(temperature_2m),
            precipitation_sum = sum(precipitation)) %>%
  rename(datetime = ymd)

# read target data
# # Retrieving the original targets data
# targets <- read_csv("https://data.ecoforecast.org/neon4cast-targets/beetles/beetles-targets.csv.gz")
# write_csv(targets, "beetles-targets.csv")
targets <-  read_csv("beetles-targets.csv")  

# filter to my site
targets <- targets %>% dplyr::filter(site_id == my_site)

targets_future <- targets %>%
  filter(datetime >= forecast_date,
         datetime < "2023-01-01")

targets <- targets %>%
  filter(datetime < forecast_date)



# Future climate
# climate drivers - only have clim data for OSBS at the moment
# only predicting a through 2025
future_climate_long <- read_csv("clim_dat_OSBS_test.csv") 

future <- future_climate_long %>%
  filter(datetime >= forecast_date,
         datetime < forecast_enddate) %>% #forecast through end of 2025 
  select(-unit) %>%
  pivot_wider(names_from = variable, values_from = prediction)


##############################################################################
# format data
targets_ts <- targets %>%
  as_tsibble(index = datetime, key = c(variable,site_id))

targets_future_ts <- targets_future %>%
  as_tsibble(index = datetime, key = c(variable,site_id))

# future_clim_ts <- future_climate %>%
#   as_tsibble(index = datetime, key = c(variable, model_id, site_id))
# 
# past_clim_ts <- past_climate %>%
#   as_tsibble(index = datetime, key = c(variable, model_id, site_id))


##############################################################################
# visualize data

# visualize targets
targets_ts %>% fabletools::autoplot() +
  facet_grid(variable ~ .,scales = "free_y")

targets_future_ts %>% fabletools::autoplot() +
  facet_grid(variable ~ .,scales = "free_y")

# # visualize climate data
# past_clim_ts %>% fabletools::autoplot() +
#   facet_grid(variable ~ .,scales = "free_y")
# 
# future_clim_ts %>% fabletools::autoplot() +
#   facet_grid(variable ~ .,scales = "free_y")


# combine data
past <- targets_ts %>%
  pivot_wider(names_from = variable, values_from = observation) %>%
  left_join(
    past_climate_wide) %>%
  mutate(
    model_id = "past"
  )

##############################################################################
# simple linear model 

#Fit linear model based on past data: water temperature = m * air temperature + b
past$abundance_log <- log(past$abundance + 0.001)
past$temp_scaled <- scale(past$temperature_2m_mean)
past$precip_scaled <- scale(past$precipitation_sum)

future$temp_scaled <- scale(future$temperature_2m_mean)
future$precip_scaled <- scale(future$precipitation_sum)

fit <- lm(past$abundance_log ~ past$temp_scaled + past$precip_scaled)

# use linear regression to forecast water temperature for each ensemble member
future[,"abundance_log"] <- fit$coefficients[1] + 
  fit$coefficients[2] * future$temp_scaled + 
  fit$coefficients[3] * future$precip_scaled

future[,"abundance"] <- exp(future$abundance_log) + 0.001

all_data <- bind_rows(as_tibble(past), future)

all_data %>% select(model_id, abundance) %>% group_by(model_id) %>%
  summarize(
    mean_abund = mean(abundance), 
    min_abund = min(abundance),
    max_abund = max(abundance))

all_abund_ts <- all_data %>%
  select(datetime, site_id, model_id, abundance) %>%
  as_tsibble(index = datetime, key = c(site_id, model_id))

all_abund_ts %>% autoplot()



##############################################################################
## Compute a simple mean/sd model per site... obviously silly given huge seasonal aspect
null_abundance <- past  %>%
  fabletools::model(null = fable::MEAN(abundance)) %>%
  fabletools::forecast(h = "1 year") %>%
  neon4cast::efi_format()

# plot results
future_null <- null_abundance %>%
  pivot_wider(names_from = parameter, values_from = prediction) %>%
  as_tsibble(index = datetime, key = c(variable,site_id))

data_plot_null <- past %>% 
  mutate(model_id = "input_data") %>%
  bind_rows(future_null)

data_plot_null %>% ggplot(aes(datetime, abundance, colorConverter(modelid))) +
  geom_point() +
  geom_line(aes(datetime, mu), color = "blue") +
  geom_ribbon(aes(ymin = mu-sigma, ymax = mu+sigma), 
              color = "blue", fill = "blue", alpha = .5)

past_gap_filled <- past %>% fill_gaps()

arima_abundance <- past_gap_filled  %>%
  fabletools::model(arima = fable::ARIMA(abundance)) %>%
  fabletools::forecast(h = "1 year") %>%
  neon4cast::efi_format()

# time series linear model
fable::TSLM()