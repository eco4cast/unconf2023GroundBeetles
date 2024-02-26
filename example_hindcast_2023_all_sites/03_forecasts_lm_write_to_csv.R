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




# get targets ----
# beetle targets are here
# url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1W/beetles-targets.csv.gz"

# read in the table
# bet_targets_all <- read_csv(url)
# write_csv(bet_targets_all, "bet_targets_all.csv")
bet_targets_all <- read_csv("bet_targets_all.csv")

# get list of sites available for bet forecasting
neon_bet_site_list <- bet_targets_all$site_id %>% unique()




# Set forecast start and end dates
forecast_date <- "2022-01-01" #fit up through 2021, forecast 2022 data
forecast_enddate <- "2025-01-01"

# We will only use output from one climate model for now
clim_model_id <- "CMCC_CM2_VHR4"

# where to write forecast output
out_dir <- "fc_efi_lm_CMCC_CM2_VHR4"
if(!dir.exists(out_dir)) dir.create(out_dir)





#
# Loop through sites ----
#
# choose site and forecast date ----

fc_lm_efi_all_sites <- data.frame()
mod_fit_stats <- data.frame()

vars_to_keep <- c("i","vars_to_keep","mod_fit_stats",
                  ls()) %>% unique()

for(i in 1:length(neon_bet_site_list)){
  # for(i in 1:3){
  # i <- 1
  # reset for next iteration of loop
  vars_to_rm <- ls() %>% dplyr::setdiff(vars_to_keep)
  rm(vars_to_rm)
  
  try({
    # get site id
    my_site = neon_bet_site_list[i]
    
    # filter to targets and format data for model fitting ----
    # read in and filter to my site
    targets <- bet_targets_all %>% 
      mutate(datetime = lubridate::as_date(datetime)) %>%
      dplyr::filter(site_id == my_site,
                    datetime < "2022-12-31") # excluding provisional data 
    
    # # visualize the targets dataset
    # targets %>% 
    #   as_tsibble(index = datetime, key = variable) %>%
    #   autoplot() +
    #   facet_grid(variable ~ ., scales = "free_y")
    
    # create a tsibble training dataset from targets
    targets_train <- targets %>%
      filter(datetime < forecast_date) %>%
      pivot_wider(names_from = variable, values_from = observation) %>%
      as_tsibble(index = datetime)
    
    
    # Get climate data for site ----
    path_to_clim_data <- "C:/Users/esokol/Box/00_MY_NEON/Forecasting_Beetles/future_climate_data/"
    
    # list files
    clim_file_list <- list.files(path_to_clim_data)
    
    # filter to files for target site
    clim_files_to_stack <- clim_file_list[grepl(my_site, clim_file_list) &
                                            grepl(clim_model_id, clim_file_list)]
    
    # loop through to stack the files
    clim_long <- data.frame()
    for(i_file in clim_files_to_stack){
      clim_long <- 
        bind_rows(
          clim_long,
          read_csv(paste0(path_to_clim_data,i_file)) %>%
            filter(datetime <= forecast_enddate))
    }
    
    # # make a tsibble object
    # clim_long_ts <- clim_long %>%
    #   as_tsibble(index = datetime, 
    #              key = c(variable, model_id))
    
    # make wide
    clim_wide <- clim_long %>%
      select(-unit) %>%
      pivot_wider(names_from = variable, values_from = prediction)
    
    # # visualize climate data
    # clim_long_ts %>%
    #   ggplot(aes(datetime, prediction, color = model_id)) +
    #   geom_line() +
    #   facet_grid(variable ~ ., scales = "free_y") +
    #   geom_vline(xintercept = lubridate::as_date(forecast_date),
    #              lty = 2)
    
    
    # Simple linear model ----
    
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
    mod_fits <- targets_clim_train %>%
      fabletools::model(
        lm_temp = fable::TSLM(log1p(abundance) ~ temperature_2m_mean),
        lm_precip = fable::TSLM(log1p(abundance) ~ precipitation_sum),
        lm_temp_precip = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + precipitation_sum))
    # mod_temp_season = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + season(period = "1 year")),
    # mod_precip_season = fable::TSLM(log1p(abundance) ~ precipitation_sum + season(period = "1 year")),
    # mod_season = fable::TSLM(log1p(abundance) ~ season(period = "1 year")))
    
    # fit stats
    mod_fit_stats <- fabletools::report(mod_fits) %>% 
      as.data.frame() %>%
      mutate(site_id = my_site) %>%
      bind_rows(mod_fit_stats)
    
    # # visualize model fit
    # # augment reformats model output into a tsibble for easier plotting
    # fabletools::augment(mod_fits) %>%
    #   ggplot(aes(x = datetime)) +
    #   geom_line(aes(y = abundance, color = "Obs")) +
    #   geom_line(aes(y = .fitted, color = .model))
    # 
    # # plot observed vs. fitted
    # fabletools::augment(mod_fits) %>%
    #   ggplot(aes(x = abundance, y = .fitted)) +
    #   geom_point(aes(color = .model)) +
    #   geom_abline(slope = 1) +
    #   xlim(0,.35) + ylim(0,.35)
    
    
    # make a forecast
    # filter "future" climate data to target climate model
    fc_mods <- mod_fits %>%
      fabletools::forecast(
        new_data = 
          clim_future %>%
          dplyr::filter(model_id == clim_model_id) %>%
          select(-model_id) %>%
          as_tsibble(index = datetime)) 
    
    # # visualize the forecast
    # fc_mods %>% 
    #   autoplot(targets_clim_train) +
    #   facet_grid(.model ~ .)
    
    # format for submission to EFI
    # for non-normal distributions, efi_format function draws samples to create
    # n time series to provide an estimate of uncertainty
    # https://projects.ecoforecast.org/neon4cast-ci/instructions.html
    # I'm putting "example" in the name so the model does not register as 
    # an official entry to the challenge
    
    # # update model name for submission
    # efi_model_id <- paste0("bet_example_lm_temp2m_",clim_model_id)
    
    # update dataframe of model output for submission
    fc_climate_mods_efi <- fc_mods %>% 
      # mutate(.model = efi_model_id) %>%
      mutate(site_id = my_site) %>% #efi needs a NEON site ID
      neon4cast::efi_format() %>%
      mutate(
        project_id = "neon4cast",
        model_id = paste0("bet_",model_id, "_",clim_model_id),
        reference_datetime = forecast_date,
        duration = "P1W")
    
    # # visualize the EFI-formatted submission
    # fc_climate_mods_efi %>% 
    #   as_tsibble(index = datetime,
    #              key = c(model_id, parameter)) %>%
    #   ggplot(aes(datetime, prediction, color = parameter)) +
    #   geom_line() +
    #   facet_grid(model_id ~ .) 
    
    
    fc_lm_efi_all_sites <- fc_climate_mods_efi %>%
      as.data.frame() %>%
      bind_rows(fc_lm_efi_all_sites, .)
    
  }) #END try
  
  message(i,"/",length(neon_bet_site_list))
  
}#END loop





# # visualize the EFI-formatted submission ----
# fc_lm_efi_all_sites %>%
#   as_tsibble(index = datetime,
#              key = c(model_id, parameter, site_id)) %>%
#   ggplot(aes(datetime, prediction, color = parameter)) +
#   geom_line() +
#   facet_grid(site_id ~ model_id, scales = "free")







# formatting outfiles for submission ----
# need to separate by mod
# multiple sites can be in one submission file I think
#file name format is: theme_name-year-month-day-model_id.csv
theme_name <- "beetles"
file_date <- Sys.Date()

# write out the forecast files
for(i_model_id in fc_lm_efi_all_sites$model_id %>% unique()){
  
  fc_to_write <- fc_lm_efi_all_sites %>%
    dplyr::filter(model_id == i_model_id)
  
  forecast_file <- paste0(out_dir,"/",theme_name,"-",file_date,"-",i_model_id,".csv.gz")
  write_csv(fc_to_write, forecast_file)
  
  # # submit the file
  # neon4cast::submit(forecast_file = forecast_file)
}







