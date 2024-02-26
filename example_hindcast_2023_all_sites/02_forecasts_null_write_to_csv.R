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

# where to write forecast output
out_dir <- "fc_efi_null"
if(!dir.exists(out_dir)) dir.create(out_dir)





#
# Loop through sites ----
#
# choose site and forecast date ----

fc_null_efi_all_sites <- data.frame()

vars_to_keep <- c("i","vars_to_keep",ls()) %>% unique()

for(i in 1:length(neon_bet_site_list)){
# for(i in 1:3){
    
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
        bet_example_mod_null = fable::MEAN(log1p(abundance)),
        bet_example_mod_naive = fable::NAIVE(log1p(abundance))) # random walk model, requires gapfill
    
    # make a forecast
    fc_null <- mod_fits %>%
      fabletools::forecast(h = "3 years") 
    
    # # visualize the forecast
    # fc_null %>% 
    #   autoplot(targets_train) +
    #   facet_grid(.model ~ ., scales = "free_y")
    
    
    # format for submission to EFI
    # for non-normal distributions, efi_format function draws samples to create
    # n time series to provide an estimate of uncertainty
    fc_null_efi <- fc_null %>% 
      mutate(site_id = my_site) %>% #efi needs a NEON site ID
      neon4cast::efi_format() %>%
      mutate(
        project_id = "neon4cast",
        reference_datetime = forecast_date,
        duration = "P1W")
    
    # # visualize the EFI-formatted submission
    # fc_null_efi %>% 
    #   as_tsibble(index = datetime,
    #              key = c(model_id, parameter)) %>%
    #   ggplot(aes(datetime, prediction, color = parameter)) +
    #   geom_line() +
    #   facet_grid(model_id ~ .)
    
    fc_null_efi_all_sites <- fc_null_efi %>%
      as.data.frame() %>%
      bind_rows(fc_null_efi_all_sites)
  })
  
  message(i,"/",length(neon_bet_site_list))

}#End loop




# # visualize the EFI-formatted submission ----
# fc_null_efi_all_sites %>%
#   as_tsibble(index = datetime,
#              key = c(model_id, parameter, site_id)) %>%
#   ggplot(aes(datetime, prediction, color = parameter)) +
#   geom_line() +
#   facet_grid(model_id ~ site_id, scales = "free_y")



# formatting outfiles for submission ----
# need to separate by mod
# multiple sites can be in one submission file I think
#file name format is: theme_name-year-month-day-model_id.csv
theme_name <- "beetles"
file_date <- Sys.Date()

# write out the forecast files
for(i_model_id in fc_null_efi_all_sites$model_id %>% unique()){
  
  fc_to_write <- fc_null_efi_all_sites %>%
    dplyr::filter(model_id == i_model_id)
  
  forecast_file <- paste0(out_dir,"/",theme_name,"-",file_date,"-",i_model_id,".csv.gz")
  write_csv(fc_to_write, forecast_file)
  
  # # submit the file
  # neon4cast::submit(forecast_file = forecast_file)
}






