# Load packages ----

library(tidyverse)
library(RopenMeteo) # remotes::install_github("FLARE-forecast/RopenMeteo")




# Getting ERA5 historical weather data ----

# set end date for weather data download
my_start_date <- "2012-01-01"
my_end_date <- "2023-12-31"

# get NEON site metadata file 
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")

if(!dir.exists("historical_weather_data")) dir.create("historical_weather_data")

for(i in 1:nrow(neon_site_info)){
  i_site <- neon_site_info$field_site_id[i]
  
  i_file_name <- paste0("historical_weather_data/hist_weather_2012-2023_",i_site,".csv")
  
  if(file.exists(i_file_name)) next
  
  i_past_climate <- data.frame()
  try({
    # get past weather for my_site
    i_past_climate <- RopenMeteo::get_historical_weather(
      latitude = neon_site_info$field_latitude[i],
      longitude = neon_site_info$field_longitude[i],
      start_date = my_start_date,
      end_date = my_end_date,
      site_id = i_site,
      variables = c("temperature_2m","precipitation"))
  })
  
  if(nrow(i_past_climate) > 0) write_csv(i_past_climate, file = i_file_name)
  
  message("Waiting...")
  Sys.sleep(60*60) # slow down process to avoid hitting rate limit
  message(i_site,": ",i,"/",nrow(neon_site_info))
}




# Get climate model outputs for past and future ----

# set end date for weather data download
my_start_date <- "2012-01-01"
my_end_date <- "2050-12-31"

# get NEON site metadata file 
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")

if(!dir.exists("future_climate_data")) dir.create("future_climate_data")

climate_model_list <- c("CMCC_CM2_VHR4",
                        "FGOALS_f3_H",
                        "HiRAM_SIT_HR",
                        "MRI_AGCM3_2_S",
                        "EC_Earth3P_HR",
                        "MPI_ESM1_2_XR",
                        "NICAM16_8S")

for(i in 1:nrow(neon_site_info)){
  for(i_mod in climate_model_list){
    i_site <- neon_site_info$field_site_id[i]
    
    i_file_name <- paste0("future_climate_data/future_climate_2012-2050_",
                          i_site, "_",
                          i_mod,".csv")
    
    if(file.exists(i_file_name)) next
    
    i_future_climate <- data.frame()
    try({
      # get past weather for my_site
      i_future_climate <- RopenMeteo::get_climate_projections(
        latitude = neon_site_info$field_latitude[i],
        longitude = neon_site_info$field_longitude[i],
        start_date = my_start_date,
        end_date = my_end_date,
        site_id = i_site,
        model = i_mod,
        variables = c("temperature_2m_mean","precipitation_sum"))
    })
    
    if(nrow(i_future_climate) > 0){ 
      write_csv(i_future_climate, file = i_file_name)
      message(i_site,": ",i_mod)
    }
  }
  
  message("Waiting...")
  Sys.sleep(60*60) #slow down process to avoid rate limit
  message(i_site,": ",i,"/",nrow(neon_site_info))
}
