library(tidyverse)
library(RopenMeteo) # remotes::install_github("FLARE-forecast/RopenMeteo")

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
