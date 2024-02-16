library(tidyverse)
library(tsibble)
library(neonstore)
library(neon4cast)
library(vars)
library(fpp)

# Retrieving the original targets data
if(!file.exists("beetles-targets.csv.gz")) {
  download.file("https://data.ecoforecast.org/neon4cast-targets/beetles/beetles-targets.csv.gz",
  "beetles-targets.csv.gz")
} else {
  targets <-  read_csv("beetles-targets.csv.gz")  
}

# Converting the targets csv to a time series tibble with abundance data first and followed by richness
targets_ts <- targets |> 
  as_tsibble(index = datetime, key = c(variable,site_id))

forecast_date <- Sys.Date() - months(3) # Disregarding the last 3 months of data.

past <-  targets_ts |>
  filter(datetime < forecast_date)  |>
  pivot_wider(names_from="variable", values_from="observation")

## access precip and temp data
sites <- unique(past$site_id)

# average over ensembles.  A little slow
if(!file.exists("historical_data.rds")) {
historical <- noaa_stage3() |> 
  filter(site_id %in% sites,
         variable %in% c("air_temperature", "precipitation_flux")) |>
  mutate(date = lubridate::as_date(datetime)) |>
  group_by(site_id, date, variable) |>
  summarise(mean = mean(prediction),
            min = min(prediction),
            max = max(prediction)) |>
  collect()

saveRDS(historical, file= "historical_data.rds")
} else {
  historical <- readRDS("historical_data.rds")
}

histroical_new <- historical |> 
  dplyr::select(-c(min, max)) |> 
  pivot_wider(names_from = "variable", values_from = "mean") |> 
     mutate(temp = air_temperature-273)

site = "OSBS"
df <- past |> 
  filter(site_id == site) |> 
  dplyr::left_join(
    histroical_new %>%
      filter(site_id == site), 
    by = c("datetime" = "date", "site_id" = "site_id"))

df_final <- df |> 
  filter(datetime >= "2021-03-22") |> 
  as_tibble() |> 
  dplyr::select(-c(richness, trapnights, air_temperature, datetime, site_id, iso_week))

plot(ts(df_final))
  
##Stationarity
##  H0: Series is Stationary

kpss.test(df_final$abundance,null="Level") # p-value=0.1 > 0.05. Stationary
kpss.test(df_final$precipitation_flux,null="Level") # p-value=0.1 > 0.05. Stationary 
kpss.test(df_final$temp,null="Level") # p-value=0.1 > 0.05. Stationary

##P-value=0.1<0.05, Thus H0 is not rejected. 
# All the series are stationary, we can proceed in fitting a VAR model. 


##VAR Model

df_final <- ts(df_final)

##Check for suitable number of lags

Model_1_lags=VARselect(df_final,lag.max = 8,type="none")
Model_1_lags   # Use no. of lags = 1

##Fit Model 

VARmodel1=VAR(df_final,p=1,type="none")
VARmodel1

##Check for suitable number of lags

Model_2_lags=VARselect(df_final,lag.max = 8,type="const")
Model_2_lags   # Use no. of lags = 1

##Fit Model
VARmodel2=VAR(df_final,p=1,type="const")

  
##Serial Correlation
## H0: Errors are uncorrelated

serial_VARmodel1=serial.test(VARmodel1,lags.pt = 12)
serial_VARmodel1   #p-value=0.99>0.05, Thus H0 is not rejected. Errors are not correlated. 

serial_VARmodel2=serial.test(VARmodel2,lags.pt = 12)
serial_VARmodel2   #p-value=0.99>0.05, Thus H0 is not rejected. Errors are not correlated. 

# Now we have 2 models with and without constant with errors uncorrelated

##Check for AIC values

AICmodel1=AIC(VARmodel1)
AICmodel1  #-414.1924

AICmodel2=AIC(VARmodel2)
AICmodel2  #-422.6987

###Lowest AIC is given by the model that includes the constant.


##Check for the significance of the parameters

summary(VARmodel2)
##But, Constant is significant in only one of the equations in the system. 
##Thus, it is not desirable to include the constant term. 
##so, we reject this model eventhough the AIC is the lowest

summary(VARmodel1)
##Significance of the parameters is somewhat better, so we select VARmodel1.


##Check for the roots
## Roots<1, for the equation to be stable. 
# (System should behave equally at different time points, then only the model can be used for long term predictions)

roots(VARmodel1)
##All are less than 1. Thus, a stable model. Can proceed with this. 


##Thus we select VARmodel1


##Obtain model coefficients
## We only look at the first equation as we are interested in modelling the abundance
VARmodel1

##Forecast

VARmodel1_forecast=forecast(VARmodel1, h = 10)
plot(VARmodel1_forecast$forecast$abundance)
  
