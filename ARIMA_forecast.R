
library(tidyverse)
library(fable)

# Retrieving the original targets data
download.file("https://data.ecoforecast.org/neon4cast-targets/beetles/beetles-targets.csv.gz",
              "beetles-targets.csv.gz")
targets <-  read_csv("beetles-targets.csv.gz")

# Converting the targets csv to a time series tibble with abundance data first and followed by richness
targets_ts <- targets |> 
  as_tsibble(index = datetime, key = c(variable,site_id))

forecast_date <- Sys.Date() - months(15) # Disregarding the last 15 months of data. 

past <-  targets_ts |>
  filter(datetime < forecast_date)  |>
  pivot_wider(names_from="variable", values_from="observation")

gap_filled <- tsibble::fill_gaps(past) # Fills in time gaps with NAs as we need equal time intervals to work with ARIMA.


# Forecast
arima_richness <- gap_filled  |> 
  model(arima = ARIMA(richness)) |>
  forecast(h = "1 year")

arima_abundance <- gap_filled  |>
  model(arima = ARIMA(abundance)) |>
  forecast(h = "1 year") 


# Plot the forecast only for a particular site (eg: BARR)
arima_richness |> 
  filter(site_id == "BARR") |> 
  autoplot()

arima_abundance |> 
  filter(site_id == "BARR") |> 
  autoplot()
