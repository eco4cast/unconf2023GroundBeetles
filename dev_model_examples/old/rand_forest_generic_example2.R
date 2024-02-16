# Load required libraries
library(tidymodels)
library(fable)
library(rsample)  # For rolling_origin()

# Example dataset with irregularly spaced time intervals
# Assuming 'abundance' is the response variable and 'predictor1', 'predictor2', etc. are predictor variables
# 'time' represents the time of measurement

# Generate example data (replace this with your own dataset)
set.seed(123)
n <- 1000
time <- sample(seq(as.Date('2020-01-01'), as.Date('2022-12-31'), by="day"), n, replace=FALSE)
predictor1 <- rnorm(n)
predictor2 <- rnorm(n)
abundance <- predictor1 + predictor2 + rnorm(n)
data <- data.frame(time, predictor1, predictor2, abundance)

# Convert 'time' to tsibble object
data_ts <- as_tsibble(data, index = time)
data_ts <- as_tsibble(data, key = "time")

# Define rolling origin for time series cross-validation
rolling_origin_spec <- rolling_origin(data_ts, initial = 365, assess = 14)  # 1 year for initial training, 14 days for testing

# ARIMA modeling specification
arima_spec <- ARIMA(
  log(abundance) ~ 1 + ARIMA(p = 1, d = 1, q = 1), 
  seasonal = ARIMA(P = 0, D = 1, Q = 0), stepwise = FALSE)

# Fit ARIMA model using time series cross-validation
arima_fit <- arima_spec %>%
  forecast(rolling_origin_spec, h = 14)  # Forecast horizon: two weeks

# Summarize the uncertainty around the prediction
summary(arima_fit)
