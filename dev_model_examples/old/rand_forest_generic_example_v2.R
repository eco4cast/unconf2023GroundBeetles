# Load required libraries
library(tidyverse)
library(randomForest)
library(tidymodels)


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

# Extract time-based features
data$year <- as.numeric(format(data$time, "%Y"))
data$month <- as.numeric(format(data$time, "%m"))  # Month of the year
data$day_of_week <- as.numeric(format(data$time, "%u"))  # Day of the week
data$time_of_day <- as.numeric(format(data$time, "%H"))  # Time of day

# Convert data to tibble
data <- as_tibble(data)

# Define the recipe -- package recipes
rf_recipe <- recipe(abundance ~ ., data = data) %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal(), -all_outcomes())

# Define the random forest model specification -- package parsnip
rf_spec <- rand_forest(mode = "regression") %>%
  set_engine("randomForest", importance = TRUE)

# Define the workflow -- package workflows
rf_workflow <- workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(rf_spec)

# Fit the workflow -- fit method from workflows package
rf_fit <- rf_workflow %>%
  fit(data)


# Generate predictions for one month out with estimates of uncertainty
forecast_data <- tibble(
  time = seq(as.Date('2023-01-01'), by = "day", length.out = 30),
  predictor1 = rnorm(30),
  predictor2 = rnorm(30),
  year = as.numeric(format(time, "%Y")),
  month = as.numeric(format(time, "%m")),
  day_of_week = as.numeric(format(time, "%u")),
  time_of_day = as.numeric(format(time, "%H"))
)

# Make predictions -- method for predict in tsibble? uses new_data argument
forecast_results <- rf_fit %>%
  predict(new_data = forecast_data) %>%
  bind_cols(forecast_data)

# Summarize predictions
summary(forecast_results)

library(ggplot2)

# Join forecast results with original data
plot_data <- forecast_results %>%
  bind_rows(data) %>%
  mutate(time = as.Date(time))  # Ensure time variable is recognized correctly


# Plot the forecast with uncertainty
ggplot(plot_data, aes(x = time)) +
  geom_line(aes(y = abundance), color = "blue", size = 1) +  # Actual values
  geom_line(aes(y = .pred), color = "red", linetype = "dashed", size = 1) +  # Predicted values
  # geom_ribbon(aes(ymin = .pred_lower, ymax = .pred_upper), fill = "grey", alpha = 0.3) +  # Confidence interval
  labs(x = "Time", y = "Abundance") +
  theme_minimal()


# Generate predictions from individual trees in the forest
tree_preds <- predict(rf_fit$fit, data, type = "response", predict.all = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(
  time = rep(data$time, each = ntree(rf_fit$fit)),
  abundance = rep(data$abundance, each = ntree(rf_fit$fit)),
  predicted = unlist(tree_preds)
)

# Plot the distribution of predictions
ggplot(plot_data, aes(x = time, y = predicted)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(y = abundance), color = "blue", size = 1) +  # Actual values
  labs(x = "Time", y = "Abundance") +
  theme_minimal()