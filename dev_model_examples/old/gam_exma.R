# Load required libraries
library(mgcv)

# Example dataset with irregularly spaced time intervals
# Assuming 'abundance' is the response variable and 'predictor1', 'predictor2', etc. are predictor variables
# 'time' represents the time of measurement

# Generate example data (replace this with your own dataset)
set.seed(123)
n <- 1000
time <- sample(seq(as.Date('2022-01-01'), as.Date('2022-12-31'), by="day"), n, replace=TRUE)
predictor1 <- rnorm(n)
predictor2 <- rnorm(n)
abundance <- predictor1 + predictor2 + rnorm(n)
data <- data.frame(time, predictor1, predictor2, abundance) %>% as.tibble()
data <- data %>%
  mutate(ntime = as.numeric(time))

# Fit GAM model
gam_model <- gam(abundance ~ s(ntime) + predictor1 + predictor2, data = data)

# Summary of the model
summary(gam_model)

# Plot smooth terms
plot(gam_model, select = 1)  # Plot smooth term for time
