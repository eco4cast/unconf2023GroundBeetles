# 1 About this tutorial

## 1.1 Learning objectives

-   Overview of the [Beetle
    Communities](https://projects.ecoforecast.org/neon4cast-docs/Beetles.html)
    theme for the [NEON Ecological Forecast
    Challenge](https://projects.ecoforecast.org/neon4cast-ci/)
-   How to create a simple forecast for the Beetle Communities theme.
-   How to submit/score a forecast to evaluate its accuracy.
-   How to use the NEON Forecast Challenge resources in your research
    and teaching.

## 1.2 Target user groups for this tutorial

This tutorial is intended to be used by ecological forecasters at any
stage of expertise and may be used as a learning tool as an introduction
to forecasting properties of ecological populations or communities.
Below, we provide code for introductory examples to walk through the
entire process of creating and submitting a forecast to the NEON
Ecological Forecasting challenge. This includes:

1.  Accessing target datasets of NEON ground beetle richness and
    abundance
2.  Accessing climate forecast data to use as drivers in models
    predicting beetle data
3.  How to use the `fable` package for R to specify and fit models
4.  How to submit a forecast to the forecast challenge

Upon completing this tutorial, participants should be able to create and
submit forecasts to the Beetle Communities theme of the EFI RCN NEON
Ecological Forecasting challenge.

## 1.3 Things you will need to complete this tutorial

You will need a current version of R (v4.2 or newer) to complete this
tutorial. We also recommend the RStudio IDE to work with R.

To complete the workshop via this markdown document the following
packages will need to be installed:

-   `tidyverse`
-   `lubridate`
-   `tsibble`
-   `fable`
-   `fabletools`
-   `remotes` (to install neon4cast from gitHub)
-   `neon4cast` (from github)

The following code chunk should be run to install packages.

``` r
install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
install.packages('lubridate') # working with dates and times
install.packages('tsibble') # working with timeseries data
install.packages('fable') # running forecasts
install.packages('fabletools') # helper functions for using fable
install.packages('remotes')
install.packages('tsibble') # package for dealing with time series data sets and tsibble objects
remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission
remotes::install_github("eco4cast/score4cast") # package to score forecasts
```

Then load the packages.

``` r
version$version.string
```

    ## [1] "R version 4.3.1 (2023-06-16 ucrt)"

``` r
library(tidyverse)
library(lubridate)
library(tsibble)
library(fable)
library(fabletools)
library(neon4cast)
library(score4cast)
```

# 2 Introduction

## 2.1 [The NEON Ecological Forecast Challenge](https://projects.ecoforecast.org/neon4cast-ci/)

The Challenge has been organized by the Ecological Forecasting
Initiative Research Coordination Network ([EFI
RCN](https://ecoforecast.org/)).

The Challenge asks the scientific community to produce ecological
forecasts of future observations of ecological data that will be
collected and published by the [National Ecological Observatory Network
(NEON)](https://www.neonscience.org/). The Challenge is split into [five
themes](https://projects.ecoforecast.org/neon4cast-ci/targets.html#sec-starting-sites)
that span aquatic and terrestrial systems, and population, community,
and ecosystem processes across a broad range of ecoregions. We are
excited to use this Challenge to learn more about the predictability of
ecological processes by forecasting NEON data before it is collected.

Which modeling frameworks, mechanistic processes, and statistical
approaches best capture community, population, and ecosystem dynamics?
These questions are answerable by a community generating a diverse array
of forecasts. The Challenge is open to any individual or team from
anywhere around the world that wants to submit forecasts. Learn more
about how you can participate
[here.](https://projects.ecoforecast.org/neon4cast-ci/instructions.html).

## 2.2 Goals for forecasts of ecological communities

Ecologists are interested in tracking changes in the **number of
individual organisms over time** (count data of abundance). Numbers of
each species will change due to births, deaths, and movement in
(immigration) or out (emigration) of populations. The abundance of an
ecological community is the sum of the number of individuals of each
species. For example, in a hypothetical community of beetles sampled in
Year 1 (time=1) , Species A has 10 individuals, Species B has 40
individuals and Species C has 50 individuals, giving a community
abundance of 100. In subsequent years the abundance may increase,
decrease or remain constant. A forecast may use this sequence of
observations over time to predict how many individuals will occur in the
next year (time=2), or a number of years into the future (time n). How
far into the future we predict is known as the forecast horizon. The
accuracy of the prediction is then compared to new observations and a
new prediction is made.

Ecologists are also interested in tracking changes in the **number of
species over time** (Species richness over time without species
identity) and in species turnover over time steps where the identity of
the species is known. In the example above there are three species
(A,B,C) but over time this can change if for example Species C was to
decline from 10 to zero individuals then the species richness would be
2, or increase if two previously unobserved species (D & E) arrive into
the community, the species richness would be 4. Note that the loss of
species A and the arrival of D and E gives a net species richness of 4,
but without keeping track of the identity of the species we may be
unaware of how community composition is changing over time.

Ecological communities change for many reasons and so it is important to
understand the drivers of changes in abundance or species richness by
adding environmental variables into the models. By knowing how species
change over time we can use the driving variables to predict, or
forecast, the values for the abundance and species richness variables
for the ecological communities into the future.

## 2.3 Overview of the [Beetle Communities](https://projects.ecoforecast.org/neon4cast-docs/Beetles.html) theme

**What**: Forecast abundance and/or richness of ground beetles
(carabids) collected in pitfall traps, standardized to sampling effort
(trap night). More information about the NEON data product
(DP1.10022.001, Ground beetles sampled from pitfall traps) we are
forecasting can be found
[here](https://data.neonscience.org/data-products/DP1.10022.001). Note
that we are not downloading the target dataset from the NEON data
portal. Rather, we will download a version of the dataset that has been
simplified and preformatted for this challenge by the EFI RCN.
Specifically, the targets are:

-   `abundance`: Total number of carabid individuals per trap-night,
    estimated each week of the year at each NEON site
-   `richness`: Total number of unique ‘species’ in a sampling bout for
    each NEON site each week.

**Where**: All 47 terrestrial [NEON
sites](https://www.neonscience.org/field-sites/explore-field-sites).

You can download metadata for all NEON sites as follows:

``` r
# To download the NEON site information table:
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")
```

This table has information about the field sites, including location,
ecoregion, and other useful metadata (e.g. elevation, mean annual
precipitation, temperature, and NLCD class).

Or, you can load a more targeted list of just the sites included in the
Beetle Communities theme:

``` r
site_data <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") %>%
  dplyr::filter(beetles == 1)
```

For this tutorial, we are going to focus on the [NEON site at
Ordway-Swisher Biological Station
(OSBS)](https://www.neonscience.org/field-sites/osbs), which is located
in Domain 03 (D03) in Florida.

**When**: Target data are available as early as 2013 at some sites, and
data are available at all sites from 2019 on. Because pitfall trap
samples need to be sorted and individuals counted and identified, the
latency for data publication can be nearly a year. In this tutorial we
will train our models on data from 2013-2021 and we will make forecasts
for the 2022 season so that we can score them immediately.

Check current ground beetle data availability on the [NEON Data
Portal](https://data.neonscience.org/data-products/DP1.10022.001#:~:text=last_page-,Availability%20and%20Download,-July%202013%20%E2%80%93%20January).

Determining a reference date (the first day of the forecast) and a
forecast horizon (the time window that is being forecast) are two major
challenges when forecasting population and community data. Other themes
in the NEON Ecological Forecast Challenge focus on targets that are
derived from instrument data, e.g., dissolved Oxygen or temperature in
lakes, that are collected at a high frequency (e.g., \> 1 Hz) and
available in near-real-time (e.g., latency \< 1 day). In practice,
forecasts for these types of data have horizons that are typically a
week to a month because they use output from weather forecasts (e.g.,
NOAA GEFS forecasts) as drivers in their models. These forecasts can be
evaluated and updated as new data roll in, and new weather forecasts are
published. This approach is known as “iterative near-term ecological
forecasting” (Dietze et al. 2018).

In contrast, population and community data are often available at a much
lower frequency (e.g., bi-weekly or annual sampling bouts) with a much
higher latency (e.g., 6 months to a year) because of the effort that is
required to collect and process samples and publish the data. Thus, the
goals and applications will likely be different for forecasts of these
types of data. There is still an opportunity to iterate, and update
forecasts, but over a much longer time period.

QUESTION: What are some use-cases for forecasts of ecological
populations and communities that you are interested in pursuing?

# 3 Forecasting NEON beetle communities

## 3.1 Define spatial and temporal parameters for our forecast

Here we set some values for variables in our code to identify the NEON
site we will be working with and the forecast start and end dates. This
will allow us to easily adapt this code for future runs at different
sites and forecasting different time windows.

Choose a NEON site:

``` r
# choose site
my_site = "OSBS"
```

Choose forecast start and end dates:

``` r
# date where we will start making predictions
forecast_startdate <- "2022-01-01" #fit up through 2021, forecast 2022 data

# date where we will stop making predictions
forecast_enddate <- "2025-01-01"
```

Note that the `forecast_startdate` will be renamed `reference_datetime`
when we submit our forecast to the The Challenge. As the parameter name
indicates, this date represents the beginning of the forecast. For this
tutorial, we are using a `forecast_startdate`, or `reference_datetime`,
that is in the past so that we can evaluate the accuracy of our
forecasts at the end of this tutorial.

The `forecast_enddate` is used to determine the forecast horizon. In
this example, we are setting a horizon to extend into the future.

If you want to create a true forecast to submit to the challenge, you
will want to set your `forecast_startdate` as today’s date. However, you
will likely need to wait until next year before you can evaluate your
forecast performance because you will need to wait for the data to be
collected, processed, and published.

## 3.2 Read in the data

We begin by first looking at the historic data - called the ‘targets’.
These data are available with a latency of approximately 330 days. Here
is how you read in the data from the targets file available from the EFI
server.

``` r
# beetle targets are here
url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1W/beetles-targets.csv.gz"

# read in the table
targets <- read_csv(url) %>%
  mutate(datetime = as_date(datetime)) %>%  # set proper formatting
  dplyr::filter(site_id == my_site,  # filter to desired site
                datetime < "2022-12-31") # excluding provisional data 
```

## 3.3 Visualise the target data

Let’s take a look at the targets data!

``` r
targets[100:110,]
```

    ## # A tibble: 11 × 6
    ##    project_id site_id datetime   duration variable  observation
    ##    <chr>      <chr>   <date>     <chr>    <chr>           <dbl>
    ##  1 neon4cast  OSBS    2017-06-12 P1W      richness       7     
    ##  2 neon4cast  OSBS    2017-06-26 P1W      abundance      0.0821
    ##  3 neon4cast  OSBS    2017-06-26 P1W      richness      10     
    ##  4 neon4cast  OSBS    2017-07-10 P1W      abundance      0.0446
    ##  5 neon4cast  OSBS    2017-07-10 P1W      richness       6     
    ##  6 neon4cast  OSBS    2017-07-24 P1W      abundance      0.114 
    ##  7 neon4cast  OSBS    2017-07-24 P1W      richness       8     
    ##  8 neon4cast  OSBS    2017-08-07 P1W      abundance      0.0196
    ##  9 neon4cast  OSBS    2017-08-07 P1W      richness       4     
    ## 10 neon4cast  OSBS    2017-08-21 P1W      abundance      0.0375
    ## 11 neon4cast  OSBS    2017-08-21 P1W      richness       6

It is good practice to examine the dataset before proceeding with
analysis:

``` r
targets %>% 
  as_tsibble(index = datetime, key = variable) %>%
  autoplot() +
  facet_grid(variable ~ ., scales = "free_y") + 
  theme_bw() +
  theme(legend.position = "none")
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot targets-1.png" alt="Figure: Beetle targets data at OSBS"  />
<p class="caption">
Figure: Beetle targets data at OSBS
</p>

Note that target data are available through 2022. As of the writing of
this document, some provisional 2023 data are available. The full 2023
NEON Ground Beetle dataset will be QCed during 2024 and published as a
release with a DOI in January 2025.

## 3.4 Create the training dataset

We will train our forecast models on target data from the beginning of
the dataset until our `forecast_startdate`, which we set above.

``` r
targets_train <- targets %>%
  filter(datetime < forecast_startdate) %>%
  pivot_wider(names_from = variable, values_from = observation) %>%
  as_tsibble(index = datetime)
```

## 3.5 Example forecasts: some simple models

-   Null models
    -   `fable::MEAN()`: Historical mean and standard deviation
    -   `fable::NAIVE()`: Random walk
-   Regression models with climate drivers (accessed from
    <https://open-meteo.com/>)
    -   Temperature: Daily mean temperature from CMIP6 climate model
        runs
    -   Precipitation: Daily cumulative precipitation from CMIP6 model
        runs
    -   Temperature + Precipitation

### 3.5.1 Forecast beetle abundance: null models

``` r
# specify and fit models
# Using a log(x + 1) transform on the abundance data
mod_fits <- targets_train %>% 
  tsibble::fill_gaps() %>%
  fabletools::model(
    mod_mean = fable::MEAN(log1p(abundance)),
    mod_naive = fable::NAIVE(log1p(abundance))) # random walk model, requires gapfill

# make a forecast
fc_null <- mod_fits %>%
  fabletools::forecast(h = "3 years") 
```

``` r
# visualize the forecast
fc_null %>% 
  autoplot(targets_train) +
  facet_grid(.model ~ ., scales = "free_y") +
  theme_bw()
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot null forecast-1.png" alt="Figure: NULL forecasts of ground beetle abundance at OSBS"  />
<p class="caption">
Figure: NULL forecasts of ground beetle abundance at OSBS
</p>

### 3.5.2 Forecast beetle abundance: regression models

Regression on climate model outputs allows us to make predictions about
future field seasons based on CMIP6 projections. We downloaded climate
model outputs from <https://open-meto.com> using the `RopenMeto`
package, which you can install using:
`remotes::install_github("FLARE-forecast/RopenMeteo")`.

So we do not overwhelm the open-meteo API, we have made the climate data
used in this tutorial available at: TBD CyVerse address

Climate model we’re using for this example is CMCC_CM2_VHR4

``` r
# Get climate data
path_to_clim_data <- "C:/Users/esokol/Box/00_MY_NEON/Forecasting_Beetles/future_climate_data/future_climate_2012-2050_OSBS_CMCC_CM2_VHR4.csv"

clim_long <- read_csv(path_to_clim_data)  %>%
        filter(datetime <= forecast_enddate)

# make a tsibble object
clim_long_ts <- clim_long %>%
  as_tsibble(index = datetime, 
             key = c(variable, model_id))

# make wide
clim_wide <- clim_long %>%
  select(-unit) %>%
  pivot_wider(names_from = variable, values_from = prediction)
```

``` r
# visualize climate data
clim_long_ts %>%
  ggplot(aes(datetime, prediction)) + 
  geom_line() +
  facet_grid(variable ~ ., scales = "free_y") +
  geom_vline(xintercept = lubridate::as_date(forecast_startdate),
             lty = 2) + 
  theme_bw() +
  theme(legend.position = "none")
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot climate data-1.png" alt="Figure: modeled climate data at OSBS"  />
<p class="caption">
Figure: modeled climate data at OSBS
</p>

Pick output from one model from the climate ensemble:

``` r
# subset into past and future datasets, based on forecast_startdate
clim_past <- clim_wide %>%
  filter(datetime < forecast_startdate,
         datetime > "2012-01-01")

clim_future <- clim_wide %>%
  filter(datetime >= forecast_startdate,
         datetime <= forecast_enddate)
```

Combine target and climate data to make a training dataset:

``` r
# combine target and climate data into a training dataset
targets_clim_train <- targets_train %>%
  left_join(clim_past)
```

Specify and fit simple linear regression models using `fable::TSLM()`,
examine model fit statistics.

``` r
# specify and fit model
mod_fit_candidates <- targets_clim_train %>%
  fabletools::model(
    mod_temp = fable::TSLM(log1p(abundance) ~ temperature_2m_mean),
    mod_precip = fable::TSLM(log1p(abundance) ~ precipitation_sum),
    mod_temp_precip = fable::TSLM(log1p(abundance) ~ temperature_2m_mean + precipitation_sum))

# look at fit stats
fabletools::report(mod_fit_candidates)
```

    ## # A tibble: 3 × 15
    ##   .model   r_squared adj_r_squared  sigma2 statistic p_value    df log_lik   AIC
    ##   <chr>        <dbl>         <dbl>   <dbl>     <dbl>   <dbl> <int>   <dbl> <dbl>
    ## 1 mod_temp  0.0332         0.0245  0.00416    3.78    0.0544     2    149. -610.
    ## 2 mod_pre…  0.000797      -0.00829 0.00430    0.0877  0.768      2    147. -606.
    ## 3 mod_tem…  0.0333         0.0156  0.00420    1.88    0.158      3    149. -608.
    ## # ℹ 6 more variables: AICc <dbl>, BIC <dbl>, CV <dbl>, deviance <dbl>,
    ## #   df.residual <int>, rank <int>

Plot the predicted versus observed abundance data:

``` r
# visualize model fit
# augment reformats model output into a tsibble for easier plotting
fabletools::augment(mod_fit_candidates) %>%
  ggplot(aes(x = datetime)) +
  geom_line(aes(y = abundance, lty = "Obs"), color = "dark gray") +
  geom_line(aes(y = .fitted, color = .model, lty = "Model")) +
  facet_grid(.model ~ .) +
  theme_bw()
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot tslm modeled vs observed-1.png" alt="Figure: TSLM predictions of beetle abundances at OSBS compared against observed data"  />
<p class="caption">
Figure: TSLM predictions of beetle abundances at OSBS compared against
observed data
</p>

We could use all of these models to make an ensemble forecast, but for
simplicity, we will just take the best model (lowest AICc), and use that
to create a forecast:

``` r
# focus on temperature model for now
mod_best_lm <- mod_fit_candidates %>% select(mod_temp)
report(mod_best_lm)
```

    ## Series: abundance 
    ## Model: TSLM 
    ## Transformation: log1p(abundance) 
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.07759 -0.04002 -0.01051  0.02192  0.28585 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)         -0.031335   0.051783  -0.605   0.5463  
    ## temperature_2m_mean  0.003947   0.002029   1.945   0.0544 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0645 on 110 degrees of freedom
    ## Multiple R-squared: 0.03324, Adjusted R-squared: 0.02445
    ## F-statistic: 3.782 on 1 and 110 DF, p-value: 0.054357

``` r
# make a forecast
# filter "future" climate data to target climate model
fc_best_lm <- mod_best_lm %>%
  fabletools::forecast(
    new_data = 
      clim_future %>%
      as_tsibble(index = datetime)) 
```

Visualize the forecast.

``` r
# visualize the forecast
fc_best_lm %>% 
  autoplot(targets_train) +
  facet_grid(.model ~ .) + 
  theme_bw()
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot tslm forecast-1.png" alt="Figure: TSLM forecast of beelte abundance at OSBS"  />
<p class="caption">
Figure: TSLM forecast of beelte abundance at OSBS
</p>

``` r
# format for submission to EFI
# for non-normal distributions, efi_format function draws samples to create
# n time series to provide an estimate of uncertainty
# https://projects.ecoforecast.org/neon4cast-ci/instructions.html
# I'm putting "example" in the name so the model does not register as 
# an official entry to the challenge
```

## 3.6 How to submit a forecast to the NEON Forecast Challenge

Detailed guidelines on how to submit a forecast to the NEON Forecast
Challenge can be found
[here](https://projects.ecoforecast.org/neon4cast-ci/instructions.html).
The [forecast file
format](https://projects.ecoforecast.org/neon4cast-ci/instructions.html#forecast-file-format)
requires the following columns:

-   `project_id`: use “neon4cast”

-   `model_id`: the short name of the model defined as the model_id in
    your registration. The model_id should have no spaces. model_id
    should reflect a method to forecast one or a set of target variables
    and must be unique to the neon4cast challenge.

-   `datetime`: forecast timestamp. Format `%Y-%m-%d %H:%M:%S` with UTC
    as the time zone. Forecasts submitted with a `%Y-%m-%d` format will
    be converted to a full datetime assuming UTC mid-night.

-   `reference_datetime`: The start of the forecast; this should be 0
    times steps in the future. There should only be one value of
    `reference_datetime` in the file. Format is `%Y-%m-%d %H:%M:%S` with
    UTC as the time zone. Forecasts submitted with a `%Y-%m-%d` format
    will be converted to a full datetime assuming UTC mid-night.

-   `duration`: the time-step of the forecast. Use the value of `P1D`
    for a daily forecast, `P1W` for a weekly forecast, and `PT30M` for
    30 minute forecast. This value should match the duration of the
    target variable that you are forecasting. Formatted as [ISO 8601
    duration](https://en.wikipedia.org/wiki/ISO_8601#Durations)

-   `site_id`: code for NEON site.

-   `family`: name of the probability distribution that is described by
    the parameter values in the parameter column (see list below for
    accepted distribution). An ensemble forecast as a family of
    ensemble. See note below about family

-   `parameter`: the parameters for the distribution (see note below
    about the parameter column) or the number of the ensemble members.
    For example, the parameters for a normal distribution are called
    `mu` and `sigma`.

-   `variable`: standardized variable name. It must match the variable
    name in the target file.

-   `prediction`: forecasted value for the parameter in the parameter
    column

To submit our example forecast, we can take the output from the
`fabletools::forecast()` function we used above and feed it into
`neon4cast::efi_format()` to format the output for submissions to The
Challenge. We also need to add a few additional columns.

Make sure all the required columns are included in the forecast output.

``` r
# update dataframe of model output for submission
fc_best_lm_efi <- fc_best_lm %>% 
  mutate(site_id = my_site) %>% #efi needs a NEON site ID
  neon4cast::efi_format() %>%
  mutate(
    project_id = "neon4cast",
    model_id = "bet_abund_example_tslm_temp",
    reference_datetime = forecast_startdate,
    duration = "P1W")
```

What does the content of the submission look like?

``` r
head(fc_best_lm_efi)
```

    ## # A tibble: 6 × 10
    ##   datetime   site_id parameter model_id    family variable prediction project_id
    ##   <date>     <chr>   <chr>     <chr>       <chr>  <chr>         <dbl> <chr>     
    ## 1 2022-01-01 OSBS    1         bet_abund_… ensem… abundan…   -0.0484  neon4cast 
    ## 2 2022-01-01 OSBS    2         bet_abund_… ensem… abundan…   -0.0270  neon4cast 
    ## 3 2022-01-01 OSBS    3         bet_abund_… ensem… abundan…    0.0215  neon4cast 
    ## 4 2022-01-01 OSBS    4         bet_abund_… ensem… abundan…    0.0269  neon4cast 
    ## 5 2022-01-01 OSBS    5         bet_abund_… ensem… abundan…    0.0834  neon4cast 
    ## 6 2022-01-01 OSBS    6         bet_abund_… ensem… abundan…    0.00969 neon4cast 
    ## # ℹ 2 more variables: reference_datetime <chr>, duration <chr>

``` r
# visualize the EFI-formatted submission
fc_best_lm_efi %>% 
  as_tsibble(index = datetime,
             key = c(model_id, parameter)) %>%
  ggplot(aes(datetime, prediction, color = parameter)) +
  geom_line() +
  facet_grid(model_id ~ .) +
  theme_bw()
```

<img src="NEON_forecast_challenge_beetle_tutorial_ESA2024_files/figure-markdown_github/plot tslm submission-1.png" alt="Figure: TSLM forecast submission file for OSBS, parameter indicates emsemble member"  />
<p class="caption">
Figure: TSLM forecast submission file for OSBS, parameter indicates
emsemble member
</p>

Here is how to actually submit the file to The Challenge:

``` r
# file name format is: theme_name-year-month-day-model_id.csv

# set the theme name
theme_name <- "beetles"

# set your submission date
file_date <- Sys.Date()

# make sure the model_id in the filename matches the model_id in the data
# NOTE: having the text string "example" in the file name will prevent this 
# submission from being displayed on the challenge dashboard or being included
# in statistics about challenge participation. 
efi_model_id <- "bet_abund_example_tslm_temp"

# format the file name
forecast_file <- paste0(theme_name,"-",file_date,"-",efi_model_id,".csv.gz")

# write the file to your working directory
write_csv(fc_best_lm_efi, forecast_file)

# submit the file
neon4cast::submit(forecast_file = forecast_file)
```

# 4 Evaluating your forecast

## 4.1 How your submission will be scored

The Challenge implements methods from the scoringRules R package to
calculate the Continuous Rank Probability Score (CRPS) via the
`score4cast` package, where a lower CRPS score indicates higher forecast
accuracy. CRPS uses information about the variance of the forecasts as
well as the estimated mean to calculate the score by comparing it with
the observation. There is some balance between accuracy and precision.
The forecasts will also be compared with ‘null’ models (RW and
climatology). More information can be found in the
[documentation](https://projects.ecoforecast.org/neon4cast-docs/Evaluation.html)
or the `score4cast` package from EFI organizers
[here](https://github.com/eco4cast/score4cast).

You can view past submissions to the Beetle Communities theme
[here:](https://projects.ecoforecast.org/neon4cast-dashboard/beetles.html).

You can also download the raw scores from the bucket directly, for
example:

``` r
# This example requires the `arrow` package
# install.packages("arrow")
library(arrow)

# what is your model_id?
# my_mod_id <- "bet_abund_example_tslm_temp"
my_mod_id <- "bet_example_mod_naive"
my_mod_id <- "bet_example_mod_null"

# format the URL
my_url <- paste0(
  "s3://anonymous@bio230014-bucket01/challenges/scores/parquet/project_id=neon4cast/duration=P1W/variable=abundance/model_id=",
  my_mod_id,
  "?endpoint_override=sdsc.osn.xsede.org")

# bind dataset
ds_mod_results <- arrow::open_dataset(my_url)

# get recs for dates that are scored
my_scores <- ds_mod_results %>%
  filter(!is.na(crps)) %>% 
  collect()

head(my_scores)
```

## 4.2 How to score your own forecast

For immediate feedback, we can use the targets data from 2022 to score
our forecast for the 2022 field season at OSBS.

``` r
# filter to 2022 because that is the latest release year
# 2023 is provisional and most sites do not yet have data reported
targets_2022 <- targets %>% 
  dplyr::filter(
    datetime >= "2022-01-01", 
    datetime < "2023-01-01",
    variable == "abundance",
    observation > 0)

# list of target site dates for filtering mod predictions
target_site_dates_2022 <- targets_2022 %>%
  select(site_id, datetime) %>% distinct()

# filter model forecast data to dates where we have observations
mod_results_to_score <- fc_best_lm_efi %>%
  left_join(target_site_dates_2022,.) %>%
  dplyr::filter(!is.na(parameter))

# score the forecasts
mod_scores <- score(
  forecast = mod_results_to_score,
  target = targets_2022) 

head(mod_scores)
```

    ## # A tibble: 6 × 17
    ##   model_id     reference_datetime site_id datetime   family variable observation
    ##   <chr>        <chr>              <chr>   <date>     <chr>  <chr>          <dbl>
    ## 1 bet_abund_e… 2022-01-01         OSBS    2022-04-04 sample abundan…      0.102 
    ## 2 bet_abund_e… 2022-01-01         OSBS    2022-04-18 sample abundan…      0.188 
    ## 3 bet_abund_e… 2022-01-01         OSBS    2022-04-25 sample abundan…      0.0877
    ## 4 bet_abund_e… 2022-01-01         OSBS    2022-05-02 sample abundan…      0.0857
    ## 5 bet_abund_e… 2022-01-01         OSBS    2022-05-16 sample abundan…      0.0786
    ## 6 bet_abund_e… 2022-01-01         OSBS    2022-05-30 sample abundan…      0.133 
    ## # ℹ 10 more variables: crps <dbl>, logs <dbl>, mean <dbl>, median <dbl>,
    ## #   sd <dbl>, quantile97.5 <dbl>, quantile02.5 <dbl>, quantile90 <dbl>,
    ## #   quantile10 <dbl>, horizon <drtn>

``` r
# compare scores from best_lm, mod_mean, and mod_naive
```
