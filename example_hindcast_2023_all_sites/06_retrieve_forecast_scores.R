
# libraries ----
library(tidyverse)
library(arrow)

# helper functions for format for EFI submission
# library(neon4cast) # remotes::install_github("eco4cast/neon4cast", dep=TRUE)


# my model names
list_of_mod_names <- c("bet_example_mod_naive", "bet_example_mod_null",
  "bet_lm_precip_clim_ensemble", "bet_lm_temp_clim_ensemble", 
  "bet_lm_temp_precip_clim_ensemble")

# example forecast results url:
# https://radiantearth.github.io/stac-browser/#/external/raw.githubusercontent.com/eco4cast/neon4cast-ci/main/catalog/scores/models/model_items/bet_lm_temp2m_mod_temp.json?.language=en&.asset=asset-3



# loop go get crps scores
all_results <- data.frame()
for(my_mod_id in list_of_mod_names){
  try({
    # structure url for mod scores dataset bucket
    my_url <- paste0("s3://anonymous@bio230014-bucket01/challenges/scores/parquet/project_id=neon4cast/duration=P1W/variable=abundance/model_id=",
                     my_mod_id,
                     "?endpoint_override=sdsc.osn.xsede.org")
    # bind dataset
    ds_mod_results <- arrow::open_dataset(my_url)
    
    # get recs for dates that are scored
    all_results <- ds_mod_results %>%
      filter(!is.na(crps)) %>% 
      collect() %>%
      mutate(model_id = my_mod_id) %>%
      bind_rows(all_results, .)
  })
  message("loaded: ",my_mod_id)
}


# write out results
write_csv(all_results, "forecast_scores.csv")


# plot results
plot_CRPS_allsites <- all_results %>%
  ggplot(aes(datetime, log10(crps), color = model_id)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  ggtitle("NEON Beetle Abundance 2023 Hindcast Scores")

# # loop through plots of CRPS scores for each site
# site_list <- all_results$site_id %>% unique()
# for(i_site in site_list){
#   p <- all_results %>%
#     dplyr::filter(site_id == i_site) %>%
#     ggplot(aes(datetime, log10(crps), color = model_id)) +
#     geom_point() +
#     geom_smooth() +
#     theme_bw() + 
#     ggplot2::ggtitle(i_site)
#   
#   print(p)
#   message(i_site, " plotted... ")
#   scan()
# }

# MOAB - lm mods with temp do well in off season, probably underestimate during peak season

