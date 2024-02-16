remotes::install_github("FLARE-forecast/RopenMeteo")

#SCBI or SERC
# SCBI 38.892925, -78.139494

my_site_id <- "OSBS"

# get NEON site list
neon_site_info <- read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv")

my_neon_site <- neon_site_info %>% dplyr::filter(field_site_id == my_site_id)

# Quinn's R package to get met and climate data

clim_mod_names <- c(
  "CMCC_CM2_VHR4",
  "FGOALS_f3_H",
  "HiRAM_SIT_HR",
  "MRI_AGCM3_2_S",
  "EC_Earth3P_HR",
  "MPI_ESM1_2_XR",
  "NICAM16_8S")

library(RopenMeteo)

clim_dat_all <- data.frame()
for(i_mod in clim_mod_names){
  clim_dat <- data.frame()
  try({
    clim_dat <- RopenMeteo::get_climate_projections(
      latitude = my_neon_site$field_latitude,
      longitude = my_neon_site$field_longitude,
      site_id = my_site_id,
      start_date = "2012-01-01",end_date = "2050-01-01",
      model = i_mod,
      variables = c("precipitation_sum","temperature_2m_mean")) #"temperature_2m_mean", "precipitation_sum"
  })
  
  clim_dat_all <- bind_rows(
    clim_dat_all, clim_dat)
  
  message(i_mod)
}

clim_dat_to_save <- clim_dat_all

clim_dat_to_save$model_id %>% unique()
write_csv(clim_dat_to_save, "clim_dat_OSBS_test.csv")
