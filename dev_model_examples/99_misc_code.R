# beetle targets are here
url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1W/beetles-targets.csv.gz"

bet_targets <- read_csv(url)

bet_targets$site_id %>% unique()
