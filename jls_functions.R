
dependencies <- function () {

  c("neonstore", "tidyverse", "ISOweek", "dplyr", "runjags", "stringi", "tsibble", "lubridate", "fable", "neon4cast")


}

load_dependencies <- function (packages = dependencies()) {

  ndepends <- length(packages)
  present  <- installed.packages()[ , "Package"]
  needed   <- packages[!(packages %in% present)] 
  nneeded  <- length(needed)


  if ("neon4cast" %in% needed) {

    if ("remotes" %in% needed) {

      install.packages("remotes")

    }

    remotes::install_github("eco4cast/neon4cast")    

  }

  present  <- installed.packages()[ , "Package"]
  needed   <- packages[!(packages %in% present)] 
  nneeded  <- length(needed)


  if (nneeded > 0) {

    install.packages(needed)

  }

  for (dep in 1:ndepends) {

    eval(bquote(library(.(packages[dep]))))

  }

}




clean_names <- function (x) {
  s <- stringi::stri_split_regex(x, "/", simplify = TRUE)[,1]
  s <- stringi::stri_extract_all_words(s, simplify = TRUE)
  if (dim(s)[2] > 1)
    stringi::stri_trim(paste(s[, 1], s[, 2]))
  else stringi::stri_trim(s[, 1])
}

resolve_taxonomy <- function(sorting, para, expert){

  taxonomy <-
    left_join(sorting,
              select(para, subsampleID, individualID, scientificName, taxonRank, taxonID, morphospeciesID),
              by = "subsampleID")  %>%
    ## why are there so many other shared columns (siteID, collectDate, etc?  and why don't they match!?)
    ## we use `select` to avoid these
    left_join(
      select(expert, -uid, -namedLocation, -domainID, -siteID, -collectDate, -plotID, -setDate, -collectDate),
      by = "individualID") %>%
    distinct() %>%
     ## Prefer the para table cols over the sorting table cols only for sampleType=="other carabid"
    mutate(taxonRank.x = ifelse(is.na(taxonRank.y) | sampleType != "other carabid", taxonRank.x, taxonRank.y),
           scientificName.x = ifelse(is.na(scientificName.y) | sampleType != "other carabid", scientificName.x, scientificName.y),
           taxonID.x = ifelse(is.na(taxonID.y) | sampleType != "other carabid", taxonID.x, taxonID.y),
           morphospeciesID.x =  ifelse(is.na(morphospeciesID.y) | sampleType != "other carabid", morphospeciesID.x, morphospeciesID.y)) %>%
      ## Prefer expert values where available
    mutate(taxonRank = ifelse(is.na(taxonRank), taxonRank.x, taxonRank),
           scientificName = ifelse(is.na(scientificName), scientificName.x, scientificName),
           taxonID = ifelse(is.na(taxonID), taxonID.x, taxonID),
           morphospeciesID =  ifelse(is.na(morphospeciesID), morphospeciesID.x, morphospeciesID),
           nativeStatusCode = ifelse(is.na(nativeStatusCode.y), nativeStatusCode.x, nativeStatusCode.y),
           sampleCondition = ifelse(is.na(sampleCondition.y), sampleCondition.x, sampleCondition.y)
           ) %>%
    select(-ends_with(".x"), -ends_with(".y")) # %>%
#    select(-individualCount)
     ## WARNING: if the subsample is split into separate taxa by experts, we do not know
     ## how many of the total count should go to each taxon in the the split
     ## since only part of that subsample have been pinned.
     ## We should flag these cases in some manner.


  #### Should we add a "species" column, using morphospecies or the best available?
  ## Use morphospecies if available for higher-rank-only classifications,
  ## Otherwise, binomialize the scientific name:
  taxonomy <- taxonomy %>%
    mutate(morphospecies =
             ifelse(taxonRank %in% c("subgenus", "genus", "family", "order") & !is.na(morphospeciesID),
                    morphospeciesID,
                    clean_names(scientificName)
             )
    )

  ## Beetles must be identified as carabids by both sorting table and the taxonomists (~3 non-Carabidae slip through in sorting)
  beetles <- taxonomy %>%
    filter(grepl("carabid", sampleType)) %>%
    filter(family == "Carabidae" | is.na(family))

  beetles
}


named_null_list <- function (element_names = NULL) {

  return_if_null(element_names)

  nelements  <- length(element_names)
  out        <- vector("list", nelements)
  names(out) <- element_names

  out

}

return_if_null <- function (x, value = NULL) {
    if (is.null(x)) {
        do.call(what = "return", args = list(value), envir = sys.frame(-1))
    }
}


runjags_inits <- function (inits) {

  rngs  <- c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")

  function (data = NULL) {

    function(chain = chain) {

      model_specific_inits <- named_null_list(element_names = names(inits))
      for (i in 1:length(model_specific_inits)) {

        model_specific_inits[[i]] <- eval(parse(text = inits[[i]]))

      }

      c(.RNG.name = sample(rngs, 1),
        .RNG.seed = sample(1:1e+06, 1), 
        model_specific_inits)                  
      
    }

  }


}

runjags_model <- function (model_file) {

  scan(file  = model_file,
       what  = "character",
       quiet = TRUE)
  
}


fit_runjags_null_abundance <- function (past, 
                                        future, 
                                        control_runjags = runjags_controls( )) {



  jags_model <- runjags_model(model_file = "jags_null_abundance.txt")

  abundance_rate          <- past$abundance
  trapnights              <- past$trapnights
  abundance               <- round(abundance_rate * trapnights)
  N                       <- length(abundance)
  log_mean_abundance_rate <- log(mean(abundance_rate, na.rm = TRUE))  

  abundance_rate_future   <- future$abundance
  trapnights_future       <- future$trapnights
  abundance_future        <- abundance_rate_future * trapnights_future

  N_future                <- length(abundance_future)

  monitor                 <- c("mu", paste0("abundance_future[", 1:N_future, "]"))

  data <- list(abundance               = abundance, 
               N                       = N,
               trapnights              = trapnights,
               log_mean_abundance_rate = log_mean_abundance_rate,
               N_future                = N_future,
               trapnights_future       = trapnights_future)

  init       <- runjags_inits(inits = list(mu = rnorm(1, mean = log_mean_abundance_rate, sd = 0.1)))

  runjags.options(silent.jags    = control_runjags$silent_jags, 
                  silent.runjags = control_runjags$silent_jags)

  model_fit <- run.jags(model     = jags_model, 
                        monitor   = monitor, 
                        inits     = init(data), 
                        data      = data, 
                        n.chains  = control_runjags$nchains, 
                        adapt     = control_runjags$adapt, 
                        burnin    = control_runjags$burnin, 
                        sample    = control_runjags$sample, 
                        thin      = control_runjags$thin, 
                        modules   = control_runjags$modules, 
                        method    = control_runjags$method, 
                        factories = control_runjags$factories, 
                        mutate    = control_runjags$mutate, 
                        summarise = TRUE, 
                        plots     = FALSE)


}


fit_runjags_null_richness <- function (past, 
                                       future, 
                                       control_runjags = runjags_controls( )) {

  jags_model <- runjags_model(model_file = "jags_null_richness.txt")

  richness_rate          <- past$richness
  trapnights             <- past$trapnights
  richness               <- round(richness_rate * trapnights)
  N                      <- length(richness)
  log_mean_richness_rate <- log(mean(richness_rate, na.rm = TRUE))  

  richness_rate_future   <- future$richness
  trapnights_future      <- future$trapnights
  richness_future        <- richness_rate_future * trapnights_future

  N_future                <- length(richness_future)

  monitor                 <- c("mu", paste0("richness_future[", 1:N_future, "]"))

  data <- list(richness               = richness, 
               N                      = N,
               trapnights             = trapnights,
               log_mean_richness_rate = log_mean_richness_rate,
               N_future               = N_future,
               trapnights_future      = trapnights_future)

  init       <- runjags_inits(inits = list(mu = rnorm(1, mean = log_mean_richness_rate, sd = 0.1)))

  runjags.options(silent.jags    = control_runjags$silent_jags, 
                  silent.runjags = control_runjags$silent_jags)

  model_fit <- run.jags(model     = jags_model, 
                        monitor   = monitor, 
                        inits     = init(data), 
                        data      = data, 
                        n.chains  = control_runjags$nchains, 
                        adapt     = control_runjags$adapt, 
                        burnin    = control_runjags$burnin, 
                        sample    = control_runjags$sample, 
                        thin      = control_runjags$thin, 
                        modules   = control_runjags$modules, 
                        method    = control_runjags$method, 
                        factories = control_runjags$factories, 
                        mutate    = control_runjags$mutate, 
                        summarise = TRUE, 
                        plots     = FALSE)


}



runjags_controls <- function (nchains     = 3, 
                              adapt       = 1e4, 
                              burnin      = 1e4, 
                              sample      = 1e4, 
                              thin        = 1, 
                              modules     = "glm", 
                              method      = "interruptible", 
                              factories   = "", 
                              mutate      = NA, 
                              silent_jags = FALSE){

  list(nchains     = nchains, 
       adapt       = adapt, 
       burnin      = burnin, 
       sample      = sample,
       thin        = thin, 
       modules     = modules, 
       method      = method, 
       factories   = factories,
       mutate      = mutate, 
       silent_jags = silent_jags)

}

