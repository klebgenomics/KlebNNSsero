##------------------------------------------------------------------------------ 
#
#   Model Run File
#
##------------------------------------------------------------------------------ 

# Load required packages
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(stringr)
library(tidyverse)


## SECTION ONE: Data preparation -----------------------------------------------

# Part A:


all_files <- list.files(path = "data_paper/", recursive = TRUE)
all_files <- all_files[all_files != "site_info.csv"]
print(all_files)

pattern <- "^(Full|Carba|ESBL|Fatal)_(ALL|min10)_Neonatal_shareable_\\d{8}_(28|365)_filterN10_(K|O|OlocusType)_site_counts_all\\.csv$"

# 2) Identify files that do NOT match
bad_files <- all_files[!grepl(pattern, all_files)]

# 3) If the vector of bad_files is non-empty, flag them. Otherwise, confirm all okay.
if (length(bad_files) > 0) {
  warning("The following files do not match the required format:\n", 
          paste(bad_files, collapse = "\n"))
} else {
  message("All files match the required format.")
}


## PLACEHOLDER TO TEST
all_files <- all_files[c(1,3)]

for(i in 1:length(all_files)){
  
  ## LOOP ONE: Data preparation ------------------------------------------------
  
  # File name for this loop
  filename <- all_files[i]
  
  # Extract metadata from the filename using regex patterns
  subset     <- str_extract(filename, "^(Full|Carba|ESBL|Fatal)")
  purpose    <- ifelse(str_detect(filename, "shareable"), "paper", "BMGF")
  date_stamp <- str_extract(filename, "\\d{8}")  # YYYYMMDD
  days       <- str_extract(filename, "_(28|365)_") %>% str_replace_all("_", "")
  what_data  <- str_extract(filename, "ALL|min10")
  antigen    <- str_extract(filename, "_(K|O|OlocusType)_") %>% str_replace_all("_", "")
  
  # Check if "LOO" is in the filename. If so, extract the study name.
  if (str_detect(filename, "LOO")) {
    # The pattern looks for the text immediately following "LOO_" until the next period.
    study <- str_extract(filename, "(?<=LOO_)[^.]+")
  } else {
    study <- NA
  }
  
  # Read the data
  site_prevalence <- read.csv(paste("data_",purpose,"/",filename, sep=""))

  # Load datasets
  site_data <- read.csv(paste("data_",purpose,"/","site_info.csv", sep=""))
  
  # Ensure site_info has unique site-region mapping only
  site_info <- site_prevalence %>%
    select(Site, Region) %>%
    distinct()
  
  # Ensure site_prevalence doesn't contain multiple versions of Region before join
  site_prevalence <- site_prevalence %>%
    select(-Region) %>%  # Remove existing Region to prevent duplication
    left_join(site_info, by = "Site")  # Now join to get Region from site_info
  
  for(type in c("adj","raw")){
    
    ## LOOP TWO(a): Filter on type and expand to include 0 event sites ---------
    
    # Create datasets based on type selection
    site_prevalence_type <- site_prevalence %>%
      select(locus, Site, Region, contains(type)) %>%
      rename(subgroup = Region, event = paste0(type, "_count"), n = paste0(type, "_sum"))
    
    # First, identify each unique subgroup-Site pair with its n:
    Site_subgroup_n <- site_prevalence_type %>%
      distinct(subgroup, Site, n)
    
    # Create all possible combinations of locus, subgroup, Site
    all_combinations <- expand.grid(
      locus = unique(site_prevalence_type$locus),
      subgroup = unique(site_prevalence_type$subgroup),
      Site = unique(site_prevalence_type$Site),
      KEEP.OUT.ATTRS = FALSE  # just to keep the output cleaner
    ) %>%
      as_tibble()
    
    # Join the subgroup-Site n-values into this expanded set
    expanded_df <- all_combinations %>%
      left_join(Site_subgroup_n, by = c("subgroup", "Site"))
    
    # Bring in the original data’s event (if present)
    # Any missing event gets replaced with 0
    expanded_df <- expanded_df %>%
      left_join(
        site_prevalence_type %>% select(subgroup, Site, locus, event),
        by = c("subgroup", "Site", "locus")
      ) %>%
      mutate(event = if_else(is.na(event), 0, event))
    
    # expanded_df now has every (subgroup, Site, locus) combo with:
    #   - the correct n from Site_subgroup_n
    #   - event set to 0 if it wasn’t found in data_csv
    # This ensures that there is an explicit row for all combinations.
    
    
    # n is NA for site combinations that do not exist. Filter NAs so that these 
    # combinations are removed from the final dataset
    expanded_df <- expanded_df %>%
      tidyr::drop_na(n)
    
    data_csv <- expanded_df
    
    
    ## LOOP TWO(b): Run model --------------------------------------------------
    
    
    fit_with_interaction <- brm(
      formula = bf(
        event | trials(n) ~ 0 + locus + subgroup + 
          (1 | Site) + 
          (1 | locus:subgroup)
      ),
      data = data_csv,
      family = binomial(link = "logit"),
      prior = c(
        prior(student_t(3, 0, 5), class = "sd", group = "Site"),
        prior(student_t(3, 0, 5), class = "sd", group = "locus:subgroup")
      ),
      iter = 20000,
      chains = 4,
      cores = 4,
      control = list(adapt_delta = 0.999999, max_treedepth = 55)
    )
    
    
    # Check model:
    
    max(rhat(fit_with_interaction))
    min(neff_ratio(fit_with_interaction))
    pp_check(fit_with_interaction,type = "stat_2d")
    
    # Save model:
    
    if (!is.na(study)) {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model_LOO_", study, ".rds")
    } else {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model.rds")
    }
    
    model_filepath <- paste0("models_", purpose, "/", model_filename)
    
    saveRDS(fit_with_interaction, file = model_filepath)
  }
  
}
