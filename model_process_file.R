##------------------------------------------------------------------------------ 
#
#   Model Process File
#
##------------------------------------------------------------------------------ 


# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)
library(brms)


all_files <- list.files(path = "data_core/", recursive = TRUE) # or 'data_LOO/'
all_files <- all_files[all_files != "site_info.csv"]
print(all_files)

pattern <- "^(Full|Carba|ESBL|Fatal)_(ALL|min10)_Neonatal_shareable_\\d{8}_(28|365)_filterN10_(K|O|OlocusType)_site_counts_all\\.csv$"

# Identify files that do NOT match
bad_files <- all_files[!grepl(pattern, all_files)]

# If the vector of bad_files is non-empty, flag them. Otherwise, confirm all okay.
if (length(bad_files) > 0) {
  warning("The following files do not match the required format:\n", 
          paste(bad_files, collapse = "\n"))
} else {
  message("All files match the required format.")
}


## PLACEHOLDER TO TEST
all_files <- all_files[13:14]

# Loop over each file in all_files
for(i in seq_along(all_files)) {
  # File name for this iteration
  filename <- all_files[i]
  
  # Extract metadata from filename
  subset      <- str_extract(filename, "^(Full|Carba|ESBL|Fatal)")
  purpose     <- ifelse(str_detect(filename, "LOO"), "LOO", "core")
  date_stamp  <- str_extract(filename, "\\d{8}")         # Extracts YYYYMMDD
  days        <- str_extract(filename, "_(28|365)_") %>% str_replace_all("_", "")
  what_data   <- str_extract(filename, "ALL|min10")
  antigen     <- str_extract(filename, "_(K|O|OlocusType)_") %>% str_replace_all("_", "")
  
  # Check if "LOO" is in the filename. If so, extract the study name.
  if (str_detect(filename, "LOO")) {
    # The pattern looks for the text immediately following "LOO_" until the next period.
    study <- str_extract(filename, "(?<=LOO_)[^.]+")
  } else {
    study <- NA
  }
  
  # Loop over the types: "adj" and "raw"
  for(type in c("adj", "raw")) {
    
    ## Load the saved model ---------------------------
    # Construct the model file name
    
    if (!is.na(study)) {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model_LOO_", study, ".rds")
    } else {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model.rds")
    }
    model_filepath <- paste0("models_", purpose, "/", model_filename)
    
    
    # Load the model
    fit_with_interaction <- readRDS(model_filepath)
    
    ## Load the data and process to obtain data_csv -----------
    # Read the prevalence data file and the site information file
    site_prevalence <- read.csv(paste0("data_", purpose, "/", filename))
    site_data     <- read.csv(paste0("data_", purpose, "/", "site_info.csv"))
    
    # Create unique site-to-region mapping from the prevalence data
    site_info <- site_prevalence %>%
      select(Site, Region) %>%
      distinct()
    
    # Remove Region from the main dataset to prevent duplication and re-join with site_info
    site_prevalence <- site_prevalence %>%
      select(-Region) %>%
      left_join(site_info, by = "Site")
    
    # Process data based on the type ("adj" or "raw")
    site_prevalence_type <- site_prevalence %>%
      select(locus, Site, Region, contains(type)) %>%
      rename(subgroup = Region,
             event = paste0(type, "_count"),
             n = paste0(type, "_sum"))
    
    # Identify unique Site-subgroup pairs with their n-values
    Site_subgroup_n <- site_prevalence_type %>%
      distinct(subgroup, Site, n)
    
    # Create all possible combinations of locus, subgroup, and Site
    all_combinations <- expand.grid(
      locus    = unique(site_prevalence_type$locus),
      subgroup = unique(site_prevalence_type$subgroup),
      Site     = unique(site_prevalence_type$Site),
      KEEP.OUT.ATTRS = FALSE
    ) %>% as_tibble()
    
    # Merge the n-values into the expanded grid
    expanded_df <- all_combinations %>%
      left_join(Site_subgroup_n, by = c("subgroup", "Site"))
    
    # Merge in the original event data; set missing events to 0
    expanded_df <- expanded_df %>%
      left_join(
        site_prevalence_type %>% select(subgroup, Site, locus, event),
        by = c("subgroup", "Site", "locus")
      ) %>%
      mutate(event = if_else(is.na(event), 0, event))
    
    # Remove any rows where n is NA (these correspond to site combinations that do not exist)
    expanded_df <- expanded_df %>%
      tidyr::drop_na(n)
    
    # Final dataset ready for analysis
    data_csv <- expanded_df
    
    # You now have both the loaded model (fit_with_interaction) and data_csv ready for further analysis
    message("Loaded model and processed data for: ", model_filename)
    
    
    
    ## LOOP TWO(b): Global Processing ------------------------------------------
    
    post_epred <- posterior_epred(
      fit_with_interaction,
      newdata    = data_csv, 
      re_formula = NULL,     # includes all REs
      allow_new_levels = TRUE,
      summary    = FALSE
    )
    
    n_draws <- nrow(post_epred)
    
    # Identify all unique loci
    unique_loci <- unique(data_csv$locus)
    
    # We'll store a row for each (draw, locus), with the posterior prevalence.
    posterior_locus_list <- list()
    
    for (loc in unique_loci) {
      # Rows of data_csv for this locus
      idx <- data_csv$locus == loc
      
      # Sum of 'n' for this locus
      total_n_loc <- sum(data_csv$n[idx], na.rm = TRUE)
      
      # For each MCMC draw, compute sum of predicted events / sum of n
      prevalence_draws <- numeric(n_draws)
      for (d in seq_len(n_draws)) {
        total_events_d <- sum(post_epred[d, idx])
        prevalence_draws[d] <- total_events_d / total_n_loc
      }
      
      # Build a mini data frame with (draw, locus, prevalence_draw)
      posterior_locus_list[[loc]] <- data.frame(
        locus = loc,
        draw  = seq_len(n_draws),
        prevalence_draw = prevalence_draws
      )
    }
    
    # Combine all loci into one data frame
    posterior_locus_df <- do.call(rbind, posterior_locus_list)
    rownames(posterior_locus_df) <- NULL
    colnames(posterior_locus_df) <- c("locus", "draw_id", "prob")
    posterior_locus_df$subgroup <- "Global"
    
    
    global_estimates <- posterior_locus_df %>%
      group_by(locus) %>%
      summarise(
        mean = mean(prob),
        median = median(prob),
        lower95 = quantile(prob, 0.025, na.rm = TRUE),
        upper95 = quantile(prob, 0.975, na.rm = TRUE)
      )
    
    global_estimates$subgroup <- "Global"
    colnames(global_estimates) <- c("locus", "mean", "median", "lower", "upper", "subgroup")

    
    
    # Create the output folder based on the purpose (e.g., "outputs_paper")
    output_folder <- paste0("outputs_", purpose, "/")
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Construct file names conditionally based on the presence of LOO (i.e., study not NA)
    if (!is.na(study)) {
      posterior_filename <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_LOO_", study, "_posterior_global.csv.gz"
      )
      estimates_filename <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_LOO_", study, "_global_estimates.csv"
      )
    } else {
      posterior_filename <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_posterior_global.csv.gz"
      )
      estimates_filename <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_global_estimates.csv"
      )
    }
    
    # Adjust filenames if they start with "OlocusType"
    posterior_filename <- sub("^OlocusType", "O", posterior_filename)
    estimates_filename <- sub("^OlocusType", "O", estimates_filename)
    
    # Save the posterior_locus_df as a gzipped CSV file.
    write.csv(posterior_locus_df,
              gzfile(file.path(output_folder, posterior_filename)),
              row.names = FALSE)
    
    # Save the global_estimates as a CSV file.
    write.csv(global_estimates,
              file.path(output_folder, estimates_filename),
              row.names = FALSE)
    
    
    
    ## LOOP TWO(c): Regional Processing ----------------------------------------
    
    post_epred <- posterior_epred(
      fit_with_interaction,
      newdata    = data_csv,
      re_formula = NULL,   # include all random effects
      allow_new_levels = TRUE,
      summary    = FALSE
    )
    
    n_draws <- nrow(post_epred)  # total MCMC draws
    
    # 4) Build Posterior Distributions for Each (locus, subgroup) ------------------
    # For each draw, we sum predicted events in rows for that (locus, subgroup),
    # then divide by sum of n in those rows to get the posterior prevalence.
    
    unique_locus_subgroups <- data_csv %>%
      distinct(locus, subgroup)
    
    posterior_list <- list()
    
    for (k in seq_len(nrow(unique_locus_subgroups))) {
      
      loc  <- unique_locus_subgroups$locus[k]
      sgrp <- unique_locus_subgroups$subgroup[k]
      
      # Identify rows of data_csv for this (locus, subgroup)
      idx <- (data_csv$locus == loc & data_csv$subgroup == sgrp)
      
      # Sum of n for those rows
      total_n_locus_subgrp <- sum(data_csv$n[idx], na.rm = TRUE)
      
      # Posterior prevalence draws
      prevalence_draws <- numeric(n_draws)
      
      for (d in seq_len(n_draws)) {
        total_events_d <- sum(post_epred[d, idx])
        prevalence_draws[d] <- total_events_d / total_n_locus_subgrp
      }
      
      # Build mini data frame
      posterior_list[[k]] <- data.frame(
        locus           = loc,
        subgroup        = sgrp,
        draw            = seq_len(n_draws),
        prevalence_draw = prevalence_draws
      )
    }
    
    # Combine all into one data frame
    posterior_locus_subgroup_df <- do.call(rbind, posterior_list)
    head(posterior_locus_subgroup_df)
    colnames(posterior_locus_subgroup_df) <- c("locus", "subgroup", "draw_id", "prob")
    
    
    subgroup_estimates <- posterior_locus_subgroup_df %>%
      group_by(locus, subgroup) %>%
      summarize(
        mean   = mean(prob),
        median = median(prob),
        lower = quantile(prob, 0.025, na.rm = TRUE),
        upper = quantile(prob, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    # Construct file names conditionally based on the presence of LOO (i.e., study not NA)
    if (!is.na(study)) {
      posterior_filename_subgroup <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_LOO_", study, "_posterior_subgroup.csv.gz"
      )
      estimates_filename_subgroup <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_LOO_", study, "_subgroup_estimates.csv"
      )
    } else {
      posterior_filename_subgroup <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_posterior_subgroup.csv.gz"
      )
      estimates_filename_subgroup <- paste0(
        antigen, "_", subset, "_", what_data, "_", type, "_", days, "_subgroup_estimates.csv"
      )
    }
    
    # Replace "OlocusType" with "O" if present at the beginning of the filenames
    posterior_filename_subgroup <- sub("^OlocusType", "O", posterior_filename_subgroup)
    estimates_filename_subgroup <- sub("^OlocusType", "O", estimates_filename_subgroup)
    
    # Save the posterior_locus_subgroup_df as a gzipped CSV file
    write.csv(posterior_locus_subgroup_df,
              gzfile(file.path(output_folder, posterior_filename_subgroup)),
              row.names = FALSE)
    
    # Save the subgroup_estimates as a CSV file
    write.csv(subgroup_estimates,
              file.path(output_folder, estimates_filename_subgroup),
              row.names = FALSE)
    
    
  } # End loop for type
} # End loop for files




