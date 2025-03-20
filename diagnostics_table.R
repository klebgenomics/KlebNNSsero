library(dplyr)
library(stringr)
library(tidyr)
library(brms)

all_files <- list.files(path = "data_paper/", recursive = TRUE)
all_files <- all_files[all_files != "site_info.csv"]

pattern <- "^(Full|Carba|ESBL|Fatal)_(ALL|min10)_Neonatal_shareable_\\d{8}_(28|365)_filterN10_(K|O|OlocusType)_site_counts_all\\.csv$"

bad_files <- all_files[!grepl(pattern, all_files)]
if (length(bad_files) > 0) {
  warning("The following files do not match the required format:\n", paste(bad_files, collapse = "\n"))
} else {
  message("All files match the required format.")
}

# Placeholder to test a subset of files
#all_files <- all_files[c(1,3)]

# Initialize results table
results_table <- data.frame(
  model = character(),
  max_rhat = numeric(),
  min_neff_ratio = numeric(),
  divergent_transitions = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(all_files)) {
  filename <- all_files[i]
  subset <- str_extract(filename, "^(Full|Carba|ESBL|Fatal)")
  purpose <- "paper"
  date_stamp <- str_extract(filename, "\\d{8}")
  days <- str_extract(filename, "_(28|365)_") %>% str_replace_all("_", "")
  what_data <- str_extract(filename, "ALL|min10")
  antigen <- str_extract(filename, "_(K|O|OlocusType)_") %>% str_replace_all("_", "")
  
  if (str_detect(filename, "LOO")) {
    study <- str_extract(filename, "(?<=LOO_)[^.]+")
  } else {
    study <- NA
  }
  
  for (type in c("adj", "raw")) {
    if (!is.na(study)) {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model_LOO_", study, ".rds")
    } else {
      model_filename <- paste0(antigen, "_", subset, "_", what_data, "_", type, "_", days, "days_model.rds")
    }
    model_filepath <- paste0("models_", purpose, "/", model_filename)
    
    fit_with_interaction <- readRDS(model_filepath)
    
    model_summary <- summary(fit_with_interaction)
    diagnostics <- rhat(fit_with_interaction)
    
    max_rhat <- max(diagnostics, na.rm = TRUE)
    min_neff_ratio <- min(neff_ratio(fit_with_interaction), na.rm = TRUE)
    divergent_transitions <- sum(fit_with_interaction$fit@sim$sampler_diagnostics$divergent__)
    
    results_table <- rbind(results_table, data.frame(
      model = model_filename,
      max_rhat = max_rhat,
      min_neff_ratio = min_neff_ratio,
      divergent_transitions = divergent_transitions,
      stringsAsFactors = FALSE
    ))
  }
}

# Display the results table
print(results_table)


write.csv(results_table, "results_table.csv")

library(dplyr)
library(stringr)

results_table_new <- results_table %>%
  mutate(
    subset     = str_extract(model, "(Full|Carba|ESBL|Fatal)"),
    type       = if_else(str_detect(model, "raw"), "raw", "adj"),
    date_stamp = str_extract(model, "\\d{8}"),                   # Will be NA if no date present
    days       = str_extract(model, "(28|365)"),                  # Extracts 28 or 365
    what_data  = str_extract(model, "ALL|min10"),
    antigen    = str_extract(model, "^(K|O)")                     # Extracts a leading K or O
  ) %>%
  select(subset, type, days, what_data, antigen,
         max_rhat, min_neff_ratio)

print(results_table_new)
write.csv(results_table_new, "results_table.csv")

 