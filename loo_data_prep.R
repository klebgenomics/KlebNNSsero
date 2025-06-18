##------------------------------------------------------------------------------ 
#
#   LOO data prep
#
##------------------------------------------------------------------------------ 

all_files <- list.files(path = "data_core/", recursive = TRUE)
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


all_files <- all_files[8]

filename <- all_files

# Extract metadata from filename
subset      <- str_extract(filename, "^(Full|Carba|ESBL|Fatal)")
purpose     <- ifelse(str_detect(filename, "shareable"), "paper", "BMGF")
date_stamp  <- str_extract(filename, "\\d{8}")         # Extracts YYYYMMDD
days        <- str_extract(filename, "_(28|365)_") %>% str_replace_all("_", "")
what_data   <- str_extract(filename, "ALL|min10")
antigen     <- str_extract(filename, "_(K|O|OlocusType)_") %>% str_replace_all("_", "")

# Read the data
site_prevalence <- read.csv(paste("data_core","/",filename, sep=""))
# Load datasets
site_data <- read.csv(paste("data_core","/","site_info_250529.csv", sep=""))
#site_data <- read.csv(paste("data_core","/","site_info.csv", sep=""))

site_prevalence <- site_prevalence %>%
  left_join(site_data %>% select(Site, Study), by = "Site")

# Get unique studies from the dataset
unique_studies <- unique(site_prevalence$Study)
print(unique_studies)

# For Subsetting
unique_studies <- c("NeoOBS","AIIMS", "MBIRA")


for(study in unique_studies) {
  # Remove rows corresponding to the current study
  loo_data <- site_prevalence %>% filter(Study != study)
  
  # Construct the new filename by appending "_LOO_[study]" before the .csv extension
  new_filename <- sub("\\.csv$", paste0("_LOO_", study, ".csv"), filename)
  
  # Save the new dataset
  write.csv(loo_data, file = paste0("data_LOO", "/", new_filename), row.names = FALSE)
  
  message("Saved file: ", new_filename)
}

