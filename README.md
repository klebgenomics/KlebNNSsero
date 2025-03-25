# README

## Figures and tables from the paper

Directories `figures` and `tables` contain all figures and tables from the paper.

The code to generate these are in the following R markdown files, which draw on functions in `seroepi_functions.R`:

`DataAnalysis.Rmd` - R code to generate figures and tables based on line-list data

`DataAnalysis_ModelledEstimates.Rmd` - R code to generate figures and tables based on modelled estimates of K/O prevalence (see details of Bayesian modelling below)

`DataAnalysis_Longitudinal.Rmd` - R code to generate Figure S6, based on longitudinal line-list data from 3 sites

## Model Run File Explanation

Loads the file: `model_run_file.R`

### Overview

This R script automates the process of preparing and analysing CSV data by loading essential libraries, listing and validating files in designated directories against a specific naming convention, extracting key metadata from the filenames, and merging the data with site information to incorporate region data; it then expands the dataset to include every combination of locus, subgroup, and site—to account for zero recordings of a given locus at a given site—and fits a Bayesian generalised linear mixed effects model using the `brms` package with a binomial link function, generic priors, and control parameters, before saving the resulting model as an RDS file with a name derived from the extracted metadata.

### Dependencies

Ensure you have the following R packages installed: -

`brms`

`dplyr`

`tidyr`

`ggplot2`

`ggrepel`

`tibble`

`stringr`

`tidyverse`

### File and Directory Structure

-   **Data Files**: CSV files are located in `data_core/` or `data_LOO/`. A `site_info.csv` file must also reside in that directory.
-   **Models**: The resulting models will be saved in `models_core/` or `models_LOO/`, depending on whether the script detects “core” or “LOO” in the file path.

The file begins by loading the above R packages needed for Bayesian modeling and data manipulation, then:

1.  **Lists and Validates Data Files:**
    -   Searches the data_core/ (or data_LOO/) folder for CSV files.
    -   Checks whether these files match a particular naming pattern (stored in the variable "pattern").
    -   Flags any files that don’t conform to this format.
2.  **Metadata Extraction:**
    -   For each file, it extracts key information:
        -   The type (e.g., "Full", "Carba", "ESBL", "Fatal").
        -   Whether the file is "core" or "LOO" (Leave-One-Out).
        -   The date stamp (YYYYMMDD).
        -   The number of days (28 or 365).
        -   The data type (e.g., "ALL" or "min10").
        -   The antigen type (e.g., "K", "O", or "OlocusType").
    -   The purpose of this section is to ensure models are saved appropriately corresponding to the data they were run on.
3.  **Data Preparation:**
    -   Reads each CSV file into R.
    -   Merges each file with site-level metadata from a site_info.csv file to include region data.
    -   Expands the dataset to ensure every combination of locus, subgroup, and site is explicitly listed—even if there are zero events.
4.  **Bayesian Modeling:**
    -   Uses the brms package to fit a logistic regression model:

        `event | trials(n) ~ 0 + locus + subgroup + (1 | Site) + (1 | locus:subgroup)`

    -   Employs a binomial link function, generic priors, and runs the model with a 6000 iterations and specific control parameters (adapt_delta = 0.999999, max_treedepth = 55).
5.  **Saving:**
    -   Checks model diagnostics such as maximum Rhat values and minimum effective sample size ratios.
    -   Saves the fitted model as an RDS file, naming it based on the extracted metadata.

## Model Process File Explanation

Loads the file: `model_process_file.R`

### Overview

This R script automates the process of reloading and processing Bayesian models previously fitted by the `model_run_file.R`. It begins by loading the required libraries and listing CSV data files from the designated directories, validating these files filenames. It then extracts key metadata from the filenames—including type, purpose (core or LOO), date stamp, duration, data type, and antigen—and uses this information to construct the corresponding model filenames. The script then loads the pre-fitted models from the appropriate model directory. Using the loaded model, the script generates posterior predictions via the `posterior_epred` function and computes global prevalence estimates for each locus as well as regional estimates for each (locus, subgroup) pair by aggregating the posterior draws into summary statistics (mean, median, and 95% credible intervals). Finally, the results are saved as CSV files (with the global posterior data compressed) in an output folder named according to the purpose (core or LOO).

### Dependencies

Ensure you have the following R packages installed:

`dplyr`\
`stringr`\
`tidyr`\
`brms`

### File and Directory Structure

-   **Data Files**: CSV files are located in `data_core/` or `data_LOO/`. A `site_info.csv` file must also reside in the same directory.
-   **Models**: Pre-fitted models are loaded from `models_core/` or `models_LOO/`, depending on whether the file path indicates “core” or “LOO”.
-   **Outputs**: Processed posterior data and summary estimates are saved in `outputs_core/` or `outputs_LOO/`, based on the detected purpose.

The file begins by loading the necessary libraries for data manipulation, string operations, and Bayesian modelling, then:

1.  **Lists and Validates Data Files:**\
    Searches the `data_core/` (or `data_LOO/`) folder for CSV files, checks them against a predefined naming pattern, and flags any files that do not conform.

2.  **Metadata Extraction:**\
    For each file, it extracts key information:

    -   The type (e.g., "Full", "Carba", "ESBL", "Fatal").
    -   Whether the file is "core" or "LOO" (Leave-One-Out).
    -   The date stamp (YYYYMMDD).
    -   The number of days (28 or 365).
    -   The data type (e.g., "ALL" or "min10").
    -   The antigen type (e.g., "K", "O", or "OlocusType").

3.  **Model Loading and Data Processing:**\
    For each file and for each data type ("adj" and "raw"), the script constructs the appropriate model filename from the metadata and loads the corresponding pre-fitted model. It then reads the prevalence data and site information, merges them to incorporate regional data, and expands the dataset to ensure every combination of locus, subgroup, and site is present.

4.  **Global Processing:**\
    The script uses the loaded model to generate posterior predictions with the `posterior_epred` function, and then calculates global prevalence estimates for each locus by aggregating the predictions (computing summary statistics such as mean, median, and 95% credible intervals). These global results are saved as CSV files, with the posterior draws compressed as .gz files.

5.  **Regional Processing:**\
    Similarly, it computes posterior prevalence estimates for each (locus, subgroup) pair, aggregates the results into summary statistics, and saves both the detailed posterior draws and the regional summary estimates as CSV files, with the posterior draws compressed as .gz files.

## 

## Diagnostic Table File Explanation

Loads the file: `diagnostics_table.R`

### Overview

This R script automates the process of assessing the diagnostic performance of pre-fitted Bayesian models and summarising key model diagnostics. It begins by loading the required libraries for data manipulation, string operations, and Bayesian modeling. For each file, the script extracts essential metadata—such as subset type, date stamp, duration (28 or 365 days), data type, and antigen—from the filename, and constructs the corresponding model filename for both adjusted ("adj") and raw ("raw") analyses. The script then loads the appropriate pre-fitted model from the directory, computes diagnostic metrics (including maximum Rhat values, minimum effective sample size ratios, and the number of divergent transitions), and appends these values to a results table. Finally, the script further processes the results table by extracting additional metadata from the model filenames and selecting relevant columns, before saving the summarised diagnostics as a CSV file (`results_table.csv`).

### Dependencies

Ensure you have the following R packages installed:

**`dplyr`**\
**`stringr`**\
**`tidyr`**\
**`brms`**

### File and Directory Structure

-   **Data Files**: CSV files are located in `data_paper/`. A `site_info.csv` file must also reside in that directory.
-   **Models**: Pre-fitted models are loaded from the `models_paper/` directory.
-   **Results**: The summarised model diagnostics are saved in a CSV file (`results_table.csv`).

The file begins by loading the necessary libraries for data manipulation and Bayesian modelling, then:

1.  **Lists and Validates Data Files:**\
    The script lists all CSV files in the `data_paper/` directory (excluding `site_info.csv`), validates them against a predefined naming pattern, and warns if any files do not match the required format.

2.  **Initialization of the Results Table:**\
    An empty results table is initialized to store key diagnostic metrics: model filename, maximum Rhat, minimum effective sample size ratio, and the number of divergent transitions.

3.  **Model Loading and Diagnostic Extraction:**\
    For each file, the script extracts metadata (subset, date stamp, duration, data type, antigen, and study information if available) from the filename and constructs the corresponding model filename for both "adj" and "raw" data types. It then loads each pre-fitted model from the `models_paper/` directory and computes diagnostic metrics, including the maximum Rhat, the minimum effective sample size ratio, and the total count of divergent transitions from the sampler diagnostics. These values are added as new rows to the results table.

4.  **Post-Processing of the Results Table:**\
    After processing all models, the results table is further refined by extracting additional metadata (such as subset, type, days, data type, and antigen) from the model filenames, and selecting only the relevant columns to create a concise summary of the model diagnostics.

5.  **Saving the Results:**\
    The final summarized results are printed to the console and saved as a CSV file (`results_table.csv`) for further review and analysis.
