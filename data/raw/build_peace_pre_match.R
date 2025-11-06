# ================================================================
# File: data-raw/build_peace_pre_match.R
# Purpose: Build the peace_pre_match.rds dataset used in the 
#          "Building a Design-Based Matching Pipeline" guide.
#
# Source: Official replication materials downloaded from the
#         supplementary information of the article available at:
#         https://doi.org/10.1561/100.00007051
#
# Description:
#   Reads the original Stata replication file (peace_pre_match.dta),
#   converts Stata monthly date formats to R Date objects, 
#   recodes region dummy variables into a single factor,
#   rescales the ethnic fractionalization index, and saves the
#   cleaned dataset as data/peace_pre_match.rds for reproducible use.
# ================================================================

# ---- Load required packages ----
if (!requireNamespace("haven", quietly = TRUE)) install.packages("haven")       # Read Stata files
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")       # Data manipulation
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate") # Handle dates
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")       # Read/write RDS files
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")         # Stable file paths

library(haven)
library(dplyr)
library(lubridate)
library(readr)
library(here)

# ---- Define input and output paths ----
# Paths are anchored to the repository root via the `.here` marker file.
raw_dta <- here("data", "raw", "qjps_8051_supp", "peace_pre_match.dta")
out_dir <- here("data")
out_rds <- file.path(out_dir, "peace_pre_match.rds")

# ---- Check for presence of the raw Stata file ----
if (!file.exists(raw_dta)) {
  stop(
    "Stata file not found at: ", raw_dta, "\n",
    "Download and unzip the replication materials (see DOI above)\n",
    "into data-raw/qjps_8051_supp/ before running this script."
  )
}

# Create output directory if it does not exist
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read and transform the data ----
# Read the Stata dataset (.dta)
d <- read_dta(file = raw_dta)

# Convert Stata monthly date variables and recode key fields
# - Stata %tm format counts months since Jan 1960
# - Combine regional dummies into one factor variable
# - Rescale ethnic fractionalization index from 0–100 to 0–1
base_month <- lubridate::ymd("1960-01-01")

d_out <- d %>%
  mutate(
    # Convert date0/date1 from Stata monthly numeric to actual year-month
    date0 = base_month %m+% months(as.integer(.data$date0)),
    date1 = base_month %m+% months(as.integer(.data$date1)),
    
    # Collapse regional dummy variables into one labeled factor
    region = case_when(
      .data$eeurop   == 1 ~ "eeurop",
      .data$lamerica == 1 ~ "lamerica",
      .data$asia     == 1 ~ "asia",
      .data$ssafrica == 1 ~ "ssafrica",
      .data$nafrme   == 1 ~ "nafrme",
      TRUE                ~ NA_character_
    ),
    region = factor(region),
    
    # Convert ethnic fractionalization index to 0-1 scale
    ethfrac = as.numeric(.data$ethfrac) / 100
  )

# ---- Save cleaned dataset as .rds ----
write_rds(x = d_out, file = out_rds, compress = "gz")

message(" Successfully wrote cleaned dataset to: ", out_rds)
message(" Rows: ", nrow(d_out), " | Columns: ", ncol(d_out))