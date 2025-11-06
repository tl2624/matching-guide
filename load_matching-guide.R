# load_matching-guide.R
# Helper loader for the Matching Guide repository
# Sources all R scripts from the main branch of the GitHub repo

base <- "https://raw.githubusercontent.com/tl2624/matching-guide/main/R/"
scripts <- c(
  "utils.R",
  "data_helpers.R"
)

for (f in scripts) {
  source(paste0(base, f))
}

message("Matching Guide helper functions loaded successfully.")