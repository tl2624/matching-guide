# ---- chunk ----
library(knitr)  # For knitting R Markdown and including external figures or files

# Insert the matching pipeline flowchart (PDF) from the project’s fig/ directory
include_graphics(path = "fig/matching_pipeline_flowchart.pdf")


# ---- chunk ----
# Define base URL for the Matching Guide GitHub repository
base_url <- "https://raw.githubusercontent.com/tl2624/matching-guide/main"

# Load the cleaned dataset
data <- readRDS(url(paste0(base_url, "/data/peace_pre_match.rds")))


# ---- chunk ----
# Define character vector of the 9 covariate names in the dataset
covs <- c("lwdeaths", "lwdurat", "ethfrac", "pop", "lmtnest", "milper", "bwgdp",
          "bwplty2", "region")


# ---- chunk ----
# Install "dplyr" package (only run if you don't already have it installed)
# Install.packages("dplyr")

# Load dplyr package for data manipulation (mutate, group_by, summarize, etc.)
library(dplyr)

# Convert categorical variable 'region' into 0/1 dummy indicators
data <- data |>  # Pipe (|>) to pass left-hand result into next function call
  mutate(
    eeurop   = ifelse(test = region == "eeurop",   yes = 1, no = 0),
    lamerica = ifelse(test = region == "lamerica", yes = 1, no = 0),
    asia     = ifelse(test = region == "asia",     yes = 1, no = 0),
    ssafrica = ifelse(test = region == "ssafrica", yes = 1, no = 0),
    nafrme   = ifelse(test = region == "nafrme",   yes = 1, no = 0)
  )

# Remove 'region' and replace with dummy indicators
expanded_covs <- c(
  # Setdiff() returns elements in 'x' that are not in 'y'
  setdiff(x = covs, y = "region"),  # Drops "region" from the covariate list
  "eeurop", "lamerica", "asia", "ssafrica", "nafrme"
)


# ---- chunk ----
# Extract covariate values for Liberia and Guinea-Bissau (cname = country name)
liberia <- data[data$cname == "Liberia", expanded_covs]
guinea_bissau <- data[data$cname == "Guinea-Bissau", expanded_covs]

# Compute Euclidean distance between the two countries on these covariates
sqrt(sum((liberia - guinea_bissau)^2))  # Display the Euclidean distance


# ---- chunk ----
# Create a formula: UN (treatment indicator) ~ covariates
# Note: we keep "region" in covs as a factor
cov_fmla <- reformulate(termlabels = covs,
                        response = "UN")

# Install optmatch if not already installed
# Install.packages("optmatch")

# Load optmatch, which provides the match_on() function
library(optmatch)

# Compute Euclidean distance matrix between treated (UN = 1) and control (UN = 0)
dist_mat_euc <- match_on(x = cov_fmla,                 # Formula for covariates
                         data = data,                  # Dataset used
                         standardization.scale = NULL, # No rescaling of covariates
                         method = "euclidean")         # Use Euclidean distance

# Add country names (cname) as row/column labels for clarity
dimnames(dist_mat_euc) <- list(data$cname[data$UN == 1], data$cname[data$UN == 0])

# Display a submatrix of distances: treated units 11–15 vs control units 27–30
round(x = dist_mat_euc[11:15, 27:30],
      digits = 2)  # Number of decimal places to round


# ---- chunk ----
# Compute Mahalanobis distance matrix between treated (UN = 1) and control (UN = 0)
dist_mat_mah <- match_on(
  x = cov_fmla,
  data = data,
  standardization.scale = NULL,
  method = "mahalanobis" # Use Mahalanobis distance
)


# ---- chunk ----
# Formula for UN ~ covariates (excluding "region")
psm_cov_fmla <- reformulate(termlabels = setdiff(x = covs, y = "region"),
                            response = "UN")

# Fit logistic regression for propensity score model
psm <- glm(
  formula = psm_cov_fmla,            # Treatment ~ covariates
  family = binomial(link = "logit"), # Logistic regression (logit link)
  data = data                        # Dataset used for model fitting
)


# ---- chunk ----
# Extract logit propensity scores (linear predictors from fitted model)
lin_cov_inds <- psm$linear.predictors  # Same as model.matrix(psm) %*% coef(psm)


# ---- chunk ----
# Extract estimated propensity scores (predicted probabilities of UN = 1)
p_scores <- psm$fitted.values

# Convert propensity scores to log-odds (logit scale)
log(p_scores / (1 - p_scores))

# Check that p_scores equals logistic(lin_cov_inds)
all.equal(p_scores, 1 / (1 + exp(-lin_cov_inds)))


# ---- chunk ----
# Install "ggplot2" (only run if not already installed)
# Install.packages("ggplot2")

# Load ggplot2 package for visualization
library(ggplot2)

# Boxplot of linear covariate index by treatment status
ggplot(data = data,  # Dataset
       mapping = aes(x = as.factor(UN),    # Treatment indicator (0/1 as factor)
                     y = lin_cov_inds)) +  # Linear covariate index (from logistic regression)
  geom_boxplot() +                         # Draw boxplots
  xlab(label = "UN intervention") +        # X-axis label for treatment
  ylab(label = "Linear covariate index") + # Y-axis label for lin cov index
  theme_bw() +                             # Apply black-and-white theme
  coord_flip()                             # Flip axes for readability


# ---- chunk ----
# Create distance structure: 0 if units are in the same region, Inf otherwise
em_region <- exactMatch(x = UN ~ region,
                        data = data)


# ---- chunk ----
# Euclidean distance on Polity score (bwplty2) with caliper = 2
# Pairs differing by >2 are set to Inf
euc_dist_polity_cal_2 <- match_on(x = UN ~ bwplty2,
                                  caliper = 2,  # Set caliper
                                  data = data,
                                  standardization.scale = NULL,
                                  method = "euclidean")


# ---- chunk ----
# Create overall distance matrix by element-wise addition of two distance matrices
overall_dist_mat <- em_region + euc_dist_polity_cal_2


# ---- chunk ----
# Apply a caliper of width 3 to the polity Euclidean distance matrix
euc_dist_polity_cal_3 <- match_on(x = UN ~ bwplty2,
                                  caliper = 3,
                                  data = data,
                                  standardization.scale = NULL,
                                  method = "euclidean")

# Combine regional exact match distance with polity distance
em_region + euc_dist_polity_cal_3


# ---- chunk ----
# Full matching using overall distance matrix; allows 0.5 - 2 controls per treated
fullmatch(x = overall_dist_mat, min.controls = 0.5, max.controls = 2, data = data)


# ---- chunk ----
# Add linear predictors from logistic regression (psm$linear.predictors) to dataset
data$logit_p_score <- lin_cov_inds

# Population standard deviation of logit_p_score (divides by n, not n - 1)
pop_sd_logit <- sqrt(mean((data$logit_p_score - mean(data$logit_p_score))^2))

# Distance matrix from propensity score (logit of estimated treatment probability)
ps_mat <- match_on(x = UN ~ logit_p_score,
                   caliper = 0.5 * pop_sd_logit,
                   data = data,
                   standardization.scale = NULL,
                   method = "euclidean")


# ---- chunk ----
# Rank-based Mahalanobis distance on covariates
# Covs was defined earlier as the set of covariate names; here we drop "region"
rank_mah_mat <- match_on(
  x      = reformulate(termlabels = setdiff(x = covs, y = "region"),
                       response   = "UN"),
  data   = data,
  standardization.scale = NULL,
  method = "rank_mahalanobis" # Use rank-based Mahalanobis distance
)


# ---- chunk ----
# Compute Euclidean distance matrix for ethnic fractionalization
eth_mat <- match_on(
  x = UN ~ ethfrac,
  caliper = 0.35,
  data = data,
  standardization.scale = NULL,
  method = "euclidean"
)

# Compute Euclidean distance matrix for logged GDP per capita
bwgdp_mat <- match_on(
  x = UN ~ bwgdp,
  caliper = 2,
  data = data,
  standardization.scale = NULL,
  method = "euclidean"
)


# ---- chunk ----
# Full matching on PS + rank-based Mahalanobis + separate Euclidean distances
# (ethfrac, bwgdp) + region exact match
fm <- fullmatch(
  x            = ps_mat + rank_mah_mat + eth_mat + bwgdp_mat + em_region,
  data         = data,
  max.controls = 4  # Up to 4 controls per treated; min.controls = 0 by default
)


# ---- chunk ----
# Effective sample size of matched sets
effectiveSampleSize(fm)


# ---- chunk ----
# Summarize matched sets (set sizes, structure) and report effective sample size
summary(fm)


# ---- chunk ----
# Add matched set ID to data for each unit
data$fm <- fm

# Look at one matched set ("ssafrica.3") for illustration
data |>
  filter(fm == "ssafrica.3") |>  # Keep only units in set "ssafrica.3"
  # Display selected variables
  select(cname, UN, region, logit_p_score, ethfrac, bwgdp, bwplty2)


# ---- chunk ----
# Install "RItools" package (only run if you don't already have it installed)
# Install.packages("RItools")

# Load RItools package for balance diagnostics (xBalance)
library(RItools)

# Compute covariate balance statistics with xBalance
cov_bal <- xBalance(
  fmla = cov_fmla,  # Formula: treatment ~ covariates
  strata = list(
    unstrat = NULL, # Overall balance without stratification
    fm = ~ fm       # Assess balance within matched sets (stratify by fm)
  ),
  data = data,
  report = c("adj.means", "std.diffs")
  # Return adjusted means (weighted by effective sample size)
  # And standardized differences (mean differences scaled by pooled SD)
)


# ---- chunk ----
# Install packages (only run these lines if not already installed)
# Install.packages("tibble")
# Install.packages("kableExtra")

# Load tibble for cleaner data frame printing/handling
library(tibble)

# Load kableExtra for formatting tables for LaTeX/HTML output
library(kableExtra)

# Extract xBalance 3-D results:
# [vars, stats (Control,Treatment,std.diff,p), strata (unstrat, fm)]
arr  <- cov_bal$results
vars <- dimnames(arr)$vars

# Get readable labels from data (haven-style), else fall back to variable name
label_from_data <- function(v) {
  if (v %in% names(data)) {
    lb <- attr(data[[v]], "label", exact = TRUE)
    if (!is.null(lb) && nzchar(lb)) return(lb)
  }
  v
}

# Map region dummy names → readable labels
region_map <- c(
  regioneeurop   = "Eastern Europe",
  regionlamerica = "Latin America",
  regionnafrme   = "North Africa & Middle East",
  regionssafrica = "Sub-Saharan Africa",
  regionasia     = "Asia"
)

# Final covariate labels
# Use region_map if present; else dataset label; else raw name
cov_labels <- sapply(
  X = vars,
  FUN = function(v) { ifelse(test = v %in% names(region_map),
                             yes = region_map[[v]],
                             no = label_from_data(v))
  },
  USE.NAMES = FALSE
)

# Helper: slice a 2-D matrix (vars × stats) for one stratum
slice_mat <- function(a, stratum) {
  a[, c("Control","Treatment","std.diff","p"),
    stratum,
    drop = FALSE][,,1, drop = TRUE]
}

# Before/after matrices
before <- slice_mat(arr, "unstrat")
after  <- slice_mat(arr, "fm")

# Build tibble with distinct internal names
cov_tab <- tibble(
  Covariate      = cov_labels,
  Control_before = round(x = before[, "Control"], digits = 2),
  Treated_before = round(x = before[, "Treatment"], digits = 2),
  StdDiff_before = before[, "std.diff"],
  Control_after  = round(x = after[,  "Control"], digits = 2),
  Treated_after  = round(x = after[,  "Treatment"], digits = 2),
  StdDiff_after  = after[,  "std.diff"],
  p_before = before[, "p"],  # Keep for star annotation
  p_after        = after[,  "p"]
)

# Star standardized differences when p <= 0.05
fmt_sd <- function(x, p) {ifelse(test = p <= 0.05,
                                 yes = sprintf("%.2f*", x),
                                 no =  sprintf("%.2f", x))
}

# Apply formatting; drop helper p-values
cov_tab <- cov_tab |>
  mutate(
    StdDiff_before = fmt_sd(StdDiff_before, p_before),
    StdDiff_after  = fmt_sd(StdDiff_after,  p_after)
  ) |>
  select(-p_before, -p_after)

# Print LaTeX table:
# - grouped headers: Before matching / After matching
# - printed column names omit "(Before)/(After)" (clean labels)
kbl(
  cov_tab,
  booktabs  = TRUE,
  align     = c("l","c","c","c","c","c","c"),
  col.names = c("Covariate", "Control mean", "Treated mean", "Std. diff",
                "Control mean", "Treated mean", "Std. diff"),
  linesep = ""  # No extra spacing after rows
) |>
  add_header_above(c(" " = 1, "Before matching" = 3, "After matching" = 3),
                   bold = TRUE) |>
  kable_styling(latex_options = c("hold_position", "scale_down"))


# ---- chunk ----
xBalance(
  fmla = cov_fmla,
  strata = list(fm = ~ fm),  # Assess balance within matched sets (stratify by fm)
  data = data,
  report = "chisquare.test"  # Return chi-square test results for overall balance
)


# ---- chunk ----
# ---- define matched set and subset data ----
set_id <- "ssafrica.3"  # Matched set name
sdat <- data |> filter(fm == set_id) |> select(cname, UN)
n <- nrow(sdat)  # Number of total units
m <- sum(sdat$UN)  # Number of units treated

# ---- enumerate all possible assignments ----
# Install "randomizr" (only run if not already installed)
# Install.packages("randomizr")
# Load randomizr for generating random assignments
library(randomizr)
Z <- obtain_permutation_matrix(declaration = declare_ra(N = n, m = m))
# Label assignment columns
colnames(Z) <- paste0("Assignment ", seq_len(ncol(Z)))  

# ---- build LaTeX table manually (no unwanted vertical lines) ----
assign_cols <- ncol(Z)
# One vline only after first col
preamble <- paste0("l|", paste(rep("c", assign_cols), collapse = ""))
header <- paste(c("", colnames(Z)), collapse = " & ")  # Blank header for unit

# ---- create one LaTeX row per unit ----
row_lines <- vapply(
  seq_len(n),
  function(i) {
    cells <- as.integer(Z[i, ])
    paste(c(sdat$cname[i], cells), collapse = " & ")
  },
  character(1)
)

# ---- wrap table in center environment ----
tex <- c(
  "\\vspace{1em}",  # Space before table
  "\\begin{center}",
  paste0("\\begin{tabular}{", preamble, "}"),
  "\\hline",
  header, " \\\\",
  "\\hline",
  paste0(row_lines, " \\\\"),
  "\\hline",
  "\\end{tabular}",
  "\\end{center}",
  "\\vspace{1em}"  # Space after table
)

# ---- print final LaTeX code ----
cat(paste(tex, collapse = "\n"))


# ---- chunk ----
# Keep only rows assigned to a matched set (drop NA in fm)
data_matched <- filter(.data = data, !is.na(fm))

# Null hypothesis value
tau_h <- 0

# Reconstruct outcomes under sharp null (tau_h = 0)
data_matched <- data_matched |>
  mutate(ldur_tilde = ldur - tau_h * UN)


# ---- chunk ----
# Load the hm_stat_rescale() function from the GitHub repo
# Base_url (defined earlier as
# "https://raw.githubusercontent.com/tl2624/matching-guide/main")
# Points to the main GitHub repo URL
source(paste0(base_url, "/R/hm_stat_rescale.R"))

# Apply the rescaling function: adds a new column (.hm_scaled)
# And returns the full matched dataset with this rescaled outcome
data_matched <- hm_stat_rescale(
  data = data_matched,  # Set dataset containing matched observations
  outcome = ldur_tilde, # Set outcome variable to be rescaled within matched sets
  treat = UN,           # Set name of treatment indicator variable
  strata = fm           # Set name of matched strata (block) variable
)

# Observed HM-weighted diff-in-means statistic
obs_stat <- sum(data_matched$ldur_tilde_hm_scaled[data_matched$UN == 1])
# Equivalent to xBalance(fmla = UN ~ ldur, strata = list(fm = ~ fm),
# Data = data_matched, report = "adj.mean.diffs")


# ---- chunk ----
# For each matched set (fm), record:
# N   = total units in the set
# M = number treated (UN == 1)
block_ns <- data_matched |>
  group_by(fm) |>    # Group results by key variables
  summarise(         # Aggregate to one row per group
    n = n(),         # Row count per group
    m = sum(UN),
    .groups = "drop" # Drop grouping after summarise
  )

# Total possible treatment assignments = product of binomial coefficients
# (choose n_s units for treatment in each set and multiply across sets)
prod(choose(n = block_ns$n, k = block_ns$m))


# ---- chunk ----
# Install "randomizr" (only run if not already installed)
# Install.packages("randomizr")

# Load randomizr for generating random assignments
library(randomizr)

exact_assigns <- obtain_permutation_matrix(
  declaration = declare_ra(   # Declare assignment procedure
    N = nrow(data_matched),   # Total number of units
    blocks = data_matched$fm, # Matched set membership
    block_m = block_ns$m      # Number treated in each set
  ),
  # Total number of feasible assignments across all matched sets
  maximum_permutations = prod(choose(n = block_ns$n, k = block_ns$m))
)


# ---- chunk ----
# Set RNG seed for reproducibility
set.seed(11242017)

# Generate permutation matrix of treatment assignments
# Each column = one possible assignment consistent with block structure
sim_assigns <- obtain_permutation_matrix(
  declaration = declare_ra(
    N = nrow(data_matched),
    blocks = data_matched$fm,
    block_m = block_ns$m
  ),
  maximum_permutations = 10^4  # Cap at 10,000 random draws
)


# ---- chunk ----
# Randomization distribution under sharp null of no effect:
# Apply sum statistic to each assignment column in 'assigns'
# Outcome has been transformed so that
# Sum statistic = harmonic-mean weighted diff in means
sim_sharp_null_dist <- apply(
  X = sim_assigns,    # Matrix of treatment assignments
  MARGIN = 2,         # Iterate over columns (assignments)
  FUN = function(x) {
    # Sum transformed outcomes among treated
    sum(data_matched$ldur_tilde_hm_scaled[x == 1])
  }
)


# ---- chunk ----
# One-sided, upper p-value: proportion of simulated randomization stats >= observed
round(x = mean(sim_sharp_null_dist >= obs_stat), digits = 4)


# ---- chunk ----
# Randomization distribution under sharp null of no effect:
# Apply sum statistic to each assignment column in 'assigns'
# Outcome has been transformed so that
# Sum statistic = harmonic-mean weighted diff in means
exact_sharp_null_dist <- apply(
  X = exact_assigns,
  MARGIN = 2,
  FUN = function(x) {
    sum(data_matched$ldur_tilde_hm_scaled[x == 1])
  }
)

# Exact one-sided, upper p-value: proportion of randomization stats >= observed
round(x = mean(exact_sharp_null_dist >= obs_stat), digits = 4)


# ---- chunk ----
# Install "senstrat" package (only run if you don't already have it installed)
# Install.packages("senstrat")

# Load senstrat for computing stratum-level null expectations/variances
# (Rosenbaum & Krieger 1990) and later sensitivity analysis
library(senstrat)

# Compute per-block null expectations and variances
per_block_moms <- data_matched |>
  group_by(fm) |>
  summarize(
    expect   = ev(
      sc = ldur_tilde_hm_scaled, # Transformed outcomes for the stratum
      z = UN,                    # Treatment indicator
      m = 1,                     # Number of "1"s in vector of hidden confounder
      # Irrelevant here since Gamma = 1
      g = 1,                     # Sensitivity parameter Gamma
      method = "RK"              # Use formula from Rosenbaum and Krieger (1990)
    )$expect,
    variance = ev(
      sc     = ldur_tilde_hm_scaled,
      z      = UN,
      m = 1,
      g      = 1,
      method = "RK"
    )$vari,
    .groups  = "drop"
  )

# Sum across blocks to get overall null expectation and variance
null_ev  <- sum(per_block_moms$expect)
null_var <- sum(per_block_moms$variance)

# Standardized test statistic and one-sided Normal p-value (upper tail)
norm_upper_p_value <- pnorm(
  q = (obs_stat - null_ev) / sqrt(null_var), # Standardized statistic
  lower.tail = FALSE                         # Compute upper-tail probability
)


# ---- chunk ----
# Install blkvar (only run if not already installed)
# Install.packages("remotes")
# Remotes::install_github("lmiratrix/blkvar")

# Load blkvar for block randomization variance estimators
library(blkvar)

# Compute results with hybrid_p method
res <- block_estimator(
  Yobs = ldur,         # Observed outcomes
  Z = UN,              # Treatment indicator
  B = fm,              # Block (matched set) membership
  data = data_matched, # Dataset
  method = "hybrid_p"  # P-value method
)

# Extract variance estimate
res$var_est


# ---- chunk ----
# Load the fine_strat_var_est() function from the GitHub repo
source(paste0(base_url, "/R/fine_strat_var_est.R"))


# ---- chunk ----
# Compute stratum sizes and stratum-specific differences in means
strat_stats <- data_matched |>
  group_by(fm) |>
  summarize(
    n = n(),                               # Stratum size
    diff_in_means = mean(ldur[UN == 1L]) - # Treated mean minus
      mean(ldur[UN == 0L]),                # control mean
    .groups = "drop"
  )

# Apply Fogarty (2018/2023) variance estimator
fine_strat_var_est(
  strat_ns = strat_stats$n,              # Vector of stratum sizes
  strat_ests = strat_stats$diff_in_means # Vector of stratum estimates
)


# ---- chunk ----
pnorm(
  q = (res$ATE_hat - 0) / sqrt(res$var_est),
  lower.tail = FALSE
  # By default: mean = 0 and sd = 1 (standard normal distribution)
)


# ---- chunk ----
# Set significance (alpha) level
alpha <- 0.05

# Grid of Gamma values
Gamma_vals <- seq(from = 1, to = 1.5, by = 0.0001)

# Collect results for each Gamma
sens_results <- lapply(
  X = Gamma_vals,
  FUN = function(g){
    out <- senstrat(
      sc = data_matched$ldur_tilde_hm_scaled, # Outcome variable
      z = data_matched$UN,                    # Treatment indicator
      st = data_matched$fm,                   # Matched set (block) identifiers
      gamma = g,                              # Sensitivity parameter (Gamma)
      alternative = "greater",                # One-sided (upper-tail) test
      detail = TRUE                           # Output computation details
    )
    data.frame(
      Gamma = g,
      p_sep = as.numeric(out$Separable["P-value"]),        # Separable approx
      p_tay = as.numeric(out$LinearBoundResult["P-value"]) # Taylor approx
    )
  })

# Bind all rows into one data frame
sens_df <- do.call(what = rbind, args = sens_results)

# Smallest Gamma where the separable p-value is >= alpha
sens_value_sep <- min(sens_df$Gamma[sens_df$p_sep >= alpha])

# Smallest Gamma where the Taylor p-value is >= alpha
sens_value_tay <- min(sens_df$Gamma[sens_df$p_tay >= alpha])


# ---- chunk ----
# Install tidyr (only run if not already installed)
# Install.packages("tidyr")

# Load tidyr for reshaping data
library(tidyr)

# Convert to long format and add readable facet labels
sens_df_long <- pivot_longer(
  data = sens_df,
  cols = c(p_sep, p_tay), # Columns to pivot
  names_to = "method",    # New column for method name
  values_to = "p_value"   # New column for p-values
) |>
  mutate(
    method = factor(
      x = method,                               # Variable to convert to factor
      levels = c("p_sep", "p_tay"),             # Original category values
      labels = c("Separable approximation",     # Factor label for p_sep
                 "Taylor series approximation") # Factor label for p_tay
    )
  )

# Colorblind-friendly palette
plot_cols <- c("Separable approximation" = "#0072B2",     # blue
               "Taylor series approximation" = "#D55E00") # vermillion

# Plot sensitivity results: Gammas vs. p-values for each method
ggplot(data = sens_df_long,
       mapping = aes(x = Gamma,
                     y = p_value,
                     color = method)) +        # Line color by method
  geom_hline(yintercept = alpha,               # Y-intercept at alpha level
             linetype = "solid",               # Solid horizontal line
             color = "grey") +                 # Line color
  geom_line(linewidth = 0.5) +                 # P-value curves
  scale_color_manual(values = plot_cols,
                     name = NULL) +            # Remove legend title
  theme_bw() +
  scale_x_continuous(breaks = seq(from = 1,
                                  to = 1.5,
                                  by = 0.1)) + # X-axis ticks
  scale_y_continuous(breaks = seq(from = 0,
                                  to = 0.1,
                                  by = 0.01),
                     limits = c(0, 0.1)) +     # Y-axis ticks and limits
  xlab(expression(Gamma)) +                    # X-axis label for Gamma
  ylab("One-sided upper p-value") +            # Y-axis label for p-value
  theme(legend.position = "right")             # Legend position


# ---- chunk ----
# Load the worst_case_IPW() function from the GitHub repo
source(paste0(base_url, "/R/worst_case_IPW.R"))


# ---- chunk ----
# For each Gamma, compute one-sided weak-null p-value
weak_results_list <- lapply(
  X = Gamma_vals,
  FUN = function(G) {
    
    # Compute per-set worst-case IPW; keep as list to preserve attributes
    strat_stats <- data_matched |>
      group_by(fm) |>
      summarise(
        n = n(),
        prop_n = n / nrow(data_matched),      # Weight n_s / n
        est = list(worst_case_IPW(            # Keep as list()
          z = UN,                             # Set treatment variable
          y = ldur,                           # Set outcome variable
          Gamma = G,                          # Set Gamma value
          tau_h = 0,                          # Set null hypothesis
          alternative = "greater"             # One-sided (upper-tail) test
        )),
        .groups = "drop"
      ) |>
      mutate(
        alt = attr(est[[1]], "alternative"),  # Read attribute once
        est = as.numeric(est)
      )
    
    # Weighted statistic and variance
    num     <- sum(strat_stats$prop_n * strat_stats$est)
    var_hat <- fine_strat_var_est(
      strat_ns   = strat_stats$n,
      strat_ests = strat_stats$est
    )
    
    # Z: only if variance finite and positive; else NA
    denom <- sqrt(var_hat)
    z <- if (is.finite(denom) && denom > 0) num / denom else NA_real_
    
    # Tail from preserved attribute (same for all rows)
    alt <- strat_stats$alt[1]
    lower_tail <- if (identical(alt, "greater")) FALSE else TRUE
    
    # One-sided p-value
    pval <- if (is.na(z)) NA_real_ else pnorm(z, lower.tail = lower_tail)
    
    data.frame(Gamma = G, p_value = pval, stringsAsFactors = FALSE)
  }
)

# Bind rows: final results data.frame for plotting/reporting
weak_sens_df <- do.call(what = rbind, args = weak_results_list)

# Sensitivity value: smallest Gamma with p-value >= alpha
weak_sens_value <- min(weak_sens_df$Gamma[weak_sens_df$p_value >= alpha])


# ---- chunk ----
# Plot weak-null sensitivity results: Gammas vs. p-values
ggplot(data = weak_sens_df,
       mapping = aes(x = Gamma,
                     y = p_value)) +
  geom_hline(yintercept = alpha,
             linetype = "solid",
             color = "grey") +
  geom_line(color = "black") +  # P-value curve
  theme_bw() +
  scale_x_continuous(breaks = seq(from = 1,
                                  to = 1.5,
                                  by = 0.1)) +
  scale_y_continuous(breaks = seq(from = 0,
                                  to = 0.1,
                                  by = 0.01),
                     limits = c(0, 0.1)) +
  xlab(expression(Gamma)) +
  ylab("One-sided upper p-value")