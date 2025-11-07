#' Rescale outcomes for harmonic-mean weighting
#'
#' Given a dataset with treatment indicator and matched-set ID,
#' this function rescales the outcome so that summing the rescaled
#' values among treated units yields the harmonic-mean–weighted
#' difference-in-means statistic.
#'
#' @param data   Data frame containing the variables.
#' @param outcome Outcome variable (unquoted, e.g. ldur_tilde).
#' @param treat  Treatment indicator (unquoted, e.g. UN; must be 0/1).
#' @param strata Matched-set ID (unquoted, e.g. fm).
#' @return A copy of `data` with an added column named
#'         `<outcome>_hm_scaled` containing the rescaled outcome values.

hm_stat_rescale <- function(data, outcome, treat, strata) {
  
  # Ensure the rlang package is available (used for tidy evaluation)
  if (!requireNamespace("rlang", quietly = TRUE)) install.packages("rlang")
  
  # Capture variable names so they can be used inside dplyr functions
  outcome <- rlang::enquo(outcome)
  treat   <- rlang::enquo(treat)
  strata  <- rlang::enquo(strata)
  
  # Drop unmatched units (NA strata)
  d <- dplyr::filter(data, !is.na(!!strata))
  
  # Compute total harmonic weight H = sum_s h_s
  H <- d |>
    dplyr::group_by(!!strata) |>
    dplyr::summarise(
      m = sum((!!treat) == 1L),             # number treated in each set
      n_minus_m = sum((!!treat) == 0L),     # number of controls in each set
      h = 1 / (1 / m + 1 / n_minus_m),      # harmonic mean weight for the set
      .groups = "drop"
    ) |>
    dplyr::summarise(H = sum(.data$h), .groups = "drop") |>
    dplyr::pull(.data$H)                    # extract scalar total weight H
  
  # Construct the new column name dynamically: <outcome>_hm_scaled
  new_col <- paste0(rlang::as_name(outcome), "_hm_scaled")
  
  # Rescale outcome so the treated sum equals the HM-weighted diff-in-means
  d |>
    dplyr::group_by(!!strata) |>
    dplyr::mutate(
      m = sum((!!treat) == 1L),             # treated count in the set
      n_minus_m = sum((!!treat) == 0L),     # control count in the set
      h = 1 / (1 / m + 1 / n_minus_m),      # set-specific harmonic weight
      sum_outcome = sum(!!outcome),         # sum of outcomes within set
      !!new_col := (                        # dynamically named new column
        (!!outcome) - h * sum_outcome / (m * n_minus_m)
      ) / H                                 # normalize by total H
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-m, -n_minus_m, -h, -sum_outcome)  # remove helper vars
}