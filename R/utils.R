#' Conservative variance estimator for finely stratified designs (Fogarty 2018, 2023)
#'
#' @param strat_ns Integer vector of stratum sizes n_s.
#' @param strat_ests Numeric vector of stratum-specific estimates \eqn{\hat{\tau}_s}.
#' @return A single numeric value: conservative variance estimate.
#' @export
fine_strat_var_est <- function(strat_ns, strat_ests) {
  if (length(strat_ns) != length(strat_ests)) {
    stop("strat_ns and strat_ests must have the same length.")
  }
  if (any(is.na(strat_ns)) || any(is.na(strat_ests))) {
    stop("strat_ns and strat_ests must not contain NA.")
  }
  if (any(strat_ns <= 0)) {
    stop("All stratum sizes in strat_ns must be positive.")
  }
  strat_ns <- as.numeric(strat_ns)
  strat_ests <- as.numeric(strat_ests)
  N <- sum(strat_ns)
  B <- length(strat_ns)
  W <- diag(B * (strat_ns / N))
  Q <- matrix(B * (strat_ns / N), ncol = 1)
  S_w <- sum(Q^2)
  H_Q <- Q %*% (1 / S_w) %*% t(Q)
  leverages <- pmin(pmax(diag(H_Q), 0), 1)
  denom <- sqrt(pmax(1 - leverages, .Machine$double.eps))
  y_Gamma <- matrix(strat_ests / denom, ncol = 1)
  I <- diag(nrow(H_Q))
  var_hat <- as.numeric((1 / B^2) * t(y_Gamma) %*% W %*% (I - H_Q) %*% W %*% y_Gamma)
  return(var_hat)
}

#' Worst-case IPW contribution for a single finely stratified set (Fogarty 2023)
#'
#' @param z 0/1 vector of treatment for one stratum (exactly 1 treated OR 1 control).
#' @param y Numeric vector of outcomes for the same units.
#' @param Gamma Sensitivity parameter \eqn{\Gamma \ge 1}.
#' @param tau_h Null value for the ATE (default 0).
#' @param alternative "greater" or "less" for the one-sided test direction.
#' @return Numeric: worst-case IPW-weighted (diff-in-means - tau_h)/|\Omega_s|.
#' @export
worst_case_IPW <- function(z, y, Gamma, tau_h = 0, alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  if (!is.numeric(y) || length(y) != length(z)) stop("y must be numeric and same length as z.")
  if (any(is.na(z)) || any(is.na(y))) stop("z and y must not contain NAs.")
  if (!all(z %in% c(0, 1))) stop("z must be a 0/1 vector.")
  if (!is.numeric(Gamma) || length(Gamma) != 1 || Gamma < 1) stop("Gamma must be a numeric scalar >= 1.")
  n_treated <- sum(z)
  n_control <- sum(1 - z)
  N <- length(z)
  if (!(n_treated == 1L || n_control == 1L)) {
    stop("Each matched set must have exactly 1 treated OR exactly 1 control.")
  }
  card_Omega <- choose(N, n_treated)
  diff_means <- mean(y[z == 1]) - mean(y[z == 0])
  d <- (diff_means - tau_h)
  prob_lb <- 1 / (Gamma * (N - 1) + 1)
  prob_ub <- Gamma / ((N - 1) + Gamma)
  if (alternative == "greater") {
    w <- if (d > 0) 1 / prob_ub else 1 / prob_lb
  } else {
    w <- if (d < 0) 1 / prob_ub else 1 / prob_lb
  }
  (1 / card_Omega) * d * w
}
