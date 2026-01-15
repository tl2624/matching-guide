# Load senstrat for evall() and senstrat()
library(senstrat)

#------------------------------------------------------------------
# Build table of null expectations and variances by stratum and j
#
# Intuition:
#   Within each stratum, Rosenbaum & Krieger (1990) show that worst-case
#   hidden bias vectors take a simple form: j leading 1s followed by 0s
#   after ordering units by outcome.
#
#   evall() computes, for each such j, the null expectation and variance
#   contribution of that stratum under a given Gamma.
#------------------------------------------------------------------
build_ev_table <- function(sc, z, st, Gamma, method = "RK") {
  
  sc <- as.numeric(sc)    # Outcome scores
  z  <- as.integer(z)     # Treatment indicators (0/1)
  st <- as.factor(st)     # Matched set labels
  
  sc_split <- split(sc, st)   # Split outcomes by stratum
  z_split  <- split(z,  st)   # Split treatment by stratum
  
  ev_list <- lapply(seq_along(sc_split), function(s) {
    
    sc_s <- sc_split[[s]]     # Outcomes in stratum s
    z_s  <- z_split[[s]]      # Treatment in stratum s
    
    # evall(): for each j = 1,...,n_s - 1, compute
    #   expectation and variance under Gamma
    ev_s <- evall(sc_s, z_s, g = Gamma, method = method)
    
    data.frame(
      stratum = levels(st)[s],
      j       = ev_s$m,       # j = number of leading 1s
      mu_sj   = ev_s$expect,  # null expectation contribution
      var_sj  = ev_s$var      # null variance contribution
    )
  })
  
  do.call(rbind, ev_list)
}

#------------------------------------------------------------------
# Separable approximation: choose j in each stratum independently
#
# Intuition:
#   The separable approximation picks, in each stratum, the j that
#   makes the null expectation as large as possible.
#   If there is a tie, it chooses the j with larger variance.
#
#   This gives a single "baseline" configuration b that is close
#   to the global worst case and serves as the expansion point.
#------------------------------------------------------------------
separable_index <- function(ev_tab) {
  
  split(x = ev_tab,
        f = ev_tab$stratum) |>
    lapply(FUN = function(df_s) {
      # Sort j values from most to least adverse: higher expectation first, then higher variance
      df_s_ord <- df_s[order(df_s$mu_sj,
                             df_s$var_sj,
                             decreasing = TRUE), ]
      df_s_ord$j[1]           # b_s for this stratum
    }) |>
    unlist()
}

#------------------------------------------------------------------
# lambda(j) for a given configuration of j values
#
# Intuition:
#   For a fixed configuration j (one j per stratum), lambda(j)
#   measures how far the observed statistic exceeds an upper
#   Normal cutoff under the null.
#------------------------------------------------------------------
lambda_from_indices <- function(ev_tab, indices, T_obs, alpha = 0.05) {
  
  ev_tab$chosen <- with(ev_tab, j == indices[stratum])
  
  mu_null  <- sum(ev_tab$mu_sj[ev_tab$chosen])    # total null expectation
  var_null <- sum(ev_tab$var_sj[ev_tab$chosen])   # total null variance
  
  crit <- qnorm(1 - alpha)                        # upper one-sided cutoff
  
  mu_null + crit * sqrt(var_null) - T_obs
}

#------------------------------------------------------------------
# Build zeta table at the separable configuration
#
# Intuition:
#   zeta_sj represents the slope of the lambda function in stratum s
#   when we move from the separable j to another j.
#
#   These slopes define the tangent (linear) upper bound.
#------------------------------------------------------------------
zeta_table <- function(ev_tab, b_indices, alpha = 0.05) {
  
  ev_tab$is_b <- with(ev_tab, j == b_indices[stratum])
  
  var_b <- sum(ev_tab$var_sj[ev_tab$is_b])        # total variance at separable point
  crit  <- qnorm(1 - alpha)
  
  ev_tab$zeta_sj <- ev_tab$mu_sj +
    crit * ev_tab$var_sj / (2 * sqrt(var_b))
  
  ev_tab
}

#------------------------------------------------------------
# Taylor / tangent-line approximation at configuration j
#   - Rosenbaum (2018), equation (10)
#   - lambda_Taylor(j) = lambda(b) + sum_s [ zeta_{s,j_s} - zeta_{s,b_s} ]
#------------------------------------------------------------
lambda_taylor_from_indices <- function(ev_zeta_tab,
                                       indices,
                                       b_indices,
                                       T_obs,
                                       alpha = 0.05) {
  # ev_zeta_tab : output from zeta_table(), contains zeta_sj by stratum and j
  # indices     : configuration j (one j per stratum) at which to evaluate tangent line
  # b_indices   : separable configuration b (one j per stratum)
  # T_obs       : observed sum statistic
  # alpha       : one-sided test level
  
  # Exact lambda at the separable configuration b
  lambda_b <- lambda_from_indices(
    ev_tab  = ev_zeta_tab,
    indices = b_indices,
    T_obs   = T_obs,
    alpha   = alpha
  )
  
  # Extract zeta_{s,j_s} for the candidate configuration j
  zeta_at_j <- ev_zeta_tab$zeta_sj[
    ev_zeta_tab$j == indices[ev_zeta_tab$stratum]
  ]
  
  # Extract zeta_{s,b_s} for the separable configuration b
  zeta_at_b <- ev_zeta_tab$zeta_sj[
    ev_zeta_tab$j == b_indices[ev_zeta_tab$stratum]
  ]
  
  # Tangent-line approximation:
  # lambda_Taylor(j) = lambda(b) + sum over strata of (zeta_{s,j_s} - zeta_{s,b_s})
  lambda_b + sum(zeta_at_j - zeta_at_b)
}