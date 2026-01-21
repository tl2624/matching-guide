#' Construct block-respecting RK-style candidate u vectors
#'
#' @param set      Vector of block / matched-set IDs (e.g., data_matched$fm), length n.
#' @param y_scaled Numeric vector of outcomes used for ordering within blocks, length n.
#'
#' @return A list with:
#'   - U:      n x K matrix of candidate u vectors (columns), in original row order.
#'   - m_grid: data frame of (m_1, ..., m_S) choices corresponding to each column of U.
#'
construct_cand_u <- function(set, y_scaled) {
  ## --------------------------
  ## Basic checks and setup
  ## --------------------------
  
  # Coerce block IDs to a factor so strata have fixed, known levels
  block    <- as.factor(set)
  y_scaled <- as.numeric(y_scaled)
  
  # Length of vectors defines the original unit order
  n <- length(block)
  if (length(y_scaled) != n) {
    stop("set and y_scaled must have the same length.")
  }
  
  ## --------------------------
  ## Block-level bookkeeping
  ## --------------------------
  
  strata <- levels(block)        # unique block labels
  S      <- length(strata)       # number of blocks
  
  # For each block we store:
  # - row indices (in original order)
  # - within-block ordering of those indices by y_scaled (descending)
  # - candidate m_s values = 1, ..., n_s - 1
  idx_list   <- vector("list", S)
  order_list <- vector("list", S)
  m_list     <- vector("list", S)
  
  for (s in seq_len(S)) {
    # Row indices in ORIGINAL order for block s
    idx_s <- which(block == strata[s])
    n_s   <- length(idx_s)
    
    if (n_s <= 1L) {
      stop("Block ", strata[s], " has size <= 1; cannot construct candidate u's.")
    }
    
    # Outcomes for this block, in original order
    y_s <- y_scaled[idx_s]
    
    # Within-block order (1..n_s) sorted by descending y_scaled
    ord_s <- order(y_s, decreasing = TRUE)
    
    idx_list[[s]]   <- idx_s
    order_list[[s]] <- ord_s
    m_list[[s]]     <- 1:(n_s - 1L)
  }
  
  ## --------------------------
  ## Enumerate all (m_1,...,m_S)
  ## --------------------------
  
  # Full grid of all combinations of m_s ∈ {1,...,n_s-1}
  m_grid <- do.call(
    what = expand.grid,
    args = c(m_list, KEEP.OUT.ATTRS = FALSE)
  )
  colnames(m_grid) <- strata
  
  K <- nrow(m_grid)  # number of candidate u vectors
  
  ## --------------------------
  ## Build all candidate u vectors
  ## --------------------------
  
  # Initialize n x K matrix; column k will be candidate u_k
  U <- matrix(0L, nrow = n, ncol = K)
  
  # Loop over each combination (m_1,...,m_S)
  for (k in seq_len(K)) {
    # Current m-vector for this candidate
    m_vec <- as.numeric(m_grid[k, ])
    
    # Start with all zeros for u
    u_vec <- integer(n)
    
    # Fill u block by block
    for (s in seq_len(S)) {
      idx_s <- idx_list[[s]]    # original row indices in block s
      ord_s <- order_list[[s]]  # within-block order by y_scaled
      m_s   <- m_vec[s]         # number of 1s to assign in block s
      
      if (m_s > 0L) {
        # Take top m_s positions in this block (by y_scaled)
        top_rows_in_block <- idx_s[ord_s[seq_len(m_s)]]
        # Set u = 1 for those rows in ORIGINAL order
        u_vec[top_rows_in_block] <- 1L
      }
    }
    
    # Store candidate u_k as column k
    U[, k] <- u_vec
  }
  
  ## --------------------------
  ## Return result
  ## --------------------------
  
  list(
    U      = U,      # n x K matrix of candidate u vectors
    m_grid = m_grid  # (m_1,...,m_S) associated with each column
  )
}