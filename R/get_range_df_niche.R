#' Generate reasonable ranges for starting parameters (math scale)
#'
#' @param env_data Data frame of environmental values.
#' @param quant_vec Numeric vector of length 3: lower, middle, upper quantiles.
#' @return A data frame with rows corresponding to parameter names and columns
#' lower, center, upper.
#' @keywords internal
get_range_df_niche <- function(env_data, quant_vec = c(0.1, 0.5, 0.9)) {
  p <- ncol(env_data)
  n_par <- 2 * p + p * (p - 1) / 2
  # Create names in order: mu1..p, log_sigma1..p, v1..q
  mu_names <- paste0("mu", 1:p)
  log_sigma_names <- paste0("log_sigma", 1:p)
  s_par <- if (p > 1) paste0("s_par", 1:(p*(p-1)/2)) else character(0)
  all_names <- c(mu_names, log_sigma_names, s_par)
  
  # Initialize data frame
  ranges <- data.frame(
    lower = numeric(n_par),
    center = numeric(n_par),
    upper = numeric(n_par),
    row.names = all_names
  )
  
  # mu: based on quantiles of observed values
  for (i in 1:p) {
    vals <- env_data[, i]
    qq <- stats::quantile(vals, probs = quant_vec, na.rm = TRUE)
    ranges[mu_names[i], ] <- as.numeric(qq)
  }
  
  # log_sigma: based on log of sd, with expansion factor
  sigma_obs <- apply(env_data, 2, stats::sd, na.rm = TRUE)
  sigma_obs <- pmax(sigma_obs, 1e-6)
  log_sigma_obs <- log(sigma_obs)
  # Expand by ±1 on log scale (factor of e ≈ 2.718)
  for (i in 1:p) {
    center <- log_sigma_obs[i]
    ranges[log_sigma_names[i], ] <- c(center - 1, center, center + 1)
  }
  
  # v: symmetric around 0, with reasonable bounds (e.g., -3 to 3)
  if (length(s_par) > 0) {
    for (nm in s_par) {
      ranges[nm, ] <- c(-3, 0, 3)
    }
  }
  
  ranges
}