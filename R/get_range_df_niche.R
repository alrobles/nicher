#' Generate reasonable ranges for starting parameters (math scale)
#'
#' Ensures that env_data is strictly numeric before computing ranges.
#' This prevents downstream issues in Sobol designs and multi-start
#' optimization (e.g., character / factor leakage).
#'
#' @param env_data Numeric matrix or data frame of environmental values.
#' @param quant_vec Numeric vector of quantiles (length 3).
#'
#' @return A data frame with rows = parameter names and columns
#'         lower, center, upper.
#' @keywords internal
get_range_df_niche <- function(env_data, quant_vec = c(0.1, 0.5, 0.9)) {

  # --------------------------------------------------------------------
  # (1) COERCE env_data TO STRICTLY NUMERIC MATRIX (CRITICAL PATCH)
  # --------------------------------------------------------------------
  if (is.data.frame(env_data))
    env_data <- as.matrix(env_data)

  if (!is.numeric(env_data))
    stop("env_data must be numeric; convert before calling get_range_df_niche().")

  storage.mode(env_data) <- "double"

  # --------------------------------------------------------------------
  # (2) Compute basic dimensions and names
  # --------------------------------------------------------------------
  p <- ncol(env_data)
  n_par <- 2 * p + p * (p - 1) / 2

  mu_names        <- paste0("mu",        seq_len(p))
  log_sigma_names <- paste0("log_sigma", seq_len(p))
  s_par_names     <- if (p > 1) paste0("s_par", seq_len(p * (p - 1) / 2)) else character(0)

  all_names <- c(mu_names, log_sigma_names, s_par_names)

  ranges <- data.frame(
    lower  = numeric(n_par),
    center = numeric(n_par),
    upper  = numeric(n_par),
    row.names = all_names
  )

  # --------------------------------------------------------------------
  # (3) mu ranges based on quantiles
  # --------------------------------------------------------------------
  for (i in seq_len(p)) {
    vals <- env_data[, i]
    qq   <- stats::quantile(vals, probs = quant_vec, na.rm = TRUE)
    qq   <- as.numeric(qq)  # ensure numeric, drop attributes
    ranges[mu_names[i], ] <- qq
  }

  # --------------------------------------------------------------------
  # (4) log_sigma: log(sd) ± 1 range
  # --------------------------------------------------------------------
  sigma_obs     <- apply(env_data, 2, stats::sd, na.rm = TRUE)
  sigma_obs     <- pmax(sigma_obs, 1e-6)  # avoid zeros
  log_sigma_obs <- log(sigma_obs)

  for (i in seq_len(p)) {
    center <- log_sigma_obs[i]
    ranges[log_sigma_names[i], ] <- c(center - 1, center, center + 1)
  }

  # --------------------------------------------------------------------
  # (5) s_par ranges: [-3, 0, 3]
  # --------------------------------------------------------------------
  if (length(s_par_names) > 0) {
    for (nm in s_par_names) {
      ranges[nm, ] <- c(-3, 0, 3)
    }
  }

  # --------------------------------------------------------------------
  # (6) Final global sanity check
  # --------------------------------------------------------------------
  if (any(!is.finite(as.matrix(ranges)))) {
    stop("Non-finite values generated in get_range_df_niche(). Check env_data.")
  }

  return(ranges)
}
