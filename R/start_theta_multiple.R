#' Generate multiple starting points for niche model optimization
#'
#' Creates a set of starting parameter vectors on the math scale for use with
#' the log-likelihood functions. Ensures all starting values are strictly
#' numeric and finite, even if env_data is supplied as a data frame.
#'
#' @param env_data Environmental data (matrix or data frame). Must be numeric.
#' @param num_starts Integer, number of starting points.
#' @param quant_vec Quantiles for mu ranges.
#' @param method "sobol" (Sobol design) or "uniform".
#'
#' @return A data frame of dimension num_starts × num_parameters.
#' @export
start_theta_multiple <- function(env_data, num_starts = 100,
                                 quant_vec = c(0.1, 0.5, 0.9),
                                 method = "sobol") {
  # ------------------------------------------------------------------
  # 1. COERCE TO NUMERIC MATRIX (THIS IS THE KEY FIX)
  # ------------------------------------------------------------------
  if (is.data.frame(env_data)) {
    env_data <- as.matrix(env_data)
  }

  if (!is.numeric(env_data)) {
    stop("env_data must be numeric. Convert your data before calling.")
  }

  storage.mode(env_data) <- "double"

  # ------------------------------------------------------------------
  # 2. Compute ranges (now always numeric)
  # ------------------------------------------------------------------
  ranges <- get_range_df_niche(env_data, quant_vec)

  lower <- as.numeric(ranges$lower)
  upper <- as.numeric(ranges$upper)
  param_names <- rownames(ranges)

  names(lower) <- param_names
  names(upper) <- param_names

  # ------------------------------------------------------------------
  # 3. Generate design
  # ------------------------------------------------------------------
  if (method == "sobol") {
    if (!requireNamespace("pomp", quietly = TRUE)) {
      stop("Package 'pomp' is required for method='sobol'.")
    }

    start_mat <- pomp::sobol_design(
      lower = lower,
      upper = upper,
      nseq  = num_starts
    )

    df <- as.data.frame(start_mat)
  } else if (method == "uniform") {
    n_par <- length(lower)
    mat <- matrix(stats::runif(num_starts * n_par),
      nrow = num_starts, ncol = n_par
    )

    for (j in seq_len(n_par)) {
      mat[, j] <- lower[j] + mat[, j] * (upper[j] - lower[j])
    }

    df <- as.data.frame(mat)
    names(df) <- param_names
  } else {
    stop("method must be either 'sobol' or 'uniform'")
  }

  # ------------------------------------------------------------------
  # 4. Final guarantee: all numeric, no NA, no Inf
  # ------------------------------------------------------------------
  df[] <- lapply(df, function(col) {
    col <- as.numeric(col)
    if (any(!is.finite(col))) {
      stop("Non-finite values detected in generated start parameters.")
    }
    col
  })

  return(df)
}
