#' Generate multiple starting points for niche model optimization
#'
#' Creates a set of starting parameter vectors on the math scale for use with
#' [loglik_niche_math()]. Ranges for each parameter are derived from the
#' environmental data of the accessibility area M using [get_range_df_niche()].
#'
#' @param env_data Data frame of environmental values from the accessibility
#'   area M. Used to define plausible ranges for the parameters.
#' @param num_starts Integer. Number of starting points to generate.
#' @param quant_vec Numeric vector of length 3 giving the quantiles used to
#'   define the lower, central and upper bounds for the `mu` parameters.
#'   Default `c(0.1, 0.5, 0.9)`.
#' @param method Character. Either `"sobol"` (requires the \pkg{pomp} package)
#'   or `"uniform"` for simple random sampling.
#'
#' @return A data frame with `num_starts` rows and one column per parameter.
#'   Column names follow the math‑scale convention: `mu1, mu2, ...`,
#'   `log_sigma1, log_sigma2, ...`, `v1, v2, ...` (the latter only if
#'   `p > 1`). Each row is a valid starting point.
#'
#' @export
#' @examples
#' \dontrun{
#' start_theta_multiple(example_env_m_2d, num_starts = 10, method = "uniform")
#' }
start_theta_multiple <- function(env_data, num_starts = 100,
                                 quant_vec = c(0.1, 0.5, 0.9),
                                 method = "sobol") {
  # Obtain parameter ranges
  ranges <- get_range_df_niche(env_data, quant_vec)
  lower <- ranges$lower
  upper <- ranges$upper
  names(lower) <- names(upper) <- rownames(ranges)
  
  if (method == "sobol") {
    if (!requireNamespace("pomp", quietly = TRUE)) {
      stop("Package 'pomp' is required for Sobol design. ",
           "Please install it or use method = 'uniform'.")
    }
    # pomp::sobol_design returns a matrix with rows as points
    start_mat <- pomp::sobol_design(lower = lower, upper = upper, nseq = num_starts)
    as.data.frame(start_mat)
  } else if (method == "uniform") {
    n_par <- length(lower)
    # Uniform random numbers in [0, 1]
    mat <- matrix(stats::runif(num_starts * n_par), nrow = num_starts, ncol = n_par)
    # Scale to [lower, upper] for each column
    for (j in seq_len(n_par)) {
      mat[, j] <- lower[j] + mat[, j] * (upper[j] - lower[j])
    }
    df <- as.data.frame(mat)
    names(df) <- names(lower)
    df
  } else {
    stop("method must be either 'sobol' or 'uniform'")
  }
}