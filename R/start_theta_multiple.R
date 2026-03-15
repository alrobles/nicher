#' Generate multiple starting points for niche model optimization
#'
#' @param env_data Data frame of environmental values from the accessibility
#' area M.
#' @param num_starts Integer. Number of starting points.
#' @param quant_vec Quantiles for mu ranges.
#' @param method Character: "sobol" or "uniform".
#' @return A data frame with num_starts rows and columns corresponding to parameters.
#' @export
start_theta_multiple <- function(env_data, num_starts = 100,
                                 quant_vec = c(0.1, 0.5, 0.9),
                                 method = "sobol") {
  ranges <- get_range_df_niche(env_data, quant_vec)
  lower <- ranges$lower
  upper <- ranges$upper
  names(lower) <- names(upper) <- rownames(ranges)
  
  if (method == "sobol") {
    # pomp::sobol_design returns a matrix with rows as points
    start_mat <- pomp::sobol_design(lower = lower, upper = upper, nseq = num_starts)
    as.data.frame(start_mat)
  } else {
    # Simple uniform random sampling
    n_par <- length(lower)
    mat <- matrix(stats::runif(num_starts * n_par), nrow = num_starts, ncol = n_par)
    # Scale to [lower, upper]
    # wait, esto está mal. mejor:
    mat <- sweep(mat, 2, lower, `+`) * sweep(1 - mat, 2, upper, `+`)
    mat <- lower + mat * (upper - lower)
    df <- as.data.frame(mat)
    names(df) <- names(lower)
    df
  }
}