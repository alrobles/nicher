#' Starting values for niche model on math scale
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @param skew Logical. If \code{TRUE}, append \code{p} zeros for the
#'   skew parameters (\code{alpha_1, ..., alpha_p}). The Gaussian-only
#'   start (\code{skew = FALSE}, the default) is used for
#'   \code{"presence_only"}, \code{"weighted"}, and
#'   \code{"kde_bias_corrected"}; the skew start is used for
#'   \code{"skew_normal"} and \code{"skew_normal_weighted"}.
#' @return Numeric vector of starting values for `theta`.
#' @export
#' @examples
#' start_theta(example_env_occ_2d)
#' start_theta(example_env_occ_2d, skew = TRUE)
start_theta <- function(env_occ, skew = FALSE) {
  p <- ncol(env_occ)
  mu0 <- colMeans(env_occ, na.rm = TRUE)
  sigma0 <- apply(env_occ, 2, stats::sd, na.rm = TRUE)
  sigma0 <- pmax(sigma0, 1e-6) # avoid zero
  log_sigma0 <- log(sigma0)
  v0 <- rep(0, p * (p - 1) / 2)
  alpha0 <- if (isTRUE(skew)) rep(0, p) else numeric(0)
  c(mu0, log_sigma0, v0, alpha0)
}
