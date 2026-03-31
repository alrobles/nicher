#' Starting values for niche model on math scale
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @return Numeric vector of starting values for `theta`.
#' @export
#' @examples
#' start_theta(example_env_occ_2d)
start_theta <- function(env_occ) {
  p <- ncol(env_occ)
  mu0 <- colMeans(env_occ, na.rm = TRUE)
  sigma0 <- apply(env_occ, 2, stats::sd, na.rm = TRUE)
  sigma0 <- pmax(sigma0, 1e-6) # avoid zero
  log_sigma0 <- log(sigma0)
  v0 <- rep(0, p * (p - 1) / 2)
  c(mu0, log_sigma0, v0)
}
