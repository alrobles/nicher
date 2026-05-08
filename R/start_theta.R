#' Starting values for niche model on math scale
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @param skew Logical. If \code{TRUE}, append \code{p} zeros for the
#'   skew parameters (\code{alpha_1, ..., alpha_p}). The Gaussian-only
#'   start (\code{skew = FALSE}, the default) is used for
#'   \code{"presence_only"}, \code{"weighted"}, and
#'   \code{"ip_weighted"}; the skew start is used for
#'   \code{"skew_normal"} and \code{"skew_normal_weighted"}.
#' @param skew_t Logical. If \code{TRUE}, append \code{log(10)} for the
#'   skew-t degrees-of-freedom parameter \code{log_r} after the alpha
#'   block (used for \code{"skew_t"} and \code{"skew_t_weighted"}).
#'   Implies \code{skew = TRUE}.
#' @return Numeric vector of starting values for `theta`.
#' @export
#' @examples
#' start_theta(example_env_occ_2d)
#' start_theta(example_env_occ_2d, skew = TRUE)
#' start_theta(example_env_occ_2d, skew_t = TRUE)
start_theta <- function(env_occ, skew = FALSE, skew_t = FALSE) {
  if (isTRUE(skew_t)) skew <- TRUE
  p <- ncol(env_occ)
  mu0 <- colMeans(env_occ, na.rm = TRUE)
  sigma0 <- apply(env_occ, 2, stats::sd, na.rm = TRUE)
  sigma0 <- pmax(sigma0, 1e-6) # avoid zero
  log_sigma0 <- log(sigma0)
  v0 <- rep(0, p * (p - 1) / 2)
  alpha0 <- if (isTRUE(skew)) rep(0, p) else numeric(0)
  log_r0 <- if (isTRUE(skew_t)) log(10) else numeric(0)
  c(mu0, log_sigma0, v0, alpha0, log_r0)
}
