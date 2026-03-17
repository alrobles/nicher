#' Construct a nicher S3 object
#'
#' Internal constructor that wraps MLE (or future Bayesian) optimization results
#' into a classed \code{nicher} object.  End users should call \code{\link{nicher}}
#' rather than this function directly.
#'
#' @param loglik    Numeric scalar. Log-likelihood at the optimum (positive scale,
#'   i.e., \emph{not} negated).
#' @param math_params Named numeric vector of optimized parameters in
#'   mathematical scale (\code{mu1..mup}, \code{L1..Lk} with log-scale
#'   diagonal entries).
#' @param bio_params List of biological-scale parameters:
#'   \describe{
#'     \item{mu}{Numeric vector of means (niche optimum).}
#'     \item{S}{Covariance matrix (p x p).}
#'     \item{sigma}{Vector of standard deviations (sqrt of S diagonal).}
#'     \item{R}{Correlation matrix derived from S.}
#'   }
#' @param optim_table Data frame of all optimization results from
#'   \code{\link[optimx]{optimx}}.
#' @param model Character. One of \code{"unweighted"}, \code{"weighted"},
#'   or \code{"presence_only"}.
#' @param method Character. Estimation method used, e.g. \code{"mle"}.
#'
#' @return An object of class \code{"nicher"}.
#'
#' @seealso \code{\link{nicher}}, \code{\link{print.nicher}},
#'   \code{\link{summary.nicher}}
#'
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' L   <- t(chol(par$S))
#' p   <- length(par$mu)
#' k   <- p * (p + 1L) / 2L
#' L_log        <- L
#' diag(L_log)  <- log(diag(L))
#' math_params  <- c(par$mu, L_log[lower.tri(L_log, diag = TRUE)])
#' names(math_params) <- c(paste0("mu", seq_len(p)), paste0("L", seq_len(k)))
#' S     <- L %*% t(L)
#' sigma <- sqrt(diag(S))
#' D_inv <- diag(1 / sigma, nrow = p)
#' R     <- (D_inv %*% S %*% D_inv + t(D_inv %*% S %*% D_inv)) / 2
#' bio_params <- list(mu = par$mu, S = S, sigma = sigma, R = R)
#' ll <- loglik_unweighted_math(spOccPnts, samMPts, par$mu, L)
#' obj <- nicher_object(ll, math_params, bio_params, data.frame(), "unweighted", "mle")
#' class(obj)
nicher_object <- function(loglik, math_params, bio_params, optim_table,
                          model, method) {
  structure(
    list(
      loglik      = loglik,
      math_params = math_params,
      bio_params  = bio_params,
      optim_table = optim_table,
      model       = model,
      method      = method
    ),
    class = "nicher"
  )
}

#' Test if an object is of class \code{nicher}
#'
#' @param x Any R object.
#' @return Logical scalar.
#' @export
#'
#' @examples
#' is.nicher(42)
is.nicher <- function(x) inherits(x, "nicher")
