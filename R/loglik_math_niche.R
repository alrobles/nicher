#' Log-likelihood for the unweighted Gaussian niche model (math-scale parameterization)
#'
#' Optimizer-facing wrapper that accepts parameters in mathematical scale
#' (log-scale Cholesky diagonals) and delegates to \code{loglik_unweighted_math}.
#'
#' @param param_vector Named numeric vector of parameters in mathematical scale.
#'   Must contain \code{mu1..mup} and \code{L1..Lk} (k = p*(p+1)/2) in
#'   column-major lower-triangular order. Diagonal entries of L are on the log scale.
#' @param sam1 Data frame of presence points.
#' @param sam2 Data frame of background points.
#' @param mask Optional named numeric vector of fixed parameter values. Names must
#'   match those in \code{param_vector}. Masked parameters are held constant.
#' @param negative Logical; if \code{TRUE} (default) the negative log-likelihood is returned.
#' @return (Negative) log-likelihood value (scalar).
#' @export
#'
#' @examples
#' par <- get_ellip_par(spOccPnts)
#' L <- t(chol(par$S))
#' p <- length(par$mu)
#' k <- p * (p + 1) / 2
#' param_vector <- c(par$mu, log(diag(L)), L[lower.tri(L)])
#' names(param_vector) <- c(paste0("mu", 1:p), paste0("L", 1:k))
#' loglik_math_niche(param_vector, spOccPnts, samMPts)
loglik_math_niche <- function(param_vector, sam1, sam2, mask = NULL, negative = TRUE) {
  mu_names <- grep("^mu", names(param_vector), value = TRUE)
  if (length(mu_names) == 0) stop("No mu parameters found in param_vector")
  p <- length(mu_names)

  k <- p * (p + 1) / 2
  full_names <- c(paste0("mu", 1:p), paste0("L", 1:k))

  if (!is.null(mask)) {
    full <- stats::setNames(rep(NA_real_, length(full_names)), full_names)
    if (any(!names(mask) %in% full_names)) {
      stop("Invalid names in mask: ",
           paste(setdiff(names(mask), full_names), collapse = ", "),
           ". Expected names are: ", paste(full_names, collapse = ", "))
    }
    full[names(mask)] <- mask
    if (any(!names(param_vector) %in% full_names)) {
      stop("Invalid names in param_vector: ",
           paste(setdiff(names(param_vector), full_names), collapse = ", "),
           ". Expected names are: ", paste(full_names, collapse = ", "))
    }
    overlap <- intersect(names(param_vector), names(mask))
    if (length(overlap) > 0) {
      stop("param_vector and mask names must not overlap: ",
           paste(overlap, collapse = ", "))
    }
    full[names(param_vector)] <- param_vector
    if (anyNA(full)) {
      stop("Missing parameters: ", paste(names(full)[is.na(full)], collapse = ", "))
    }
    param_vector <- full
  }

  if (!all(full_names %in% names(param_vector))) {
    stop("param_vector does not contain all required parameter names")
  }

  bio <- math_to_bio_niche(param_vector)

  logL <- loglik_unweighted_math(sam1, sam2, bio$mu, bio$L)

  if (negative) return(-logL) else return(logL)
}
