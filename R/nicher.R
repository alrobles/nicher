#' Build a nicher object from optimization results
#'
#' Internal constructor that wraps optimization output into an S3 object of
#' class \code{"nicher"}.  Users typically obtain a \code{nicher} object via
#' the high-level \code{\link{nicher}} function rather than calling this
#' constructor directly.
#'
#' @param loglik Scalar numeric; the maximized log-likelihood value.
#' @param math_params Named numeric vector of optimized parameters in
#'   mathematical (log-Cholesky) scale.
#' @param bioscale_params List with biological-scale parameters:
#'   \describe{
#'     \item{mu}{Numeric vector of optimized environmental optima.}
#'     \item{L}{Lower-triangular Cholesky factor matrix.}
#'     \item{S}{Covariance matrix (S = L \%*\% t(L)).}
#'     \item{R}{Correlation matrix derived from S.}
#'     \item{variances}{Numeric vector of marginal variances (diagonal of S).}
#'   }
#' @param optimx_table Data frame of optimization results (all starting
#'   conditions and methods tried).
#' @param model Character string; one of \code{"unweighted"},
#'   \code{"weighted"}, or \code{"presence_only"}.
#' @param method Character string; optimization method used (e.g.,
#'   \code{"mle"}).
#' @param env_names Optional character vector of variable names for the
#'   environmental axes.
#' @return An object of class \code{"nicher"}.
#' @keywords internal
new_nicher <- function(loglik, math_params, bioscale_params,
                       optimx_table, model, method, env_names = NULL) {
  obj <- structure(
    list(
      loglik         = loglik,
      math_params    = math_params,
      bioscale_params = bioscale_params,
      optimx_table   = optimx_table,
      model          = model,
      method         = method,
      env_names      = env_names
    ),
    class = "nicher"
  )
  check_nicher(obj)
  obj
}

#' Fit an ecological niche model
#'
#' High-level wrapper that fits a Gaussian ellipsoidal niche model to presence
#' and background environmental data using maximum likelihood estimation.
#' Returns a \code{nicher} object that encapsulates the optimized parameters in
#' both mathematical and biological scale, the log-likelihood, and the full
#' optimization table.
#'
#' @param env_occ Data frame or matrix of environmental values at presence
#'   points (rows = observations, columns = environmental variables).
#' @param env_m Data frame or matrix of environmental values at background
#'   (M-region) points.
#' @param method Optimization framework.  Currently only \code{"mle"} is
#'   supported; a Bayesian backend is planned for a future release.
#' @param model Niche model variant.  One of:
#'   \describe{
#'     \item{\code{"unweighted"}}{Full Gaussian model with uniform background
#'       weights; maximises a log-likelihood difference between presence and
#'       background.}
#'     \item{\code{"weighted"}}{Full Gaussian model with KDE-derived background
#'       weights; corrects for environmental sampling bias.}
#'     \item{\code{"presence_only"}}{Semi-log (Poisson process) model that
#'       requires only presence and background points without assuming a known
#'       prevalence.}
#'   }
#' @param itnmax Integer; maximum number of optimizer iterations (passed to
#'   \code{\link[optimx]{optimx}}).
#' @param ... Additional arguments (currently unused; reserved for future
#'   extensions).
#'
#' @return An object of class \code{"nicher"} containing:
#'   \describe{
#'     \item{\code{loglik}}{Maximized log-likelihood (scalar numeric).}
#'     \item{\code{math_params}}{Named numeric vector of optimized parameters
#'       in mathematical (log-Cholesky) scale.}
#'     \item{\code{bioscale_params}}{List with biological-scale parameters:
#'       \code{mu} (optima), \code{L} (Cholesky factor), \code{S} (covariance
#'       matrix), \code{R} (correlation matrix), and \code{variances} (marginal
#'       variances).}
#'     \item{\code{optimx_table}}{Data frame of all optimization runs.}
#'     \item{\code{model}}{Model variant used.}
#'     \item{\code{method}}{Optimization method used.}
#'     \item{\code{env_names}}{Environmental variable names.}
#'   }
#'
#' @export
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' print(fit)
#' summary(fit)
nicher <- function(env_occ, env_m,
                   method = "mle",
                   model  = c("presence_only", "unweighted", "weighted"),
                   itnmax = 200,
                   ...) {
  method <- match.arg(method, choices = "mle")
  model  <- match.arg(model)

  env_occ <- as.data.frame(env_occ)
  env_m   <- as.data.frame(env_m)

  env_names <- colnames(env_occ)
  p         <- ncol(env_occ)

  if (ncol(env_m) != p) {
    stop("'env_occ' and 'env_m' must have the same number of columns")
  }

  if (method == "mle") {
    result <- .nicher_mle(env_occ, env_m, model, itnmax, env_names)
  } else {
    stop("method '", method, "' is not yet implemented")
  }

  result
}

# ---------------------------------------------------------------------------
# Internal MLE fitting dispatcher
# ---------------------------------------------------------------------------

.nicher_mle <- function(env_occ, env_m, model, itnmax, env_names) {
  par_init <- get_ellip_par(env_occ)
  p        <- length(par_init$mu)
  k        <- p * (p + 1) / 2

  L_init       <- t(chol(par_init$S))
  L_log        <- L_init
  diag(L_log)  <- log(diag(L_init))
  param_init   <- c(par_init$mu, L_log[lower.tri(L_log, diag = TRUE)])
  param_names  <- c(paste0("mu", seq_len(p)), paste0("L", seq_len(k)))
  names(param_init) <- param_names

  like_fn <- .make_like_fn(model, env_occ, env_m, p, param_names)

  methods_try <- c("Nelder-Mead", "BFGS", "L-BFGS-B")
  suppressWarnings({
    opt_table <- optimx::optimx(
      par    = param_init,
      fn     = like_fn,
      method = methods_try,
      itnmax = itnmax
    )
  })

  best_idx   <- which.min(opt_table$value)
  best_row   <- opt_table[best_idx, , drop = FALSE]
  best_vec   <- as.numeric(best_row[, seq_len(length(param_init))])
  names(best_vec) <- param_names

  loglik <- -best_row$value

  bio   <- math_to_bio_niche(best_vec)
  S_bio <- bio$L %*% t(bio$L)
  D_inv <- diag(1 / sqrt(diag(S_bio)), nrow = nrow(S_bio))
  R_bio <- D_inv %*% S_bio %*% D_inv

  bioscale_params <- list(
    mu        = bio$mu,
    L         = bio$L,
    S         = S_bio,
    R         = R_bio,
    variances = diag(S_bio)
  )

  if (!is.null(env_names)) {
    names(bioscale_params$mu)        <- env_names
    names(bioscale_params$variances) <- env_names
    rownames(bioscale_params$L)      <- env_names
    colnames(bioscale_params$L)      <- env_names
    rownames(bioscale_params$S)      <- env_names
    colnames(bioscale_params$S)      <- env_names
    rownames(bioscale_params$R)      <- env_names
    colnames(bioscale_params$R)      <- env_names
  }

  new_nicher(
    loglik          = loglik,
    math_params     = best_vec,
    bioscale_params = bioscale_params,
    optimx_table    = opt_table,
    model           = model,
    method          = "mle",
    env_names       = env_names
  )
}

# ---------------------------------------------------------------------------
# Internal helper: build the objective function for a given model
# ---------------------------------------------------------------------------

.make_like_fn <- function(model, env_occ, env_m, p, param_names) {
  k <- p * (p + 1) / 2

  if (model == "unweighted") {
    function(param_vec) {
      names(param_vec) <- param_names
      bio <- math_to_bio_niche(param_vec)
      -loglik_unweighted_math(env_occ, env_m, bio$mu, bio$L)
    }
  } else if (model == "presence_only") {
    function(param_vec) {
      names(param_vec) <- param_names
      bio <- math_to_bio_niche(param_vec)
      S   <- bio$L %*% t(bio$L)
      loglik_presenceonly_math(env_occ, env_m, bio$mu, S)
    }
  } else if (model == "weighted") {
    weights <- exp(kde_eval_cached(env_m, env_m))
    weights <- weights / sum(weights)
    function(param_vec) {
      names(param_vec) <- param_names
      bio <- math_to_bio_niche(param_vec)
      S   <- bio$L %*% t(bio$L)
      loglik_weighted_math(env_occ, env_m, bio$mu, S, weights = weights)
    }
  } else {
    stop("Unknown model: ", model)
  }
}
