#' Fit an ecological niche model
#'
#' High-level wrapper that fits a Gaussian ellipsoidal niche model to presence
#' and background environmental data.  Currently only Maximum Likelihood
#' Estimation (\code{"mle"}) is supported; a Bayesian backend (via Stan /
#' cmdstanr) is planned for a future release.
#'
#' @param env_occ Data frame or matrix of environmental variables at presence
#'   locations (rows = observations, columns = variables).
#' @param env_m   Data frame or matrix of environmental variables at background
#'   (M-hypothesis) locations.
#' @param method  Character.  Estimation method.  Currently only \code{"mle"}
#'   is accepted.
#' @param model   Character.  Niche model type.  One of:
#'   \describe{
#'     \item{\code{"unweighted"}}{Full Gaussian model fitted to both presence
#'       and background (Cholesky-parameterised MLE).}
#'     \item{\code{"weighted"}}{Presence-only model where the background sample
#'       is weighted by a Gaussian KDE of environmental space.}
#'     \item{\code{"presence_only"}}{Presence-only semi-log-likelihood (equal
#'       background weights).}
#'   }
#' @param ... Additional arguments passed to \code{\link[optimx]{optimx}}, for
#'   example \code{itnmax = 200} to increase the iteration budget.
#'
#' @return A \code{\link{nicher_object}} of class \code{"nicher"} with
#'   components:
#'   \describe{
#'     \item{\code{loglik}}{Log-likelihood at the optimum.}
#'     \item{\code{math_params}}{Named numeric vector of parameters in
#'       mathematical scale (\code{mu1..mup}, \code{L1..Lk}).}
#'     \item{\code{bio_params}}{List with biological-scale parameters:
#'       \code{mu}, \code{S}, \code{sigma}, \code{R}.}
#'     \item{\code{optim_table}}{Data frame with all optimizer runs and their
#'       results.}
#'     \item{\code{model}}{Model type used.}
#'     \item{\code{method}}{Estimation method used.}
#'   }
#'
#' @seealso \code{\link{nicher_object}}, \code{\link{print.nicher}},
#'   \code{\link{summary.nicher}}
#'
#' @export
#'
#' @examples
#' fit <- nicher(spOccPnts, samMPts, model = "presence_only")
#' print(fit)
#'
#' fit_unw <- nicher(spOccPnts, samMPts, model = "unweighted")
#' summary(fit_unw)
nicher <- function(env_occ, env_m,
                   method = "mle",
                   model  = c("unweighted", "weighted", "presence_only"),
                   ...) {
  method <- match.arg(method, choices = "mle")
  model  <- match.arg(model)

  if (method == "mle") {
    .nicher_mle(env_occ, env_m, model = model, ...)
  } else {
    stop("Method '", method, "' is not yet implemented.")
  }
}

# ---------------------------------------------------------------------------
# Internal: MLE dispatcher
# ---------------------------------------------------------------------------

.nicher_mle <- function(env_occ, env_m, model = "unweighted", ...) {
  env_occ <- as.matrix(env_occ)
  env_m   <- as.matrix(env_m)

  # Starting values: sample moments from presence data
  par0 <- get_ellip_par(env_occ)
  p    <- length(par0$mu)
  k    <- p * (p + 1L) / 2L

  # Lower-triangular Cholesky factor  L  s.t.  S = L %*% t(L)
  L0           <- t(chol(par0$S))
  L0_log       <- L0
  diag(L0_log) <- log(diag(L0))   # log-scale diagonal (math parameterisation)

  # Canonical math-scale parameter names
  mu_names <- paste0("mu", seq_len(p))
  L_names  <- paste0("L",  seq_len(k))

  start        <- c(par0$mu, L0_log[lower.tri(L0_log, diag = TRUE)])
  names(start) <- c(mu_names, L_names)

  # Build model-specific objective function (minimises negative log-likelihood)
  obj_fn <- .nicher_obj_fn(model, env_occ, env_m, mu_names, L_names)

  # Run multi-method optimisation
  opt_result <- suppressWarnings(
    optimx::optimx(
      par    = start,
      fn     = obj_fn,
      method = c("Nelder-Mead", "BFGS", "L-BFGS-B"),
      ...
    )
  )

  # Select best (lowest NLL)
  best_idx <- which.min(opt_result$value)
  best_row <- opt_result[best_idx, , drop = FALSE]

  # Extract optimal math-scale parameters
  n_par        <- length(start)
  math_vec     <- as.numeric(best_row[1L, seq_len(n_par)])
  names(math_vec) <- names(start)

  # Convert to biological scale
  bio <- .nicher_math_to_bio(math_vec, p, k)

  # Log-likelihood (not negated)
  loglik <- as.numeric(-best_row$value)

  nicher_object(
    loglik      = loglik,
    math_params = math_vec,
    bio_params  = bio,
    optim_table = as.data.frame(opt_result),
    model       = model,
    method      = "mle"
  )
}

# ---------------------------------------------------------------------------
# Internal: build model-specific objective (returns NLL to be minimised)
# ---------------------------------------------------------------------------

.nicher_obj_fn <- function(model, env_occ, env_m, mu_names, L_names) {
  all_names <- c(mu_names, L_names)

  if (model == "unweighted") {
    function(pv) {
      names(pv) <- all_names
      bio <- math_to_bio_niche(pv)
      -loglik_unweighted_math(env_occ, env_m, bio$mu, bio$L)
    }
  } else if (model == "weighted") {
    # Pre-compute KDE weights once; the background sample is fixed during optimisation
    log_weights <- kde_eval_cached(env_m, env_m)
    weights     <- exp(log_weights)
    function(pv) {
      names(pv) <- all_names
      bio <- math_to_bio_niche(pv)
      S   <- bio$L %*% t(bio$L)
      loglik_weighted_math(env_occ, env_m, bio$mu, S, weights = weights)
    }
  } else {   # "presence_only"
    function(pv) {
      names(pv) <- all_names
      bio <- math_to_bio_niche(pv)
      S   <- bio$L %*% t(bio$L)
      loglik_presenceonly_math(env_occ, env_m, bio$mu, S)
    }
  }
}

# ---------------------------------------------------------------------------
# Internal: math-scale -> full biological-scale parameter list
# ---------------------------------------------------------------------------

.nicher_math_to_bio <- function(math_vec, p, k) {
  bio <- math_to_bio_niche(math_vec)
  mu  <- bio$mu
  L   <- bio$L
  S   <- L %*% t(L)

  # Standard deviations and correlation matrix
  sigma <- sqrt(diag(S))
  D_inv <- diag(1 / sigma, nrow = p)
  R     <- D_inv %*% S %*% D_inv
  R     <- (R + t(R)) / 2   # enforce symmetry numerically

  # Propagate variable names from mu
  var_names        <- names(mu)
  dimnames(S)      <- list(var_names, var_names)
  dimnames(R)      <- list(var_names, var_names)
  names(sigma)     <- var_names

  list(mu = mu, S = S, sigma = sigma, R = R)
}
