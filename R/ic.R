# R/ic.R
# Information criteria (AIC, BIC) for nicher model objects.
#
# Exports: logLik.nicher, nobs.nicher, AIC.nicher, BIC.nicher,
#          compare_nicher.

# -----------------------------------------------------------------------
# Cheap deterministic fingerprint used to detect non-comparable inputs
# across multiple `nicher` fits passed to `compare_nicher()`.
#
# Returns a small list whose components are stable under row reordering
# of `x` only if the user reorders columns identically -- which is
# correct, since reordering columns changes the model.
# -----------------------------------------------------------------------
.env_fingerprint <- function(x) {
  if (is.null(x)) return(NULL)
  m <- as.matrix(x)
  list(
    nrow      = nrow(m),
    ncol      = ncol(m),
    col_names = colnames(m),
    col_mean  = if (nrow(m) > 0L) unname(colMeans(m)) else numeric(0L),
    col_sd    = if (nrow(m) > 1L)
                  unname(apply(m, 2L, stats::sd)) else
                  rep(NA_real_, ncol(m))
  )
}

# -----------------------------------------------------------------------
# Internal: evaluate the un-penalised log-likelihood at a converged
# theta, used by `optimize_niche()` to populate
# `best$loglik_unpenalised`. All seven likelihood families dispatch
# here; the wrappers' `*_lambda = 0, *_center = NULL` flags strip the
# ridge so the returned value is the bare data log-likelihood.
# -----------------------------------------------------------------------
.compute_unpenalised_loglik <- function(theta, likelihood, env_occ, env_m,
                                        weighted_inputs, eta) {
  occ_mat <- as.matrix(env_occ)
  p <- ncol(occ_mat)
  zero_p <- numeric(p)

  M_den <- if (!is.null(env_m) && !is.null(weighted_inputs)) {
    as.matrix(env_m)[weighted_inputs$den_idx, , drop = FALSE]
  } else {
    NULL
  }

  neg <- switch(
    likelihood,
    presence_only = loglik_niche_math_presence_only(
      theta, env_occ = occ_mat, neg = TRUE, eta = eta
    ),
    kde_bias_corrected = loglik_niche_math_kde_bias_corrected(
      theta, env_occ = occ_mat, env_m = as.matrix(env_m), neg = TRUE,
      den_idx       = weighted_inputs$den_idx,
      kde_idx       = weighted_inputs$kde_idx,
      precomp_w_den = weighted_inputs$w_den,
      eta           = eta
    ),
    weighted = loglik_niche_math_weighted_cpp(
      theta = theta, env_occ = occ_mat, M_den = M_den,
      w_occ = weighted_inputs$w_occ, w_den = weighted_inputs$w_den,
      prior_log_sigma_center = zero_p, prior_log_sigma_lambda = 0,
      eta = eta,
      prior_mu_center    = NULL, prior_mu_lambda    = 0,
      prior_alpha_lambda = 0
    ),
    skew_normal = loglik_niche_math_skew_normal_cpp(
      theta = theta, env_occ = occ_mat, eta = eta,
      prior_mu_center        = NULL, prior_mu_lambda        = 0,
      prior_log_sigma_center = NULL, prior_log_sigma_lambda = 0,
      prior_alpha_lambda     = 0
    ),
    skew_normal_weighted = loglik_niche_math_skew_normal_weighted_cpp(
      theta = theta, env_occ = occ_mat, M_den = M_den,
      w_occ = weighted_inputs$w_occ, w_den = weighted_inputs$w_den,
      prior_log_sigma_center = zero_p, prior_log_sigma_lambda = 0,
      eta = eta,
      prior_mu_center    = NULL, prior_mu_lambda    = 0,
      prior_alpha_lambda = 0
    ),
    skew_t = loglik_niche_math_skew_t_cpp(
      theta = theta, env_occ = occ_mat, eta = eta,
      prior_mu_center        = NULL, prior_mu_lambda        = 0,
      prior_log_sigma_center = NULL, prior_log_sigma_lambda = 0,
      prior_alpha_lambda     = 0
    ),
    skew_t_weighted = loglik_niche_math_skew_t_weighted_cpp(
      theta = theta, env_occ = occ_mat, M_den = M_den,
      w_occ = weighted_inputs$w_occ, w_den = weighted_inputs$w_den,
      prior_log_sigma_center = zero_p, prior_log_sigma_lambda = 0,
      eta = eta,
      prior_mu_center    = NULL, prior_mu_lambda    = 0,
      prior_alpha_lambda = 0
    )
  )
  -as.numeric(neg)
}

# -----------------------------------------------------------------------
# logLik.nicher
# -----------------------------------------------------------------------

#' Log-likelihood of a fitted nicher model
#'
#' Returns the **un-penalised** log-likelihood of a \code{"nicher"} fit at
#' the converged \eqn{\theta}. This is the appropriate quantity for
#' likelihood-based model comparison via \code{\link{AIC}} or
#' \code{\link{BIC}}: even when \code{\link{optimize_niche}} minimises a
#' ridge-penalised objective \eqn{-\log L(\theta) + \mathrm{penalty}(\theta)},
#' the comparison metric is the bare \eqn{\log L(\theta)} evaluated at
#' the optimised \eqn{\theta}.
#'
#' For penalised fits this means the value returned is **not** the same
#' as \code{x$best$loglik}: \code{best$loglik} is the negative of the
#' penalised objective at the optimum (and is monotone in fit quality
#' but not directly comparable across penalty strengths), while
#' \code{logLik(x)} is the bare data log-likelihood. \code{logLik(x)}
#' is what AIC/BIC require.
#'
#' Stored on \code{x$best$loglik_unpenalised} by
#' \code{\link{optimize_niche}}; older fits without that field are
#' detected and an informative error is raised.
#'
#' @section Effective degrees of freedom:
#'
#' \code{df} is reported naively as \code{length(theta)}: 2p + p(p-1)/2
#' for the Gaussian families, 3p + p(p-1)/2 for the skew-normal families,
#' 3p + p(p-1)/2 + 1 for the skew-t families. This **ignores** the
#' shrinkage induced by the ridge penalties (\code{prior_mu_lambda},
#' \code{prior_log_sigma_lambda}, \code{prior_alpha_lambda}). For weak
#' penalties (\eqn{\lambda \le 1}) the bias is small; for stronger
#' penalties the effective df is lower than \code{length(theta)} so the
#' penalty (\eqn{2k} or \eqn{k \log n}) is conservative — IC values
#' will under-favour the more-regularised model. Quantifying effective
#' df via the trace of the influence matrix is left to a future release.
#'
#' @param object A \code{"nicher"} object returned by
#'   \code{\link{optimize_niche}}.
#' @param ... Ignored.
#'
#' @return An object of class \code{"logLik"} with attributes
#'   \code{df} (number of free parameters) and \code{nobs} (number of
#'   occurrence rows used in the fit).
#' @export
#' @examples
#' \dontrun{
#' fit <- optimize_niche(env_occ = example_env_occ_2d,
#'                       env_m   = example_env_m_2d,
#'                       num_starts = 20L, likelihood = "weighted")
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' }
logLik.nicher <- function(object, ...) {
  if (is.null(object$best$loglik_unpenalised)) {
    stop("This nicher object lacks `best$loglik_unpenalised` (likely fit ",
         "with nicher < 3.2.0). Refit with the current version to use ",
         "logLik() / AIC() / BIC().")
  }
  if (is.null(object$best$nobs)) {
    stop("This nicher object lacks `best$nobs`. Refit with the current ",
         "version of nicher.")
  }
  ll <- structure(
    as.numeric(object$best$loglik_unpenalised),
    df    = length(object$best$theta),
    nobs  = as.integer(object$best$nobs),
    class = "logLik"
  )
  ll
}

# -----------------------------------------------------------------------
# nobs.nicher
# -----------------------------------------------------------------------

#' Number of observations used to fit a nicher model
#'
#' Returns the number of presence rows (\code{nrow(env_occ)}) that were
#' passed to \code{\link{optimize_niche}}. This is the sample size used
#' by \code{\link{BIC}} and is the conventional choice for the
#' presence-only and weighted families: even though the weighted
#' likelihood includes a denominator over background points, the
#' observed quantities being modeled are the \eqn{n_{occ}} occurrences.
#'
#' @param object A \code{"nicher"} object.
#' @param ... Ignored.
#'
#' @return Integer scalar.
#' @export
nobs.nicher <- function(object, ...) {
  if (is.null(object$best$nobs)) {
    stop("This nicher object lacks `best$nobs`. Refit with nicher >= 3.2.0.")
  }
  as.integer(object$best$nobs)
}

# -----------------------------------------------------------------------
# AIC / BIC dispatch
# -----------------------------------------------------------------------

#' Akaike information criterion for a nicher fit
#'
#' Computes \eqn{-2 \log L + 2 k}, where \eqn{\log L} is the un-penalised
#' log-likelihood (see \code{\link{logLik.nicher}}) and \eqn{k} is the
#' number of free parameters.
#'
#' @inheritParams stats::AIC
#' @param ... Additional fitted model objects.
#' @param k Numeric, the penalty per parameter; default \code{2} (AIC).
#'   Pass \code{k = log(nobs(fit))} to recover BIC.
#'
#' @return Numeric scalar, or a \code{data.frame} when multiple objects
#'   are passed.
#' @export
AIC.nicher <- function(object, ..., k = 2) {
  fits     <- c(list(object), list(...))
  call_nms <- as.character(sys.call())[-1L]
  call_nms <- call_nms[seq_along(fits)]
  lls      <- lapply(fits, logLik)
  vals     <- vapply(lls, function(x) {
    -2 * as.numeric(x) + k * attr(x, "df")
  }, numeric(1L))
  if (length(fits) == 1L) return(vals[[1L]])
  data.frame(
    df  = vapply(lls, function(x) as.integer(attr(x, "df")), integer(1L)),
    AIC = vals,
    row.names = call_nms
  )
}

#' Bayesian information criterion for a nicher fit
#'
#' Computes \eqn{-2 \log L + k \log n}, where \eqn{n} is
#' \code{nobs(object)}. Calls \code{\link[stats]{BIC}}.
#'
#' @inheritParams stats::BIC
#' @param ... Additional fitted model objects.
#'
#' @return Numeric scalar, or a \code{data.frame} when multiple objects
#'   are passed.
#' @export
BIC.nicher <- function(object, ...) {
  fits     <- c(list(object), list(...))
  call_nms <- as.character(sys.call())[-1L]
  call_nms <- call_nms[seq_along(fits)]
  lls      <- lapply(fits, logLik)
  vals     <- vapply(lls, function(x) {
    -2 * as.numeric(x) + log(attr(x, "nobs")) * attr(x, "df")
  }, numeric(1L))
  if (length(fits) == 1L) return(vals[[1L]])
  data.frame(
    df  = vapply(lls, function(x) as.integer(attr(x, "df")), integer(1L)),
    BIC = vals,
    row.names = call_nms
  )
}

# -----------------------------------------------------------------------
# compare_nicher
# -----------------------------------------------------------------------

#' Compare nicher fits side-by-side via AIC and BIC
#'
#' Builds a one-row-per-model summary table for a list of \code{"nicher"}
#' fits, including log-likelihood, parameter count, AIC, BIC, deltas
#' relative to the best model, and Akaike weights.
#'
#' All fits MUST have been produced from the same \code{env_occ}; any
#' mismatch in \code{nobs} or in the variable-name ordering raises an
#' error, since IC values are not comparable across heterogeneous data.
#'
#' @param ... Two or more \code{"nicher"} objects, optionally named, OR a
#'   single named list of such objects.
#' @param sort_by Character. One of \code{"AIC"} (default), \code{"BIC"},
#'   or \code{"loglik"}: which metric to order rows by. Best model is
#'   always at the top.
#' @param comparison_basis Character, one of \code{"unpenalised"}
#'   (default) or \code{"penalised"}. Determines which log-likelihood is
#'   used as the basis for AIC and BIC.
#' \itemize{
#'   \item \code{"unpenalised"}: the bare data log-likelihood
#'         \eqn{\ell(\hat\theta)} returned by
#'         \code{\link{logLik.nicher}}. Standard frequentist convention,
#'         and the right basis when comparing **different likelihood
#'         families** at fixed penalty knobs.
#'   \item \code{"penalised"}: the optimiser's actual objective,
#'         \eqn{\ell(\hat\theta) - \mathrm{pen}(\hat\theta)}, i.e.
#'         \code{x$best$loglik}. The right basis when comparing **the
#'         same family at different penalty strengths**, since the
#'         unpenalised log-likelihood is monotone in
#'         \eqn{(\lambda_{\mu}, \lambda_{\log\sigma}, \lambda_{\alpha})}
#'         and would always favour \eqn{\lambda = 0}. Note: this score
#'         is NOT a true frequentist IC -- it is the
#'         penalised-objective AIC, an effective-fit score.
#' }
#'
#' @return A \code{data.frame} with one row per fit and columns:
#' \describe{
#'   \item{\code{model}}{Name (or call index) of the fit.}
#'   \item{\code{likelihood}}{Likelihood family (\code{x$likelihood}).}
#'   \item{\code{loglik}}{Un-penalised log-likelihood at convergence.}
#'   \item{\code{df}}{Number of free parameters.}
#'   \item{\code{nobs}}{Sample size (occurrences).}
#'   \item{\code{AIC}, \code{BIC}}{Information criteria.}
#'   \item{\code{dAIC}, \code{dBIC}}{Deltas vs. the minimum.}
#'   \item{\code{weight_AIC}}{Akaike weights, \eqn{\exp(-\Delta_i / 2) /
#'         \sum_j \exp(-\Delta_j / 2)}.}
#'   \item{\code{convergence}}{ucminf convergence code.}
#' }
#'
#' @section Caveats:
#'
#' AIC and BIC compare *likelihoods* on the SAME observed data. When
#' comparing weighted vs. presence-only families, both families evaluate
#' \eqn{\log L} on \eqn{n_{occ}} occurrences; the weighted family
#' additionally normalises by an \code{env_m} integral, so its
#' likelihood is on a different scale. In this implementation,
#' \code{logLik(weighted)} returns the full weighted-likelihood value
#' (occurrence sum minus log-denominator), which is the right thing
#' for **within-family** comparison (e.g., \code{weighted} vs.
#' \code{skew_normal_weighted}) but should be interpreted with care
#' when contrasting weighted with presence-only.
#'
#' @export
#' @examples
#' \dontrun{
#' fit_w   <- optimize_niche(env_occ, env_m, likelihood = "weighted")
#' fit_skn <- optimize_niche(env_occ, env_m, likelihood = "skew_normal_weighted")
#' fit_skt <- optimize_niche(env_occ, env_m, likelihood = "skew_t_weighted")
#' compare_nicher(weighted = fit_w, skew_normal = fit_skn, skew_t = fit_skt)
#' }
compare_nicher <- function(..., sort_by = c("AIC", "BIC", "loglik"),
                           comparison_basis = c("unpenalised",
                                                "penalised")) {
  sort_by          <- match.arg(sort_by)
  comparison_basis <- match.arg(comparison_basis)
  args <- list(...)
  # Support both compare_nicher(a, b, c) and compare_nicher(list(a=a, b=b))
  if (length(args) == 1L && is.list(args[[1L]]) &&
      !inherits(args[[1L]], "nicher")) {
    args <- args[[1L]]
  }
  if (length(args) < 2L) {
    stop("compare_nicher() requires at least two `nicher` objects.")
  }
  is_nicher <- vapply(args, inherits, logical(1L), what = "nicher")
  if (!all(is_nicher)) {
    stop("All arguments to compare_nicher() must be `nicher` objects.")
  }
  nm <- names(args)
  if (is.null(nm)) nm <- rep("", length(args))
  nm[nm == ""] <- paste0("fit", which(nm == ""))

  # Cross-validate that fits are comparable: same nobs and var_names.
  nobs_vec <- vapply(args, function(x) {
    if (is.null(x$best$nobs)) NA_integer_ else as.integer(x$best$nobs)
  }, integer(1L))
  if (any(is.na(nobs_vec))) {
    stop("One or more fits lack `best$nobs` -- refit with nicher >= 3.2.0.")
  }
  if (length(unique(nobs_vec)) > 1L) {
    stop("Fits have different `nobs` (", paste(nobs_vec, collapse = ", "),
         "); IC values are not comparable.")
  }
  vn0 <- args[[1L]]$var_names
  for (i in seq_along(args)[-1L]) {
    if (!identical(args[[i]]$var_names, vn0)) {
      stop("Fits have different `var_names`; IC values are not comparable.")
    }
  }

  # env_occ fingerprint: must match exactly across fits (same data).
  occ_fps <- lapply(args, function(x) x$best$env_occ_fingerprint)
  if (any(vapply(occ_fps, is.null, logical(1L)))) {
    stop("One or more fits lack `best$env_occ_fingerprint` -- refit with ",
         "nicher >= 3.2.0.")
  }
  for (i in seq_along(args)[-1L]) {
    if (!identical(occ_fps[[i]], occ_fps[[1L]])) {
      stop("Fits used different `env_occ`; IC values are not comparable.")
    }
  }

  # env_m fingerprint: optional (presence_only / skew_normal / skew_t do
  # not consume env_m). Among fits that DO carry an env_m fingerprint,
  # all must match. If some fits used env_m and others did not, warn.
  m_fps   <- lapply(args, function(x) x$best$env_m_fingerprint)
  has_m   <- !vapply(m_fps, is.null, logical(1L))
  if (any(has_m)) {
    nz <- which(has_m)
    ref_m <- m_fps[[nz[1L]]]
    for (j in nz[-1L]) {
      if (!identical(m_fps[[j]], ref_m)) {
        stop("Fits used different `env_m`; IC values are not comparable.")
      }
    }
    if (!all(has_m)) {
      warning("Mixing fits with and without `env_m` (likelihood families ",
              paste(unique(vapply(args, `[[`, character(1L), "likelihood")),
                    collapse = ", "),
              "). AIC/BIC values are still computed but interpret with ",
              "care: weighted families normalise the likelihood by an ",
              "env_m integral, so their loglik is on a different scale ",
              "than presence_only / skew_normal / skew_t.",
              call. = FALSE)
    }
  }

  rows <- lapply(seq_along(args), function(i) {
    x   <- args[[i]]
    ll  <- logLik(x)
    df  <- attr(ll, "df")
    n   <- attr(ll, "nobs")
    ll_basis <- if (comparison_basis == "penalised") {
      as.numeric(x$best$loglik)
    } else {
      as.numeric(ll)
    }
    data.frame(
      model       = nm[i],
      likelihood  = x$likelihood,
      loglik      = ll_basis,
      df          = as.integer(df),
      nobs        = as.integer(n),
      AIC         = -2 * ll_basis + 2 * df,
      BIC         = -2 * ll_basis + log(n) * df,
      convergence = as.integer(x$best$convergence),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "comparison_basis") <- comparison_basis

  out$dAIC <- out$AIC - min(out$AIC, na.rm = TRUE)
  out$dBIC <- out$BIC - min(out$BIC, na.rm = TRUE)

  w_unnorm <- exp(-out$dAIC / 2)
  out$weight_AIC <- w_unnorm / sum(w_unnorm, na.rm = TRUE)

  ord <- switch(sort_by,
    AIC    = order(out$AIC),
    BIC    = order(out$BIC),
    loglik = order(-out$loglik)
  )
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out <- out[, c("model", "likelihood", "loglik", "df", "nobs",
                 "AIC", "dAIC", "BIC", "dBIC", "weight_AIC",
                 "convergence")]
  attr(out, "comparison_basis") <- comparison_basis
  out
}
