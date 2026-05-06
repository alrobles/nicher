# R/cv.R
# Cross-validation for nicher fits.
#
# Exports: cv_nicher.

# -----------------------------------------------------------------------
# Build a negative-log-likelihood closure (penalised, mirroring the
# active fit's penalty knobs) that consumes a chosen `env_occ_sub`.
# Used by `cv_nicher()` to (a) refit on each train fold using ucminf
# warm-started from the full-data theta and (b) score the held-out
# fold using `.compute_unpenalised_loglik()` (already in R/ic.R).
#
# Returns a function `fn(theta) -> neg_loglik` that closes over all
# the penalty configuration.  Keeps `weighted_inputs` simple: for the
# weighted families we recompute KDE weights for the train fold from
# scratch via `.resolve_weighted_inputs()`; this matches what
# `optimize_niche()` itself does and is cheap (one KDE per fold).
# -----------------------------------------------------------------------
.cv_make_neg_loglik <- function(likelihood, env_occ_sub, env_m, args,
                                weighted_inputs_sub) {
  pmc  <- if (!is.null(args$prior_mu_center))        as.numeric(args$prior_mu_center)        else NULL
  plsc <- if (!is.null(args$prior_log_sigma_center)) as.numeric(args$prior_log_sigma_center) else NULL

  M_den <- if (!is.null(env_m) && !is.null(weighted_inputs_sub)) {
    as.matrix(env_m)[weighted_inputs_sub$den_idx, , drop = FALSE]
  } else {
    NULL
  }
  occ_mat <- as.matrix(env_occ_sub)

  switch(
    likelihood,
    presence_only = function(theta) {
      loglik_niche_math_presence_only(
        theta, env_occ = occ_mat, neg = TRUE, eta = args$eta
      )
    },
    kde_bias_corrected = function(theta) {
      loglik_niche_math_kde_bias_corrected(
        theta, env_occ = occ_mat, env_m = as.matrix(env_m), neg = TRUE,
        den_idx       = weighted_inputs_sub$den_idx,
        kde_idx       = weighted_inputs_sub$kde_idx,
        precomp_w_den = weighted_inputs_sub$w_den,
        eta           = args$eta
      )
    },
    weighted = function(theta) {
      loglik_niche_math_weighted_cpp(
        theta = theta, env_occ = occ_mat, M_den = M_den,
        w_occ = weighted_inputs_sub$w_occ,
        w_den = weighted_inputs_sub$w_den,
        prior_log_sigma_center = if (is.null(plsc)) numeric(ncol(occ_mat))
                                 else plsc,
        prior_log_sigma_lambda = args$prior_log_sigma_lambda,
        eta = args$eta,
        prior_mu_center    = pmc,
        prior_mu_lambda    = args$prior_mu_lambda,
        prior_alpha_lambda = args$prior_alpha_lambda
      )
    },
    skew_normal = function(theta) {
      loglik_niche_math_skew_normal_cpp(
        theta = theta, env_occ = occ_mat, eta = args$eta,
        prior_mu_center        = pmc,
        prior_mu_lambda        = args$prior_mu_lambda,
        prior_log_sigma_center = plsc,
        prior_log_sigma_lambda = args$prior_log_sigma_lambda,
        prior_alpha_lambda     = args$prior_alpha_lambda
      )
    },
    skew_normal_weighted = function(theta) {
      loglik_niche_math_skew_normal_weighted_cpp(
        theta = theta, env_occ = occ_mat, M_den = M_den,
        w_occ = weighted_inputs_sub$w_occ,
        w_den = weighted_inputs_sub$w_den,
        prior_log_sigma_center = if (is.null(plsc)) numeric(ncol(occ_mat))
                                 else plsc,
        prior_log_sigma_lambda = args$prior_log_sigma_lambda,
        eta = args$eta,
        prior_mu_center    = pmc,
        prior_mu_lambda    = args$prior_mu_lambda,
        prior_alpha_lambda = args$prior_alpha_lambda
      )
    },
    skew_t = function(theta) {
      loglik_niche_math_skew_t_cpp(
        theta = theta, env_occ = occ_mat, eta = args$eta,
        prior_mu_center        = pmc,
        prior_mu_lambda        = args$prior_mu_lambda,
        prior_log_sigma_center = plsc,
        prior_log_sigma_lambda = args$prior_log_sigma_lambda,
        prior_alpha_lambda     = args$prior_alpha_lambda
      )
    },
    skew_t_weighted = function(theta) {
      loglik_niche_math_skew_t_weighted_cpp(
        theta = theta, env_occ = occ_mat, M_den = M_den,
        w_occ = weighted_inputs_sub$w_occ,
        w_den = weighted_inputs_sub$w_den,
        prior_log_sigma_center = if (is.null(plsc)) numeric(ncol(occ_mat))
                                 else plsc,
        prior_log_sigma_lambda = args$prior_log_sigma_lambda,
        eta = args$eta,
        prior_mu_center    = pmc,
        prior_mu_lambda    = args$prior_mu_lambda,
        prior_alpha_lambda = args$prior_alpha_lambda
      )
    }
  )
}

# -----------------------------------------------------------------------
# Refit on a single train fold, warm-started from `theta_init`.
# Uses ucminf::ucminf with FD gradients; for CV, the per-fit cost
# dominates and the gradient choice is not the bottleneck.
# -----------------------------------------------------------------------
.cv_refit_fold <- function(theta_init, train_occ, env_m, likelihood, args,
                           ucminf_control = list(maxeval = 500L)) {
  weighted_fams <- c("kde_bias_corrected", "weighted",
                     "skew_normal_weighted", "skew_t_weighted")
  weighted_inputs_train <- if (likelihood %in% weighted_fams) {
    .resolve_weighted_inputs(
      env_occ         = train_occ,
      env_m           = env_m,
      m_subsample     = args$m_subsample,
      m_kde_subsample = args$m_kde_subsample,
      seed            = args$seed
    )
  } else NULL

  fn <- .cv_make_neg_loglik(
    likelihood          = likelihood,
    env_occ_sub         = train_occ,
    env_m               = env_m,
    args                = args,
    weighted_inputs_sub = weighted_inputs_train
  )
  res <- ucminf::ucminf(par = theta_init, fn = fn,
                        control = ucminf_control, hessian = FALSE)
  list(theta = as.numeric(res$par), convergence = as.integer(res$convergence))
}

# -----------------------------------------------------------------------
# cv_nicher
# -----------------------------------------------------------------------

#' k-fold (or leave-one-out) cross-validation for a nicher fit
#'
#' Refits the same likelihood family + same penalty knobs as the
#' supplied \code{nicher} object on each train fold of \code{env_occ},
#' warm-started from the full-data \eqn{\hat\theta}, then scores the
#' held-out fold by evaluating the un-penalised log-likelihood at the
#' fold's converged \eqn{\hat\theta_{train}}. Returns the **summed**
#' held-out log-likelihood across folds, which is the standard
#' out-of-sample log-likelihood used to tune ridge regularisation.
#'
#' For penalised fits the un-penalised log-likelihood at the converged
#' \eqn{\hat\theta} (what \code{\link{logLik.nicher}} returns) is monotone
#' in the penalty strengths -- the in-sample loglik is always best at
#' \eqn{\lambda = 0}. Cross-validation breaks that monotonicity by scoring
#' on data the optimiser did not see, so a properly tuned penalty
#' achieves higher out-of-sample loglik than the unpenalised MLE.
#'
#' @param fit A \code{"nicher"} object returned by
#'   \code{\link{optimize_niche}}.
#' @param env_occ A \code{matrix} or \code{data.frame} of occurrences,
#'   identical (up to row order) to the one passed to the original
#'   \code{optimize_niche()} call. Validated against
#'   \code{fit$best$env_occ_fingerprint}.
#' @param env_m A \code{matrix} / \code{data.frame} of background points,
#'   or \code{NULL} for presence-only families. Validated against
#'   \code{fit$best$env_m_fingerprint}.
#' @param type Character. \code{"kfold"} (default) or \code{"loo"}.
#' @param k Integer. Number of folds when \code{type = "kfold"}; ignored
#'   for \code{"loo"}. Default \code{5L}.
#' @param seed Integer or \code{NULL}. Random seed for fold assignment
#'   (\code{type = "kfold"} only). Default \code{NULL} (no seed).
#' @param num_starts_cv Integer. Number of starts per fold. Default
#'   \code{1L} (warm-start from \code{fit$best$theta}). Set higher for
#'   robustness in \code{"kfold"} mode; for \code{"loo"} this should
#'   essentially always be \code{1L} since you are running n_occ
#'   refits.
#' @param verbose Logical. Print per-fold progress.
#' @param ucminf_control Optional list passed to \code{ucminf::ucminf}
#'   for each fold's optimiser; defaults to \code{list(maxeval = 500L)}.
#'
#' @return A list with class \code{"nicher_cv"}:
#' \describe{
#'   \item{\code{cv_loglik}}{Summed held-out log-likelihood across folds.}
#'   \item{\code{cv_loglik_mean}}{Per-occurrence average:
#'     \code{cv_loglik / nobs}.}
#'   \item{\code{per_fold}}{\code{data.frame(fold, n_test, loglik_test,
#'     convergence)}.}
#'   \item{\code{type}}{\code{"kfold"} or \code{"loo"}.}
#'   \item{\code{k}}{Number of folds (\code{n_occ} for LOO).}
#'   \item{\code{seed}}{Seed used (or \code{NA_integer_}).}
#' }
#'
#' @section Caveats:
#'
#' * Random k-fold CV assumes the occurrences are exchangeable. For
#'   spatially autocorrelated data this can be optimistic; spatial-block
#'   CV is left to a future release.
#' * LOO with \code{n_occ} large (e.g. Vicugna's 545) is feasible but
#'   non-trivial: 545 ucminf refits at ~0.1-1 s each = ~1-10 minutes.
#'   Consider \code{type = "kfold", k = 10L} as a faster proxy.
#' * Each fold's refit uses \code{ucminf::ucminf} starting from
#'   \code{fit$best$theta}. If the per-fold likelihood is not unimodal
#'   (rare on dense data, common with skew families on sparse data),
#'   bump \code{num_starts_cv} or use the original
#'   \code{optimize_niche()} multistart machinery.
#'
#' @seealso \code{\link{compare_nicher}} for AIC/BIC-based comparisons.
#' @export
#' @examples
#' \dontrun{
#' fit_w <- optimize_niche(
#'   env_occ = example_env_occ_2d, env_m = example_env_m_2d,
#'   num_starts = 20L, breadth = 0.45, likelihood = "weighted",
#'   prior_mu_lambda = 1.0
#' )
#' cv_nicher(fit_w, example_env_occ_2d, example_env_m_2d,
#'           type = "kfold", k = 5L, seed = 42L)
#' }
cv_nicher <- function(fit, env_occ, env_m = NULL,
                      type = c("kfold", "loo"),
                      k = 5L, seed = NULL,
                      num_starts_cv = 1L,
                      verbose = FALSE,
                      ucminf_control = list(maxeval = 500L)) {
  if (!inherits(fit, "nicher")) {
    stop("`fit` must be a `nicher` object returned by optimize_niche().")
  }
  if (is.null(fit$fit_args)) {
    stop("This nicher object lacks `fit_args`. Refit with the current ",
         "version to use cv_nicher().")
  }
  type <- match.arg(type)

  args        <- fit$fit_args
  likelihood  <- args$likelihood
  weighted_fams <- c("kde_bias_corrected", "weighted",
                     "skew_normal_weighted", "skew_t_weighted")
  uses_env_m  <- likelihood %in% weighted_fams

  # ---- validate env_occ / env_m against stored fingerprints ----------
  if (is.null(fit$best$env_occ_fingerprint)) {
    stop("Fit lacks `env_occ_fingerprint`; refit with the current version.")
  }
  if (!identical(.env_fingerprint(env_occ), fit$best$env_occ_fingerprint)) {
    stop("`env_occ` does not match the data the fit was trained on ",
         "(fingerprint mismatch).")
  }
  if (uses_env_m) {
    if (is.null(env_m)) {
      stop("`env_m` is required for likelihood '", likelihood, "'.")
    }
    if (!identical(.env_fingerprint(env_m), fit$best$env_m_fingerprint)) {
      stop("`env_m` does not match the data the fit was trained on ",
           "(fingerprint mismatch).")
    }
  }

  if (!is.numeric(num_starts_cv) || length(num_starts_cv) != 1L ||
      num_starts_cv < 1L || !is.finite(num_starts_cv)) {
    stop("`num_starts_cv` must be a positive integer.")
  }
  num_starts_cv <- as.integer(num_starts_cv)

  # ---- fold assignment ----------------------------------------------
  occ_mat <- as.matrix(env_occ)
  n_occ   <- nrow(occ_mat)
  if (type == "loo") {
    folds <- as.list(seq_len(n_occ))
    k_eff <- n_occ
    seed_used <- NA_integer_
  } else {
    if (!is.numeric(k) || length(k) != 1L || k < 2L || k > n_occ) {
      stop(sprintf("`k` must be an integer in [2, n_occ=%d].", n_occ))
    }
    k <- as.integer(k)
    seed_used <- if (is.null(seed)) NA_integer_ else as.integer(seed)
    if (!is.null(seed)) set.seed(seed)
    perm <- sample.int(n_occ)
    folds <- split(perm, cut(seq_along(perm), breaks = k, labels = FALSE))
    folds <- lapply(folds, sort)
    k_eff <- k
  }

  theta_init <- fit$best$theta

  # ---- per-fold refit + score ----------------------------------------
  per_fold <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(n_occ), test_idx)
    train_occ <- occ_mat[train_idx, , drop = FALSE]
    test_occ  <- occ_mat[test_idx,  , drop = FALSE]
    if (!is.null(colnames(occ_mat))) {
      colnames(train_occ) <- colnames(occ_mat)
      colnames(test_occ)  <- colnames(occ_mat)
    }

    if (verbose) {
      message(sprintf("[cv_nicher] fold %d/%d (n_test = %d) ...",
                      i, k_eff, length(test_idx)))
    }

    refit <- if (num_starts_cv == 1L) {
      .cv_refit_fold(theta_init, train_occ, env_m, likelihood, args,
                     ucminf_control = ucminf_control)
    } else {
      # Multi-start refit: use Sobol design via optimize_niche() recursively.
      # We restrict noise on the optimizer by passing the full set of
      # original args except `seed` (which we vary by fold for diversity).
      sub <- do.call(optimize_niche, c(
        list(env_occ = train_occ, env_m = env_m,
             num_starts = num_starts_cv,
             breadth    = args$breadth,
             likelihood = likelihood,
             grad       = args$grad,
             m_subsample            = args$m_subsample,
             m_kde_subsample        = args$m_kde_subsample,
             seed                   = if (is.null(seed)) NULL else seed + i,
             warm_start             = args$warm_start,
             prior_log_sigma_lambda = args$prior_log_sigma_lambda,
             prior_log_sigma_center = args$prior_log_sigma_center,
             prior_mu_lambda        = args$prior_mu_lambda,
             prior_mu_center        = args$prior_mu_center,
             prior_alpha_lambda     = args$prior_alpha_lambda,
             control                = args$control,
             verbose                = FALSE,
             eta                    = args$eta)
      ))
      list(theta = sub$best$theta, convergence = sub$best$convergence)
    }

    # Score test fold: un-penalised loglik at theta_train.
    weighted_inputs_score <- if (uses_env_m) {
      # KDE bandwidth derived from env_m as in the full fit; we reuse
      # `args$m_subsample` / `m_kde_subsample` but compute on test occ.
      .resolve_weighted_inputs(
        env_occ         = test_occ,
        env_m           = env_m,
        m_subsample     = args$m_subsample,
        m_kde_subsample = args$m_kde_subsample,
        seed            = args$seed
      )
    } else NULL

    score <- tryCatch(
      .compute_unpenalised_loglik(
        theta           = refit$theta,
        likelihood      = likelihood,
        env_occ         = test_occ,
        env_m           = env_m,
        weighted_inputs = weighted_inputs_score,
        eta             = args$eta
      ),
      error = function(e) NA_real_
    )

    per_fold[[i]] <- data.frame(
      fold        = i,
      n_test      = length(test_idx),
      loglik_test = score,
      convergence = refit$convergence,
      stringsAsFactors = FALSE
    )
  }
  per_fold <- do.call(rbind, per_fold)

  cv_total <- sum(per_fold$loglik_test, na.rm = FALSE)
  out <- list(
    cv_loglik       = cv_total,
    cv_loglik_mean  = cv_total / n_occ,
    per_fold        = per_fold,
    type            = type,
    k               = k_eff,
    seed            = seed_used,
    likelihood      = likelihood
  )
  class(out) <- "nicher_cv"
  out
}

#' @export
print.nicher_cv <- function(x, ...) {
  cat("-- nicher_cv result --\n")
  cat("  Type     :", x$type, ifelse(x$type == "kfold",
                                     paste0("(k = ", x$k, ")"), ""), "\n")
  cat("  Family   :", x$likelihood, "\n")
  cat("  Folds    :", nrow(x$per_fold), "\n")
  cat("  CV loglik:", round(x$cv_loglik, 4L),
      "(per-occ:", round(x$cv_loglik_mean, 4L), ")\n")
  if (any(is.na(x$per_fold$loglik_test))) {
    cat("  Warning  :", sum(is.na(x$per_fold$loglik_test)),
        "fold(s) returned NA loglik.\n")
  }
  invisible(x)
}
