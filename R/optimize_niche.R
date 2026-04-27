#' Optimize niche model log-likelihood with multi-start Sobol design
#'
#' Runs multi-start optimization over a Sobol low-discrepancy sequence
#' (via \pkg{pomp}) of starting points covering the parameter space implied
#' by \code{env_occ} and \code{breadth}.
#'
#' Two likelihood models are supported:
#' \itemize{
#'   \item \code{"weighted"}: weighted-normal model with KDE correction
#'         (Jiménez & Soberón 2022).
#'   \item \code{"presence_only"}: model using only presence points,
#'         no background correction.
#' }
#'
#' @section Backends:
#'
#' \describe{
#'   \item{\code{backend = "cpp"} (default)}{Optimizes via
#'     \code{ucminfcpp::ucminf_xptr()} with a compiled C++ objective built by
#'     \code{create_niche_obj_ptr()}. The full \code{(theta -> mu, log sigma, v)}
#'     unpacking, \code{cvine_cholesky}, log-likelihood, and gradient are
#'     evaluated in pure C++ with no R-callback overhead.
#'     For the weighted likelihood, gradient mode \code{"analytic"} uses
#'     closed-form derivatives over \code{mu} and \code{log sigma} and central
#'     finite differences over the C-vine partial-correlation block (hybrid
#'     analytic gradient).}
#'   \item{\code{backend = "r"} (deprecated)}{Optimizes via
#'     \code{ucminf::ucminf()} with the legacy R-level objective functions
#'     (\code{loglik_niche_math_weighted()} or
#'     \code{loglik_niche_math_presence_only()}). Provided ONLY to allow
#'     side-by-side benchmarking with the C++ backend; will be removed in
#'     the next minor release.}
#' }
#'
#' @section KDE sampling (weighted model):
#'
#' The KDE that weights the M-background depends only on the environmental
#' data (not on \code{theta}), so the weights are computed once before
#' optimization. By default, the KDE reference (\code{m_kde_subsample}) and
#' the denominator subset (\code{m_subsample}) are both capped at
#' 10\,000 distinct background combinations -- the maximum the model can
#' usefully exploit even for very large rasters. A floor of
#' \code{max(500, 50 \cdot 2^p)} rows is used as a heuristic minimum
#' representative sample (Silverman 1986; Wand & Jones 1995); below this
#' floor a \code{warning()} is emitted.
#'
#' @param env_occ Data frame of environmental values at presence points.
#' @param env_m Data frame of background environmental values. Required for
#'   \code{likelihood = "weighted"}; ignored for \code{"presence_only"}.
#' @param num_starts Integer. Number of Sobol starting points.
#' @param breadth Numeric in (0, 0.5). Controls the quantile range used to
#'   define starting bounds for \code{mu} parameters. Default \code{0.1}.
#' @param likelihood One of \code{"weighted"} (default) or
#'   \code{"presence_only"}.
#' @param backend One of \code{"cpp"} (default) or \code{"r"} (deprecated;
#'   see Backends section).
#' @param grad Gradient strategy: \code{"auto"} (default) selects
#'   \code{"analytic"} for the weighted model with \code{backend = "cpp"} and
#'   \code{"central"} otherwise. Force one of
#'   \code{c("analytic", "central", "forward")} to override.
#' @param m_subsample,m_kde_subsample Optional integer or fraction in (0, 1].
#'   Resolved to \code{min(nrow(env_m), 10000)} when \code{NULL} (default).
#' @param seed Optional integer to make subsampling deterministic.
#' @param control Named list of control parameters for
#'   \code{ucminfcpp::ucminf_xptr()} (cpp backend) or \code{ucminf::ucminf()}
#'   (r backend). Recognized entries:
#'   \describe{
#'     \item{\code{grad}}{"central" (default)}
#'     \item{\code{gradstep}}{c(1e-6, 1e-8)}
#'     \item{\code{grtol}}{1e-4}
#'     \item{\code{xtol}}{1e-8}
#'     \item{\code{stepmax}}{5}
#'     \item{\code{maxeval}}{2000}
#'   }
#' @param verbose Logical. If \code{TRUE}, print per-start progress.
#' @param ... Additional arguments forwarded to the objective function
#'   (e.g. \code{eta}).
#'
#' @return An object of class \code{"nicher"} (see \code{\link{new_nicher}}).
#'
#' @export
#' @examples
#' \dontrun{
#' result <- optimize_niche(
#'   env_occ    = example_env_occ_2d,
#'   env_m      = example_env_m_2d,
#'   num_starts = 20L,
#'   breadth    = 0.1,
#'   likelihood = "weighted"
#' )
#' print(result)
#' assess(result)
#' }
optimize_niche <- function(env_occ,
                           env_m,
                           num_starts = 100L,
                           breadth    = 0.1,
                           likelihood = c("weighted", "presence_only"),
                           backend    = c("cpp", "r"),
                           grad       = c("auto", "analytic",
                                          "central", "forward"),
                           m_subsample     = NULL,
                           m_kde_subsample = NULL,
                           seed            = NULL,
                           control = list(),
                           verbose = FALSE,
                           ...) {
  likelihood <- match.arg(likelihood)
  backend    <- match.arg(backend)
  grad       <- match.arg(grad)

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  if (likelihood != "presence_only") {
    if (missing(env_m) || is.null(env_m)) {
      stop("env_m must be provided for likelihood '", likelihood, "'")
    }
    if (!identical(sort(colnames(env_occ)), sort(colnames(env_m)))) {
      stop("env_occ and env_m must have the same",
           " variables (column names)")
    }
  }
  if (!is.numeric(breadth) || length(breadth) != 1L ||
      breadth <= 0 || breadth >= 0.5) {
    stop("breadth must be a single number in (0, 0.5)")
  }

  if (backend == "r") {
    lifecycle::deprecate_soft(
      when = "2.1.0",
      what = I('optimize_niche(backend = "r")'),
      details = paste0(
        'The legacy R-objective backend will be removed in nicher 2.2.0. ',
        'Use backend = "cpp" (default) for production runs; the "r" backend ',
        'is retained only for side-by-side benchmarking via ',
        'benchmark_optimize_niche().'
      )
    )
  }

  # Resolve grad="auto"
  resolved_grad <- if (grad == "auto") {
    if (likelihood == "weighted" && backend == "cpp") "analytic" else "central"
  } else grad

  # ------------------------------------------------------------------
  # Sobol starts
  # ------------------------------------------------------------------
  quant_vec <- c(breadth, 0.5, 1.0 - breadth)
  starts_df <- start_theta_multiple(
    env_data   = env_occ,
    num_starts = num_starts,
    quant_vec  = quant_vec,
    method     = "sobol"
  )
  starts_list <- split(starts_df, seq_len(nrow(starts_df)))
  starts_list <- lapply(starts_list, function(x) {
    v <- as.numeric(x)
    names(v) <- colnames(starts_df)
    v
  })

  # ------------------------------------------------------------------
  # Default control list
  # ------------------------------------------------------------------
  default_ctrl <- list(
    grad     = "central",
    gradstep = c(1e-6, 1e-8),
    grtol    = 1e-4,
    xtol     = 1e-8,
    stepmax  = 5,
    maxeval  = 2000L
  )
  ctrl <- utils::modifyList(default_ctrl, control)

  # ------------------------------------------------------------------
  # Resolve KDE subsamples / weights ONCE per fit (weighted only)
  # ------------------------------------------------------------------
  weighted_inputs <- NULL
  if (likelihood == "weighted") {
    weighted_inputs <- .resolve_weighted_inputs(
      env_occ         = env_occ,
      env_m           = env_m,
      m_subsample     = m_subsample,
      m_kde_subsample = m_kde_subsample,
      seed            = seed
    )
  }

  # ------------------------------------------------------------------
  # Run all starts
  # ------------------------------------------------------------------
  helper <- if (backend == "cpp") .optimize_niche_helper_cpp
            else                  .optimize_niche_helper_r

  results <- vector("list", length(starts_list))
  for (i in seq_along(starts_list)) {
    if (verbose) message(sprintf("Start %d / %d", i, length(starts_list)))
    results[[i]] <- helper(
      param           = starts_list[[i]],
      env_occ         = env_occ,
      env_m           = env_m,
      control         = ctrl,
      likelihood      = likelihood,
      grad            = resolved_grad,
      weighted_inputs = weighted_inputs,
      ...
    )
  }

  # ------------------------------------------------------------------
  # Compile solutions
  # ------------------------------------------------------------------
  solutions <- data.frame(
    start_id    = seq_along(results),
    loglik      = vapply(results,
                         function(x) as.numeric(x$loglik), numeric(1L)),
    convergence = vapply(results,
                         function(x) as.integer(x$convergence), integer(1L)),
    stringsAsFactors = FALSE
  )
  solutions$full_par <- lapply(results, function(x) x$theta)

  ord       <- order(solutions$loglik, decreasing = TRUE)
  solutions <- solutions[ord, ]
  rownames(solutions) <- NULL

  # Prefer the best CONVERGED start (ucminf codes 1 / 2); fall back to the
  # overall best if nothing converged. Mirrors the safeguard previously
  # introduced in PR #31.
  conv_mask <- solutions$convergence %in% c(1L, 2L)
  best_idx  <- if (any(conv_mask)) which(conv_mask)[1L] else 1L

  best <- list(
    theta       = solutions$full_par[[best_idx]],
    loglik      = solutions$loglik[best_idx],
    convergence = solutions$convergence[best_idx]
  )

  # ------------------------------------------------------------------
  # Internal safeguard: validate cpp-backend xptr result against ucminf
  # ------------------------------------------------------------------
  if (backend == "cpp" && best$convergence %in% c(1L, 2L)) {
    .validate_xptr_result(
      best            = best,
      env_occ         = env_occ,
      env_m           = env_m,
      likelihood      = likelihood,
      weighted_inputs = weighted_inputs,
      ctrl            = ctrl,
      ...
    )
  }

  if (verbose) {
    message(sprintf(
      "Best log-likelihood: %.6f (convergence = %d)",
      best$loglik, best$convergence
    ))
  }

  new_nicher(
    solutions  = solutions,
    best       = best,
    likelihood = likelihood,
    n_starts   = num_starts,
    var_names  = colnames(env_occ)
  )
}


# ===========================================================================
# Internal: KDE sample-size policy and weight precomputation
# ===========================================================================

#' Resolve KDE subsampling indices and precompute KDE weights ONCE per fit.
#'
#' Implements the package's KDE sampling policy:
#' \itemize{
#'   \item Hard cap of 10000 distinct background combinations.
#'   \item Heuristic minimum representative sample of
#'     \code{max(500, 50 * 2^p)} rows; below this a \code{warning()} is fired.
#'   \item Optional fraction (\code{x < 1}) or integer count for
#'     \code{m_subsample} / \code{m_kde_subsample}.
#' }
#'
#' Returns a list with: \code{den_idx}, \code{kde_idx}, \code{w_occ},
#' \code{w_den}, \code{n_m}.
#'
#' @keywords internal
.resolve_weighted_inputs <- function(env_occ, env_m,
                                     m_subsample     = NULL,
                                     m_kde_subsample = NULL,
                                     seed            = NULL) {
  occ_mat <- as.matrix(env_occ)
  m_mat   <- as.matrix(env_m)
  p       <- ncol(occ_mat)
  n_m     <- nrow(m_mat)
  cap     <- 10000L
  floor_n <- max(500L, 50L * (2L ^ p))

  pick_size <- function(x, nmax) {
    if (is.null(x)) return(min(nmax, cap))
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
      stop("m_subsample/m_kde_subsample must be a single positive number.")
    }
    n <- if (x < 1) max(1L, floor(x * nmax))
         else       min(nmax, as.integer(round(x)))
    min(n, cap)
  }

  if (!is.null(seed)) set.seed(seed)

  n_den <- pick_size(m_subsample,     n_m)
  n_kde <- pick_size(m_kde_subsample, n_m)

  if (n_kde < floor_n) {
    warning(sprintf(
      paste0("KDE reference sample (n_kde = %d) is below the recommended ",
             "minimum %d for p = %d; KDE weights may be unstable."),
      n_kde, floor_n, p
    ))
  }

  den_idx <- if (n_den < n_m) sample.int(n_m, n_den) else seq_len(n_m)
  kde_idx <- if (n_kde < n_m) sample.int(n_m, n_kde) else seq_len(n_m)

  M_kde <- m_mat[kde_idx, , drop = FALSE]
  w_occ <- as.numeric(kde_gaussian(occ_mat, M_kde))
  w_den <- as.numeric(kde_gaussian(m_mat[den_idx, , drop = FALSE], M_kde))

  list(
    den_idx = as.integer(den_idx),
    kde_idx = as.integer(kde_idx),
    w_occ   = w_occ,
    w_den   = w_den,
    n_m     = n_m
  )
}


# ===========================================================================
# Internal safeguard: cross-check xptr result with ucminf::ucminf
# ===========================================================================

#' Validate ucminfcpp::ucminf_xptr result against ucminf::ucminf.
#'
#' Runs a short \code{ucminf::ucminf} optimization from the best theta
#' and warns if the two log-likelihoods differ by more than \code{1e-3}.
#' This guards against pointer-safety issues in the C++ backend.
#'
#' @keywords internal
.validate_xptr_result <- function(best, env_occ, env_m,
                                  likelihood, ctrl,
                                  weighted_inputs = NULL, ...) {
  # Forward the SAME subsampled KDE inputs used by the optimizer so that the
  # reference ucminf::ucminf evaluates the identical objective function.
  # Otherwise (default cap = 10 000) any env_m larger than the cap would
  # cause optimizer and validator to optimize different functions and fire
  # spurious "pointer-safety" warnings.
  fn <- switch(likelihood,
    presence_only = function(theta) {
      loglik_niche_math_presence_only(theta, env_occ = env_occ,
                                      neg = TRUE, ...)
    },
    weighted = function(theta) {
      loglik_niche_math_weighted(
        theta, env_occ = env_occ, env_m = env_m, neg = TRUE,
        den_idx       = weighted_inputs$den_idx,
        kde_idx       = weighted_inputs$kde_idx,
        precomp_w_den = weighted_inputs$w_den,
        ...
      )
    }
  )
  ref <- tryCatch({
    res <- ucminf::ucminf(par = best$theta, fn = fn,
                          control = list(maxeval = 500L), hessian = FALSE)
    -res$value
  }, error = function(e) NA_real_)

  if (is.finite(ref) && abs(best$loglik - ref) > 1e-3) {
    warning(sprintf(
      paste0("ucminfcpp::ucminf_xptr and ucminf::ucminf disagree: ",
             "xptr loglik = %.6f, ucminf loglik = %.6f. ",
             "Possible pointer-safety issue."),
      best$loglik, ref
    ))
  }
  invisible(NULL)
}


# ===========================================================================
# Internal: cpp-backend per-start helper (ucminfcpp::ucminf_xptr)
# ===========================================================================

#' Run ucminfcpp::ucminf_xptr() for a single starting vector.
#' @keywords internal
.optimize_niche_helper_cpp <- function(param, env_occ, env_m, control,
                                       likelihood, grad,
                                       weighted_inputs = NULL,
                                       ...) {
  param_names <- names(param)
  param <- as.numeric(param)
  if (!is.null(param_names)) names(param) <- param_names
  if (any(!is.finite(param))) {
    stop("All starting parameters must be finite")
  }

  dots <- list(...)
  eta <- if (!is.null(dots$eta)) dots$eta else 1.0

  env_occ_mat <- as.matrix(env_occ)
  env_m_mat   <- if (!is.null(env_m)) as.matrix(env_m) else NULL

  den_idx <- kde_idx <- precomp_w_occ <- precomp_w_den <- NULL
  if (likelihood == "weighted" && !is.null(weighted_inputs)) {
    den_idx       <- weighted_inputs$den_idx
    kde_idx       <- weighted_inputs$kde_idx
    precomp_w_occ <- weighted_inputs$w_occ
    precomp_w_den <- weighted_inputs$w_den
  }

  gs <- if (!is.null(control$gradstep)) control$gradstep else c(1e-6, 1e-8)
  xptr <- create_niche_obj_ptr(
    env_occ       = env_occ_mat,
    env_m         = env_m_mat,
    eta           = eta,
    likelihood    = likelihood,
    den_idx       = den_idx,
    kde_idx       = kde_idx,
    precomp_w_occ = precomp_w_occ,
    precomp_w_den = precomp_w_den,
    grad          = grad,
    gradstep      = gs
  )

  control_args <- control
  if (is.null(control_args$grad))     control_args$grad <- "central"
  if (is.null(control_args$gradstep)) control_args$gradstep <- gs
  if (!is.null(control_args$maxeval)) control_args$maxeval <- as.integer(control_args$maxeval)
  con <- do.call(ucminfcpp::ucminf_control, control_args)

  out <- tryCatch({
    res <- ucminfcpp::ucminf_xptr(par = param, xptr = xptr, control = con)
    if (is.null(names(res$par)) && !is.null(names(param))) {
      names(res$par) <- names(param)
    }
    if (is.null(res$convergence)) res$convergence <- NA_integer_
    res
  }, error = function(e) {
    list(par = param, value = Inf, convergence = NA_integer_,
         error = conditionMessage(e))
  })

  list(theta = out$par, loglik = -out$value,
       convergence = out$convergence)
}


# ===========================================================================
# Internal: r-backend per-start helper (legacy ucminf::ucminf)
# ===========================================================================

#' Run ucminf::ucminf() with the legacy R objective function for a single
#' starting vector.
#'
#' Used only when \code{backend = "r"}; will be removed alongside that flag
#' in nicher 2.2.0.
#' @keywords internal
.optimize_niche_helper_r <- function(param, env_occ, env_m, control,
                                     likelihood, grad,
                                     weighted_inputs = NULL,
                                     ...) {
  param_names <- names(param)
  param <- as.numeric(param)
  if (!is.null(param_names)) names(param) <- param_names
  if (any(!is.finite(param))) stop("All starting parameters must be finite")

  dots <- list(...)
  eta <- if (!is.null(dots$eta)) dots$eta else 1.0

  fn <- if (likelihood == "weighted") {
    den_idx       <- if (!is.null(weighted_inputs)) weighted_inputs$den_idx else NULL
    kde_idx       <- if (!is.null(weighted_inputs)) weighted_inputs$kde_idx else NULL
    precomp_w_den <- if (!is.null(weighted_inputs)) weighted_inputs$w_den   else NULL
    function(theta) {
      loglik_niche_math_weighted(
        theta, env_occ = env_occ, env_m = env_m, eta = eta, neg = TRUE,
        den_idx = den_idx, kde_idx = kde_idx,
        precomp_w_den = precomp_w_den
      )
    }
  } else {
    function(theta) {
      loglik_niche_math_presence_only(
        theta, env_occ = env_occ, eta = eta, neg = TRUE
      )
    }
  }

  ucminf_ctrl <- list(
    grtol   = if (!is.null(control$grtol))   control$grtol   else 1e-4,
    xtol    = if (!is.null(control$xtol))    control$xtol    else 1e-8,
    stepmax = if (!is.null(control$stepmax)) control$stepmax else 5,
    maxeval = if (!is.null(control$maxeval)) as.integer(control$maxeval) else 2000L
  )

  out <- tryCatch({
    res <- ucminf::ucminf(par = param, fn = fn, hessian = FALSE,
                          control = ucminf_ctrl)
    if (is.null(res$convergence)) res$convergence <- NA_integer_
    res
  }, error = function(e) {
    list(par = param, value = Inf, convergence = NA_integer_,
         error = conditionMessage(e))
  })

  par_out <- out$par
  if (is.null(names(par_out)) && !is.null(param_names)) names(par_out) <- param_names
  list(theta = par_out, loglik = -out$value,
       convergence = out$convergence)
}
