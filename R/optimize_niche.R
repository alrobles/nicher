#' Optimize niche model log-likelihood with multi-start Sobol design
#'
#' Runs multi-start optimization using \code{ucminfcpp::ucminf_xptr()} with
#' a compiled C++ objective function. Starting points are drawn from a
#' Sobol low-discrepancy sequence (via \pkg{pomp}) over the parameter space
#' implied by \code{env_occ} and \code{breadth}.
#'
#' Three likelihood models are supported:
#' \itemize{
#'   \item \code{"unweighted"}: ellipsoid model with background correction
#'         (Jiménez & Soberón 2019).
#'   \item \code{"weighted"}: weighted-normal model with KDE correction
#'         (Jiménez & Soberón 2022).
#'   \item \code{"presence_only"}: model using only presence points.
#' }
#'
#' An internal safeguard re-optimizes the best result with
#' \code{ucminf::ucminf} and warns if the two log-likelihoods diverge
#' by more than \code{1e-3}, indicating a potential pointer-safety issue.
#'
#' @param env_occ Data frame of environmental values at presence points.
#' @param env_m Data frame of background environmental values. Ignored
#'   when \code{likelihood = "presence_only"}.
#' @param num_starts Integer. Number of Sobol starting points.
#' @param breadth Numeric in (0, 0.5). Controls the quantile range used
#'   to define starting bounds for \code{mu} parameters.
#'   Quantiles are set to \code{c(breadth, 0.5, 1 - breadth)}.
#'   Default is \code{0.1}.
#' @param likelihood Character. One of \code{"unweighted"} (default),
#'   \code{"weighted"}, or \code{"presence_only"}.
#' @param control Named list of control parameters for
#'   \code{ucminfcpp::ucminf_xptr()}. Recognized entries:
#'   \describe{
#'     \item{\code{grad}}{"central" (default)}
#'     \item{\code{gradstep}}{c(1e-6, 1e-8)}
#'     \item{\code{grtol}}{1e-4}
#'     \item{\code{xtol}}{1e-8}
#'     \item{\code{stepmax}}{5}
#'     \item{\code{maxeval}}{2000}
#'   }
#' @param verbose Logical. If \code{TRUE}, print per-start progress.
#' @param ... Additional arguments forwarded to the C++ helper, e.g.
#'   \code{eta} for the LKJ prior, or \code{m_subsample},
#'   \code{m_kde_subsample}, \code{seed} for the weighted model.
#'
#' @return An object of class \code{"nicher"} (see \code{\link{new_nicher}})
#'   with components:
#'   \describe{
#'     \item{\code{solutions}}{Data frame with columns \code{start_id},
#'       \code{loglik}, \code{convergence}, \code{full_par}.}
#'     \item{\code{best}}{List: \code{theta}, \code{loglik},
#'       \code{convergence}.}
#'     \item{\code{likelihood}}{The likelihood model used.}
#'     \item{\code{n_starts}}{Number of starting points.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' result <- optimize_niche(
#'   env_occ    = example_env_occ_2d,
#'   env_m      = example_env_m_2d,
#'   num_starts = 20L,
#'   breadth    = 0.1,
#'   likelihood = "unweighted"
#' )
#' print(result)
#' assess(result)
#' }
optimize_niche <- function(env_occ,
                           env_m,
                           num_starts = 100L,
                           breadth    = 0.1,
                           likelihood = c(
                             "unweighted",
                             "weighted",
                             "presence_only"
                           ),
                           control  = list(),
                           verbose  = FALSE,
                           ...) {
  likelihood <- match.arg(likelihood)

  # ----------------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------------
  if (likelihood != "presence_only") {
    if (missing(env_m) || is.null(env_m)) {
      stop(
        "env_m must be provided for likelihood '",
        likelihood, "'"
      )
    }
    if (!identical(
      sort(colnames(env_occ)),
      sort(colnames(env_m))
    )) {
      stop(
        "env_occ and env_m must have the same",
        " variables (column names)"
      )
    }
  }

  if (!is.numeric(breadth) ||
      length(breadth) != 1L ||
      breadth <= 0 || breadth >= 0.5) {
    stop("breadth must be a single number in (0, 0.5)")
  }

  # ----------------------------------------------------------------
  # Generate Sobol starting points
  # ----------------------------------------------------------------
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

  # ----------------------------------------------------------------
  # Default ucminfcpp control
  # ----------------------------------------------------------------
  default_ctrl <- list(
    grad     = "central",
    gradstep = c(1e-6, 1e-8),
    grtol    = 1e-4,
    xtol     = 1e-8,
    stepmax  = 5,
    maxeval  = 2000L
  )
  ctrl <- utils::modifyList(default_ctrl, control)

  # ----------------------------------------------------------------
  # Run all starts via C++ xptr backend
  # ----------------------------------------------------------------
  results <- vector("list", length(starts_list))
  for (i in seq_along(starts_list)) {
    if (verbose) {
      message(sprintf("Start %d / %d", i, length(starts_list)))
    }
    results[[i]] <- .optimize_niche_helper_cpp(
      param      = starts_list[[i]],
      env_occ    = env_occ,
      env_m      = env_m,
      control    = ctrl,
      likelihood = likelihood,
      ...
    )
  }

  # ----------------------------------------------------------------
  # Compile solutions data frame
  # ----------------------------------------------------------------
  solutions <- data.frame(
    start_id    = seq_along(results),
    loglik      = vapply(
      results, function(x) as.numeric(x$loglik), numeric(1L)
    ),
    convergence = vapply(
      results,
      function(x) as.integer(x$convergence),
      integer(1L)
    ),
    stringsAsFactors = FALSE
  )
  solutions$full_par <- lapply(results, function(x) x$theta)

  # Sort by decreasing log-likelihood
  ord       <- order(solutions$loglik, decreasing = TRUE)
  solutions <- solutions[ord, ]
  rownames(solutions) <- NULL

  best <- list(
    theta       = solutions$full_par[[1L]],
    loglik      = solutions$loglik[1L],
    convergence = solutions$convergence[1L]
  )

  # ----------------------------------------------------------------
  # Internal safeguard: validate xptr result against ucminf::ucminf
  # ----------------------------------------------------------------
  .validate_xptr_result(
    best       = best,
    env_occ    = env_occ,
    env_m      = env_m,
    likelihood = likelihood,
    ctrl       = ctrl,
    ...
  )

  if (verbose) {
    message(sprintf(
      "Best log-likelihood: %.6f (convergence = %d)",
      best$loglik,
      best$convergence
    ))
  }

  new_nicher(
    solutions  = solutions,
    best       = best,
    likelihood = likelihood,
    n_starts   = num_starts
  )
}

# -----------------------------------------------------------------------
# Internal safeguard: compare xptr result with ucminf::ucminf
# -----------------------------------------------------------------------

#' Validate ucminfcpp::ucminf_xptr result against ucminf::ucminf
#'
#' Runs a short \code{ucminf::ucminf} optimization from the best theta
#' and warns if the two log-likelihoods differ by more than \code{1e-3}.
#' This guards against pointer-safety issues in the C++ backend.
#'
#' @param best List with \code{theta} and \code{loglik}.
#' @param env_occ,env_m Data as in \code{\link{optimize_niche}}.
#' @param likelihood Character likelihood type.
#' @param ctrl Control list.
#' @param ... Forwarded to objective function.
#' @keywords internal
.validate_xptr_result <- function(best, env_occ, env_m,
                                  likelihood, ctrl, ...) {
  obj_fun <- switch(likelihood,
    unweighted    = loglik_niche_math_cpp,
    weighted      = loglik_niche_math_weighted,
    presence_only = loglik_niche_math_presence_only
  )

  ref <- tryCatch(
    {
      fn <- function(theta) {
        obj_fun(
          theta,
          env_occ = env_occ,
          env_m   = env_m,
          neg     = TRUE,
          ...
        )
      }
      validation_ctrl <- list(maxeval = 500L)
      res <- ucminf::ucminf(
        par     = best$theta,
        fn      = fn,
        control = validation_ctrl,
        hessian = FALSE
      )
      -res$value
    },
    error = function(e) NA_real_
  )

  if (is.finite(ref) && abs(best$loglik - ref) > 1e-3) {
    warning(sprintf(
      paste0(
        "ucminfcpp::ucminf_xptr and ucminf::ucminf disagree: ",
        "xptr loglik = %.6f, ucminf loglik = %.6f. ",
        "Possible pointer-safety issue."
      ),
      best$loglik, ref
    ))
  }
  invisible(NULL)
}

#' Internal helper: run ucminfcpp::ucminf_xptr for a single starting vector
#'
#' Creates a compiled C++ objective function via [create_niche_obj_ptr()] and
#' calls [ucminfcpp::ucminf_xptr()] to optimize from a single starting point.
#' This bypasses the R interpreter on every function/gradient evaluation for
#' maximum performance.
#'
#' @param param Numeric vector of starting parameters (named).
#' @param env_occ,env_m Data frames as in [optimize_niche].
#' @param control List of control parameters.
#' @param likelihood Character. One of "unweighted", "weighted", "presence_only".
#' @param ... Additional arguments (eta, m_subsample, m_kde_subsample, seed).
#' @return A list with components: theta, loglik, convergence.
#' @keywords internal
.optimize_niche_helper_cpp <- function(param, env_occ, env_m, control,
                                       likelihood, ...) {
  param_names <- names(param)
  param <- as.numeric(param)
  if (!is.null(param_names)) {
    names(param) <- param_names
  }
  if (any(!is.finite(param))) {
    stop("All starting parameters must be finite")
  }

  dots <- list(...)
  eta <- if (!is.null(dots$eta)) dots$eta else 1.0

  # Convert data frames to matrices
  env_occ_mat <- as.matrix(env_occ)
  env_m_mat <- if (!is.null(env_m)) as.matrix(env_m) else NULL

  # For the weighted likelihood, resolve subsampling indices once before
  # optimization (fixed indices ensure a smooth, consistent objective function)
  den_idx <- NULL
  kde_idx <- NULL
  precomp_w_den <- NULL

  if (likelihood == "weighted" && !is.null(env_m_mat)) {
    n_m <- nrow(env_m_mat)
    if (!is.null(dots$seed)) set.seed(dots$seed)

    pick_size <- function(x, nmax) {
      if (is.null(x)) {
        return(NULL)
      }
      if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
        stop("m_subsample/m_kde_subsample must be a single positive number.")
      }
      if (x < 1) max(1L, floor(x * nmax)) else min(nmax, as.integer(round(x)))
    }

    n_den <- pick_size(dots$m_subsample, n_m)
    den_idx <- if (!is.null(n_den) && n_den < n_m) sample.int(n_m, n_den) else NULL
    n_kde <- pick_size(dots$m_kde_subsample, n_m)
    kde_idx <- if (!is.null(n_kde) && n_kde < n_m) sample.int(n_m, n_kde) else NULL
  }

  # Create the compiled C++ objective function (external pointer).
  # Pass gradstep from control so finite-difference steps match the user's
  # configured values (defaults: c(1e-6, 1e-8), same as ucminf defaults).
  gs <- if (!is.null(control$gradstep)) control$gradstep else c(1e-6, 1e-8)
  xptr <- create_niche_obj_ptr(
    env_occ       = env_occ_mat,
    env_m         = env_m_mat,
    eta           = eta,
    likelihood    = likelihood,
    den_idx       = den_idx,
    kde_idx       = kde_idx,
    precomp_w_den = precomp_w_den,
    gradstep      = gs
  )

  # Build ucminfcpp control object by forwarding the (validated) control list.
  # This preserves existing defaults while allowing additional ucminf options.
  control_args <- control
  if (is.null(control_args$grad)) {
    control_args$grad <- "central"
  }
  if (is.null(control_args$gradstep)) {
    control_args$gradstep <- gs
  }
  if (!is.null(control_args$maxeval)) {
    control_args$maxeval <- as.integer(control_args$maxeval)
  }
  con <- do.call(ucminfcpp::ucminf_control, control_args)

  # Run optimization with error handling
  out <- tryCatch(
    {
      res <- ucminfcpp::ucminf_xptr(par = param, xptr = xptr, control = con)
      if (is.null(names(res$par)) && !is.null(names(param))) {
        names(res$par) <- names(param)
      }
      if (is.null(res$convergence)) res$convergence <- NA_integer_
      res
    },
    error = function(e) {
      list(
        par = param, value = Inf, convergence = NA_integer_,
        error = conditionMessage(e)
      )
    }
  )

  list(
    theta       = out$par,
    loglik      = -out$value,
    convergence = out$convergence
  )
}
