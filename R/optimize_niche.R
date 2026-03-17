#' Optimize niche model log-likelihood (math scale) with multiple starts (ucminf)
#'
#' Runs multiple optimizations using [ucminf::ucminf()] from different starting
#' points to maximize the likelihood of the ellipsoid niche model on the math
#' scale. The objective function can be:
#' \itemize{
#'   \item `"unweighted"`: the basic ellipsoid model with background correction
#'         (Jiménez & Soberón 2019) – uses [loglik_niche_math_cpp()].
#'   \item `"weighted"`: the weighted‑normal model with integrated KDE
#'         (Jiménez & Soberón 2022) – uses [loglik_niche_math_weighted()].
#'   \item `"presence_only"`: the simplest model using only presence points
#'         (no background correction) – uses [loglik_niche_math_presence_only()].
#' }
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m Data frame with environmental values from the accessibility
#'   area M (background). Ignored if `likelihood = "presence_only"`.
#' @param num_starts Integer. Number of starting points.
#' @param quant_vec Quantiles used to define ranges for `mu` (passed to
#'   [start_theta_multiple()]).
#' @param start_method Method for generating starting points:
#'   `"sobol"` (requires \pkg{pomp}) or `"uniform"`.
#' @param likelihood Character. One of `"unweighted"` (default), `"weighted"`,
#'   or `"presence_only"`.
#' @param control List of control parameters passed to [ucminf::ucminf()].
#'   Default values (if not supplied) are:
#'   \describe{
#'     \item{grad}{"central"}
#'     \item{gradstep}{c(1e-6, 1e-8)}
#'     \item{grtol}{1e-4}
#'     \item{xtol}{1e-8}
#'     \item{stepmax}{5}
#'     \item{maxeval}{2000}
#'   }
#' @param verbose Logical. If `TRUE`, print progress messages.
#' @param ... Additional arguments passed to the chosen objective function,
#'   e.g., `eta` for the LKJ prior, or subsampling parameters for the weighted
#'   model (`m_subsample`, `m_kde_subsample`, `seed`).
#'
#' @return A list with two components:
#'   \item{solutions}{A data frame with columns:
#'     \itemize{
#'       \item `start_id`: integer identifier of the starting point.
#'       \item `loglik`: the maximized log-likelihood (not negative).
#'       \item `convergence`: convergence code from `ucminf` (1 or 2 mean success).
#'       \item `full_par`: list column with the full parameter vector at the optimum.
#'     }}
#'   \item{best}{A list with the best solution: `theta`, `loglik`, `convergence`.}
#'
#' @export
#' @examples
#' \dontrun{
#' # Unweighted (2019) with 5 starting points
#' res1 <- optimize_niche(example_env_occ_2d, example_env_m_2d,
#'                        num_starts = 5, start_method = "uniform",
#'                        likelihood = "unweighted", eta = 1)
#' res1$best
#'
#' # Weighted (2022) with 5 starting points
#' res2 <- optimize_niche(example_env_occ_2d, example_env_m_2d,
#'                        num_starts = 5, start_method = "uniform",
#'                        likelihood = "weighted", eta = 1,
#'                        m_subsample = 2000, m_kde_subsample = 5000, seed = 123)
#' res2$best
#'
#' # Presence‑only (ultra‑light) with 5 starting points
#' res3 <- optimize_niche(example_env_occ_2d, env_m = NULL,
#'                        num_starts = 5, start_method = "uniform",
#'                        likelihood = "presence_only", eta = 1)
#' res3$best
#' }
optimize_niche <- function(env_occ, env_m,
                           num_starts = 100,
                           quant_vec = c(0.1, 0.5, 0.9),
                           start_method = "sobol",
                           likelihood = c("unweighted", "weighted", "presence_only"),
                           control = list(),
                           verbose = FALSE,
                           ...) {
  # Match likelihood
  likelihood <- match.arg(likelihood)
  
  # Validate inputs based on likelihood
  if (likelihood != "presence_only") {
    if (missing(env_m) || is.null(env_m)) stop("env_m must be provided for likelihood '", likelihood, "'")
    if (!identical(sort(colnames(env_occ)), sort(colnames(env_m)))) {
      stop("env_occ and env_m must have the same variables (column names)")
    }
  }
  
  # Determine which data to use for generating starting ranges
  range_data <- if (likelihood == "presence_only") env_occ else env_m
  
  # Generate starting points (using the appropriate data for ranges)
  starts_df <- start_theta_multiple(env_data  = range_data,
                                    num_starts = num_starts,
                                    quant_vec  = quant_vec,
                                    method     = start_method)
  
  # Convert each row to a named numeric vector
  starts_list <- split(starts_df, seq_len(nrow(starts_df)))
  starts_list <- lapply(starts_list, function(x) {
    v <- as.numeric(x)
    names(v) <- colnames(starts_df)
    v
  })
  
  # Default control for ucminf
  default_ctrl <- list(
    grad     = "central",
    gradstep = c(1e-6, 1e-8),
    grtol    = 1e-4,
    xtol     = 1e-8,
    stepmax  = 5,
    maxeval  = 2000
  )
  ctrl <- utils::modifyList(default_ctrl, control)
  
  # Choose objective function pointer
  obj_fun <- switch(likelihood,
                    unweighted    = loglik_niche_math_cpp,
                    weighted      = loglik_niche_math_weighted,
                    presence_only = loglik_niche_math_presence_only)
  
  # Runner for a single start
  runner <- function(start_vec, id) {
    if (verbose) {
      cat(sprintf("Starting point %d\n", id))
    }
    res <- .optimize_niche_helper(
      param    = start_vec,
      obj_fun  = obj_fun,
      env_occ  = env_occ,
      env_m    = env_m,          # will be ignored inside presence_only if not used
      control  = ctrl,
      ...
    )
    list(
      theta       = res$par,
      loglik      = -res$value,
      convergence = res$convergence
    )
  }
  
  # Run sequentially
  results <- vector("list", length(starts_list))
  for (i in seq_along(starts_list)) {
    results[[i]] <- runner(starts_list[[i]], i)
  }
  
  # Compile results
  solutions <- data.frame(
    start_id   = seq_along(results),
    loglik     = sapply(results, function(x) x$loglik),
    convergence = sapply(results, function(x) x$convergence),
    stringsAsFactors = FALSE
  )
  solutions$full_par <- lapply(results, function(x) x$theta)
  
  # Sort by decreasing log-likelihood
  ord <- order(solutions$loglik, decreasing = TRUE)
  solutions <- solutions[ord, ]
  rownames(solutions) <- NULL
  
  best <- list(
    theta       = solutions$full_par[[1]],
    loglik      = solutions$loglik[1],
    convergence = solutions$convergence[1]
  )
  
  if (verbose) {
    cat(sprintf("Best log-likelihood: %.6f (convergence = %d)\n",
                best$loglik, best$convergence))
  }
  
  list(solutions = solutions, best = best)
}

#' Internal helper: run ucminf for a single starting vector
#'
#' @param param Numeric vector of starting parameters (named).
#' @param obj_fun Objective function pointer: either [loglik_niche_math_cpp],
#'   [loglik_niche_math_weighted], or [loglik_niche_math_presence_only].
#' @param env_occ,env_m Data frames as in [optimize_niche].
#' @param control List of control parameters for ucminf.
#' @param ... Additional arguments passed to `obj_fun` (e.g., eta, subsampling).
#' @return A list with components: par, value, convergence, invhessian.lt (if computed).
#' @keywords internal
.optimize_niche_helper <- function(param, obj_fun, env_occ, env_m, control, ...) {
  # Ensure param is numeric and finite
  param <- as.numeric(param)
  if (any(!is.finite(param))) {
    stop("All starting parameters must be finite")
  }
  
  # Objective function (returns negative log-likelihood)
  fn <- function(theta) {
    obj_fun(theta,
            env_occ = env_occ,
            env_m   = env_m,
            neg     = TRUE,
            ...)
  }
  
  # Run ucminf with error handling
  out <- tryCatch(
    {
      res <- ucminf::ucminf(par = param,
                            fn  = fn,
                            control = control,
                            hessian = FALSE)  # Hessian not needed for optimization
      # Restore names if ucminf drops them
      if (is.null(names(res$par)) && !is.null(names(param))) {
        names(res$par) <- names(param)
      }
      if (is.null(res$convergence)) res$convergence <- NA_integer_
      res
    },
    error = function(e) {
      list(par = param,
           value = Inf,
           convergence = NA_integer_,
           error = conditionMessage(e))
    }
  )
  
  list(
    par = out$par,
    value = out$value,
    convergence = out$convergence,
    invhessian.lt = if (!is.null(out$invhessian.lt)) out$invhessian.lt else NULL
  )
}