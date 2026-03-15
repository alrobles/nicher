#' Optimize niche model log-likelihood (math scale) with multiple starts (ucminf)
#'
#' Runs multiple optimizations using [ucminf::ucminf()] from different starting
#' points to maximize the likelihood of the ellipsoid niche model on the math
#' scale. The objective function can be either the unweighted likelihood
#' [loglik_niche_math()] (Jiménez & Soberón 2019) or the weighted-normal
#' likelihood [loglik_niche_math_weighted()] (Jiménez & Soberón 2022).
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m Data frame with environmental values from the accessibility
#'   area M (background).
#' @param num_starts Integer. Number of starting points.
#' @param quant_vec Quantiles used to define ranges for `mu` (passed to
#'   [start_theta_multiple()]).
#' @param start_method Method for generating starting points:
#'   `"sobol"` (requires \pkg{pomp}) or `"uniform"`.
#' @param likelihood Character. Either `"unweighted"` (default; uses
#'   [loglik_niche_math()]) or `"weighted"` (uses
#'   [loglik_niche_math_weighted()]).
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
#' @param ... Additional arguments passed to the chosen objective function
#'   ([loglik_niche_math()] or [loglik_niche_math_weighted()]), e.g., `eta`.
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
#' # Weighted-normal (2022) with 5 starting points
#' # (requires loglik_niche_math_weighted to be defined)
#' res2 <- optimize_niche(example_env_occ_2d, example_env_m_2d,
#'                        num_starts = 5, start_method = "uniform",
#'                        likelihood = "weighted", eta = 1)
#' res2$best
#' }
optimize_niche <- function(env_occ, env_m,
                           num_starts = 100,
                           quant_vec = c(0.1, 0.5, 0.9),
                           start_method = "sobol",
                           likelihood = c("unweighted", "weighted"),
                           control = list(),
                           verbose = FALSE,
                           ...) {
  # Match likelihood
  likelihood <- match.arg(likelihood)
  
  # Check that env_occ and env_m have the same columns
  if (!identical(sort(colnames(env_occ)), sort(colnames(env_m)))) {
    stop("env_occ and env_m must have the same variables (column names)")
  }
  
  # Verify that the chosen likelihood function exists
  if (likelihood == "weighted" && !exists("loglik_niche_math_weighted", mode = "function")) {
    stop("Weighted likelihood function 'loglik_niche_math_weighted' is not defined. ",
         "Please implement it or use likelihood = 'unweighted'.")
  }
  
  # Generate starting points (using env_m as reference for ranges)
  starts_df <- start_theta_multiple(env_data       = env_m,
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
                    unweighted = loglik_niche_math,
                    weighted   = loglik_niche_math_weighted)
  
  # Runner for a single start
  runner <- function(start_vec, id) {
    if (verbose) {
      cat(sprintf("Starting point %d\n", id))
    }
    res <- .optimize_niche_helper(
      param    = start_vec,
      obj_fun  = obj_fun,
      env_occ  = env_occ,
      env_m    = env_m,
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
#' @param obj_fun Objective function pointer: either [loglik_niche_math] or
#'   [loglik_niche_math_weighted].
#' @param env_occ,env_m Data frames as in [loglik_niche_math()].
#' @param control List of control parameters for ucminf.
#' @param ... Additional arguments passed to `obj_fun` (e.g., eta).
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