#' Optimize niche model log-likelihood from multiple starting points
#'
#' @param env_occ Data frame with environmental values at presence points.
#' @param env_m Data frame with environmental values from M area.
#' @param num_starts Integer. Number of starting points.
#' @param quant_vec Quantiles for generating ranges.
#' @param start_method Starting point generation method ("sobol" or "uniform").
#' @param control List of control parameters passed to optim (e.g., list(maxit = 1000)).
#' @param ... Additional arguments passed to loglik_niche_math (e.g., eta).
#'
#' @return A list with components:
#'   \item{solutions}{Data frame with columns start_id, loglik, convergence, and list-column full_par (theta vector).}
#'   \item{best}{List with best theta, loglik, and convergence.}
#' @export
optimize_niche <- function(env_occ, env_m,
                           num_starts = 100,
                           quant_vec = c(0.1, 0.5, 0.9),
                           start_method = "sobol",
                           #parallel = FALSE,
                           control = list(),
                           ...) {
  # Generate starting points
  starts_df <- start_theta_multiple(env_m, num_starts, quant_vec, method = start_method)
  starts_list <- split(starts_df, seq_len(nrow(starts_df)))
  starts_list <- lapply(starts_list, function(x) as.numeric(unlist(x)))
  
  # Objective function (negative log-likelihood) with fixed data
  obj_fun <- function(theta) {
    loglik_niche_math(theta, env_occ = env_occ, env_m = env_m, neg = TRUE, ...)
  }
  
  # Default control for optim
  default_ctrl <- list(maxit = 5000, trace = 0)
  ctrl <- utils::modifyList(default_ctrl, control)
  
  # Run one optimization
  run_one <- function(start) {
    res <- tryCatch(
      stats::optim(par = start, fn = obj_fun, method = "BFGS", control = ctrl),
      error = function(e) list(par = start, value = Inf, convergence = NA)
    )
    list(
      theta = res$par,
      loglik = -res$value,
      convergence = if (!is.null(res$convergence)) res$convergence else NA
    )
  }
  
  # if (parallel) {
  #   if (!requireNamespace("furrr", quietly = TRUE)) {
  #     stop("Package 'furrr' required for parallel execution. Please install it.")
  #   }
  #   future::plan("multisession", workers = future::availableCores())
  #   on.exit(future::plan("sequential"), add = TRUE)
  #   results <- furrr::future_map(starts_list, run_one, .progress = TRUE)
  # } else {
  #   results <- lapply(starts_list, run_one)
  # }
  results <- lapply(starts_list, run_one)
  # Compile results
  solutions <- data.frame(
    start_id = seq_along(results),
    loglik = sapply(results, function(x) x$loglik),
    convergence = sapply(results, function(x) x$convergence),
    stringsAsFactors = FALSE
  )
  solutions$full_par <- lapply(results, function(x) x$theta)
  
  # Sort by decreasing log-likelihood
  ord <- order(solutions$loglik, decreasing = TRUE)
  solutions <- solutions[ord, ]
  rownames(solutions) <- NULL
  
  best <- list(
    theta = solutions$full_par[[1]],
    loglik = solutions$loglik[1],
    convergence = solutions$convergence[1]
  )
  
  list(solutions = solutions, best = best)
}