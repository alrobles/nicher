#' Benchmark cpp vs r backends for optimize_niche
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("experimental")}
#'
#' Runs \code{\link{optimize_niche}} from the same Sobol starts under both
#' \code{backend = "cpp"} and \code{backend = "r"} and returns timing /
#' convergence statistics for direct comparison. This function exists ONLY
#' to enable side-by-side benchmarking of the new C++ XPtr backend against
#' the legacy R-objective backend, and will be removed in the next minor
#' release alongside \code{backend = "r"}.
#'
#' @param env_occ Data frame of presence-point environmental values.
#' @param env_m Data frame of background environmental values
#'   (required for \code{likelihood = "weighted"}).
#' @param num_starts Integer. Number of Sobol starting points
#'   (default \code{10L}; small to keep the legacy backend fast).
#' @param likelihood One of \code{"weighted"} or \code{"presence_only"}.
#' @param seed Integer seed (passed to both runs for identical sampling).
#' @param control Named list of control parameters
#'   (forwarded to \code{optimize_niche()}).
#' @param ... Additional arguments forwarded to \code{optimize_niche()}.
#'
#' @return A data frame with one row per backend, columns:
#'   \describe{
#'     \item{\code{backend}}{"cpp" or "r"}
#'     \item{\code{best_loglik}}{Log-likelihood of the best converged start.}
#'     \item{\code{best_convergence}}{ucminf convergence code of the best.}
#'     \item{\code{n_converged}}{Count of starts with convergence \%in\% \{1, 2\}.}
#'     \item{\code{wall_time_sec}}{Total wall-clock seconds for the run.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' bench <- benchmark_optimize_niche(
#'   env_occ    = example_env_occ_2d,
#'   env_m      = example_env_m_2d,
#'   num_starts = 5L,
#'   likelihood = "weighted",
#'   seed       = 1L
#' )
#' print(bench)
#' }
benchmark_optimize_niche <- function(env_occ,
                                     env_m      = NULL,
                                     num_starts = 10L,
                                     likelihood = c("weighted",
                                                    "presence_only"),
                                     seed       = NULL,
                                     control    = list(),
                                     ...) {
  likelihood <- match.arg(likelihood)

  run_one <- function(backend) {
    t0 <- proc.time()[["elapsed"]]
    res <- suppressWarnings(
      optimize_niche(
        env_occ    = env_occ,
        env_m      = env_m,
        num_starts = num_starts,
        likelihood = likelihood,
        backend    = backend,
        seed       = seed,
        control    = control,
        verbose    = FALSE,
        ...
      )
    )
    t1 <- proc.time()[["elapsed"]]
    list(
      backend          = backend,
      best_loglik      = res$best$loglik,
      best_convergence = res$best$convergence,
      n_converged      = sum(res$solutions$convergence %in% c(1L, 2L)),
      wall_time_sec    = t1 - t0
    )
  }

  out <- lapply(c("cpp", "r"), run_one)
  do.call(rbind, lapply(out, as.data.frame))
}
