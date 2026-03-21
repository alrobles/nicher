#' Optimize niche model using a compiled XPtr backend
#'
#' Runs ucminfcpp::ucminf_xptr() over a compiled C++ objective function
#' created by create_niche_obj_ptr(). Supports single-start and multi-start
#' optimization. Ensures that all starting vectors are numeric doubles and finite.
#'
#' @param start Numeric vector (single start) or list of numeric vectors (multi-start).
#' @param xptr External pointer created by create_niche_obj_ptr().
#' @param control List of control parameters for ucminfcpp::ucminf_control().
#' @param multi_start Logical; TRUE if multiple starts are supplied.
#'
#' @return A list with:
#'   \itemize{
#'     \item par — best parameter vector
#'     \item value — negative log-likelihood
#'     \item conv — convergence code
#'     \item all_results — (multi-start only) table of all runs
#'   }
#'
#' @export
optimize_niche_xptr <- function(start = NULL,
                                xptr,
                                control = ucminfcpp::ucminf_control(
                                  grad = "central",
                                  gradstep = c(1e-6, 1e-8),
                                  maxeval = 200
                                ),
                                multi_start = FALSE) {

  # ----------------------------------------------------------
  # Sanity check: start provided?
  # ----------------------------------------------------------
  if (is.null(start))
    stop("Start vector must be provided.")

  # ----------------------------------------------------------
  # Helper to sanitize a single starting vector
  # ----------------------------------------------------------
  sanitize_start <- function(s) {
    if (is.null(s))
      stop("NULL start received.")

    # Convert to numeric
    s <- as.numeric(s)

    # Must be finite doubles
    if (!is.numeric(s))
      stop("Start contains non-numeric values.")

    if (any(!is.finite(s)))
      stop("Start contains non-finite (NA/Inf) values.")

    # Drop attributes
    attributes(s) <- NULL

    storage.mode(s) <- "double"
    return(s)
  }

  # ==========================================================
  # SINGLE-START CASE
  # ==========================================================
  if (!multi_start) {

    s <- sanitize_start(start)

    out <- tryCatch(
      ucminfcpp:::ucminf_xptr(
        par = s,
        xptr = xptr,
        control = control
      ),
      error = function(e) {
        stop("Optimization failed in single-start mode: ", e$message)
      }
    )

    # Ensure convergence code exists
    if (is.null(out$convergence))
      out$convergence <- NA_integer_

    return(list(
      par  = out$par,
      value = out$value,
      conv = out$convergence
    ))
  }

  # ==========================================================
  # MULTI-START CASE
  # ==========================================================
  if (!is.list(start))
    stop("multi_start = TRUE requires a list of numeric start vectors.")

  all_results <- vector("list", length(start))

  for (i in seq_along(start)) {

    s <- tryCatch(
      sanitize_start(start[[i]]),
      error = function(e) {
        # register invalid start cleanly
        all_results[[i]] <<- list(start_id = i, value = NA, conv = NA)
        return(NULL)
      }
    )

    if (is.null(s))
      next

    res <- tryCatch(
      ucminfcpp:::ucminf_xptr(
        par = s,
        xptr = xptr,
        control = control
      ),
      error = function(e) {
        all_results[[i]] <<- list(start_id = i, value = NA, conv = NA)
        return(NULL)
      }
    )

    if (is.null(res)) next

    if (is.null(res$convergence))
      res$convergence <- NA_integer_

    all_results[[i]] <- list(
      start_id = i,
      value    = res$value,
      conv     = res$convergence,
      par      = res$par
    )
  }

  # Consolidate results
  df <- do.call(rbind, lapply(all_results, function(r) {
    data.frame(
      start_id = r$start_id,
      value    = r$value,
      conv     = r$conv
    )
  }))

  # Best finite result
  idx <- which(is.finite(df$value))
  if (length(idx) == 0)
    stop("No valid start produced a finite objective value.")

  best <- idx[which.min(df$value)]

  return(list(
    par         = all_results[[best]]$par,
    value       = all_results[[best]]$value,
    conv        = all_results[[best]]$conv,
    all_results = df
  ))
}
