# R/autoplot_nicher.R
# S3 method on ggplot2::autoplot for a single 2-D nicher fit, plus an
# explicit convenience wrapper for comparing multiple compatible models.

#' Default ggplot2 visualization for a fitted \code{nicher} model
#'
#' \code{autoplot.nicher()} provides a default \code{ggplot2}-compatible
#' visualization of a 2-D fitted \code{nicher} model in environmental
#' space (E-space). It composes three exported layers
#' (\code{\link{geom_nicher_background}}, \code{\link{geom_nicher_occ}},
#' \code{\link{geom_nicher_ellipse}}) into a complete \code{ggplot}
#' object and returns it without printing.
#'
#' Visualization is restricted to two-dimensional models: a model with
#' \code{length(object$var_names) != 2} (or a non-2-D legacy fit) is
#' rejected with an informative error. No projections, marginalizations
#' or automatic reductions are performed.
#'
#' Fitted \code{nicher} objects do not store environmental data, so
#' \code{env_occ} (occurrence environment) and \code{env_m} (background
#' environment) must be supplied explicitly. They are interpreted
#' positionally (column 1 = x-axis variable, column 2 = y-axis variable),
#' regardless of column names.
#'
#' @param object A \code{nicher} object returned by
#'   \code{\link{optimize_niche}} with \code{length(var_names) == 2}.
#' @param env_occ Matrix or data frame of occurrence-environment values
#'   (rows = occurrences, first two columns = the two environmental
#'   variables in the order assumed by the fit).
#' @param env_m Matrix or data frame of background-environment values
#'   defining the domain of interpretation.
#' @param level Numeric vector of contour levels in \eqn{(0, 1]}. See
#'   \code{\link{geom_nicher_ellipse}}. Default \code{c(0.95, 0.5, 0.05)}.
#' @param level_type Either \code{"suitability"} (default, paper-aligned)
#'   or \code{"chisq"}.
#' @param n Integer number of points per ellipse. Default 200.
#' @param ... Currently ignored; reserved for future aesthetic overrides.
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{geom_nicher_background}},
#'   \code{\link{geom_nicher_occ}}, \code{\link{geom_nicher_ellipse}},
#'   \code{\link{nicher_compare_plot}}.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   fit <- optimize_niche(
#'     env_occ    = example_env_occ_2d,
#'     env_m      = example_env_m_2d,
#'     num_starts = 10L
#'   )
#'   autoplot(fit, env_occ = example_env_occ_2d, env_m = example_env_m_2d)
#' }
#'
#' @method autoplot nicher
#' @export autoplot.nicher
autoplot.nicher <- function(object,
                             env_occ,
                             env_m,
                             level = c(0.95, 0.5, 0.05),
                             level_type = c("suitability", "chisq"),
                             n = 200L,
                             ...) {
  .require_ggplot2("autoplot.nicher")
  .assert_nicher_2d(object, "object")
  if (missing(env_occ) || missing(env_m)) {
    stop("`env_occ` and `env_m` must be provided explicitly. ",
         "nicher objects do not store environmental data.",
         call. = FALSE)
  }
  level_type <- match.arg(level_type)

  vn <- object$var_names
  if (is.null(vn)) vn <- c("x1", "x2")

  ggplot2::ggplot() +
    geom_nicher_background(env_m,   var_names = vn) +
    geom_nicher_occ(       env_occ, var_names = vn) +
    geom_nicher_ellipse(   object,
                            level      = level,
                            level_type = level_type,
                            n          = n) +
    ggplot2::xlab(vn[1]) +
    ggplot2::ylab(vn[2])
}

#' Compare multiple compatible \code{nicher} models in 2-D E-space
#'
#' Convenience wrapper that composes background + occurrence + ellipse
#' layers for several 2-D \code{nicher} models in a single \code{ggplot}
#' object. Models are distinguished visually by mapping the model name to
#' the \code{colour} aesthetic on the ellipse layer (and on the
#' occurrence layer when \code{env_occ} is per-model).
#'
#' Two contexts are supported:
#' \itemize{
#'   \item \strong{Single species, multiple models} (e.g., presence-only
#'     vs. weighted). Pass a single matrix/data frame as \code{env_occ};
#'     it is shared across all models, drawn once with no colour
#'     grouping.
#'   \item \strong{Multiple species, multiple models}. Pass a named list
#'     of matrices/data frames as \code{env_occ}, with names identical
#'     to \code{names(models)}; each species' occurrences are drawn in
#'     the colour of its model.
#' }
#'
#' All models must share the same ordered set of environmental variables:
#' \code{identical(model$var_names, ref_var_names)} is enforced for
#' every model. No automatic reordering is performed.
#'
#' @param models Named list of \code{nicher} objects.
#' @param env_occ Either a single matrix/data frame (shared occurrence
#'   environment), or a named list of matrices/data frames whose names
#'   match \code{names(models)}.
#' @param env_m Matrix or data frame of background-environment values
#'   (shared across all models).
#' @param level Numeric vector of contour levels in \eqn{(0, 1]}.
#' @param level_type Either \code{"suitability"} (default) or \code{"chisq"}.
#' @param n Integer number of points per ellipse.
#' @param ... Reserved for future use; currently ignored.
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{autoplot.nicher}}, \code{\link{geom_nicher_ellipse}}.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2)
#'   fit_po <- optimize_niche(
#'     env_occ = example_env_occ_2d, env_m = example_env_m_2d,
#'     loglik = "presence_only", num_starts = 10L
#'   )
#'   fit_w  <- optimize_niche(
#'     env_occ = example_env_occ_2d, env_m = example_env_m_2d,
#'     loglik = "weighted", num_starts = 10L
#'   )
#'   nicher_compare_plot(
#'     models  = list(presence_only = fit_po, weighted = fit_w),
#'     env_occ = example_env_occ_2d,
#'     env_m   = example_env_m_2d
#'   )
#' }
#' @export
nicher_compare_plot <- function(models,
                                 env_occ,
                                 env_m,
                                 level = c(0.95, 0.5, 0.05),
                                 level_type = c("suitability", "chisq"),
                                 n = 200L,
                                 ...) {
  .require_ggplot2("nicher_compare_plot")
  if (!is.list(models) || !length(models)) {
    stop("`models` must be a non-empty (named) list of nicher objects.",
         call. = FALSE)
  }
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("`models` must be a named list (every element must have a name).",
         call. = FALSE)
  }
  for (i in seq_along(models)) {
    .assert_nicher_2d(models[[i]], sprintf("models[[\"%s\"]]",
                                            names(models)[i]))
  }
  ref_vn <- .assert_var_names_compatible(models)
  level_type <- match.arg(level_type)

  vn <- if (is.null(ref_vn)) c("x1", "x2") else ref_vn

  # Build occurrence layer(s).
  if (is.list(env_occ) && !is.data.frame(env_occ)) {
    if (is.null(names(env_occ)) ||
        !setequal(names(env_occ), names(models))) {
      stop("Names of `env_occ` must match `names(models)` exactly.",
           call. = FALSE)
    }
    occ_df <- do.call(rbind, lapply(names(models), function(nm) {
      d <- .coerce_xy(env_occ[[nm]], vn, name = sprintf("env_occ[[\"%s\"]]",
                                                          nm))
      d$model <- nm
      d
    }))
    occ_nm <- names(occ_df)[1:2]
    occ_layer <- ggplot2::layer(
      data        = occ_df,
      mapping     = ggplot2::aes(x = .data[[ occ_nm[1] ]],
                                  y = .data[[ occ_nm[2] ]],
                                  colour = .data[["model"]]),
      geom        = "point",
      stat        = "identity",
      position    = "identity",
      inherit.aes = FALSE,
      params      = list(size = 1.2)
    )
  } else {
    # Shared occurrences: single layer, no colour grouping.
    occ_layer <- geom_nicher_occ(env_occ, var_names = vn)
  }

  # Build one ellipse-data table for ALL models, with model + level
  # columns so a single geom_path with colour = model and group =
  # interaction(model, level) does the work.
  ell_parts <- lapply(names(models), function(nm) {
    ms <- .recover_mu_sigma(models[[nm]])
    d  <- .ellipse_paths(ms$mu, ms$Sigma, level, level_type, n,
                         var_names = vn)
    if (!nrow(d)) return(NULL)
    d$model <- nm
    d
  })
  ell_parts <- ell_parts[!vapply(ell_parts, is.null, logical(1L))]
  if (!length(ell_parts)) {
    stop("No ellipse paths could be constructed for any model ",
         "(all Sigma matrices were rank-deficient?).", call. = FALSE)
  }
  ell_df <- do.call(rbind, ell_parts)
  ell_df$.grp <- paste(ell_df$model, ell_df$level, sep = "_")
  ell_nm <- names(ell_df)[1:2]

  lw <- .linewidth_param(0.6)
  ell_layer <- ggplot2::layer(
    data        = ell_df,
    mapping     = ggplot2::aes(x = .data[[ ell_nm[1] ]],
                                y = .data[[ ell_nm[2] ]],
                                colour = .data[["model"]],
                                group  = .data[[".grp"]]),
    geom        = "path",
    stat        = "identity",
    position    = "identity",
    inherit.aes = FALSE,
    params      = lw
  )

  ggplot2::ggplot() +
    geom_nicher_background(env_m, var_names = vn) +
    occ_layer +
    ell_layer +
    ggplot2::xlab(vn[1]) +
    ggplot2::ylab(vn[2]) +
    ggplot2::labs(colour = "model")
}
