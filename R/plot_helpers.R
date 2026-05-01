# R/plot_helpers.R
# Private helpers shared by autoplot.nicher / nicher_compare_plot /
# geom_nicher_*. None are exported.

# ---------------------------------------------------------------------------
# 20-node standard Gauss-Laguerre quadrature (nodes / weights from
# numpy.polynomial.laguerre.laggauss(20)). Used by
# geom_nicher_isosuitability() to evaluate the NCST suitability integral
# on a 2-D plotting grid. The optimization-side C++ kernel uses 32 nodes
# for tighter accuracy on higher-dim fits; 20 nodes is more than enough
# for the (smooth) plotting integrand.
# ---------------------------------------------------------------------------
.gl_quad_20 <- function() {
  nodes <- c(
    7.0539889691988739e-02, 3.7212681800161157e-01, 9.1658210248327376e-01,
    1.7073065310283435e+00, 2.7491992553094322e+00, 4.0489253138508881e+00,
    5.6151749708616174e+00, 7.4590174536710627e+00, 9.5943928695810978e+00,
    1.2038802546964316e+01, 1.4814293442630740e+01, 1.7948895520519375e+01,
    2.1478788240285009e+01, 2.5451702793186904e+01, 2.9932554631700611e+01,
    3.5013434240479000e+01, 4.0833057056728570e+01, 4.7619994047346502e+01,
    5.5810795750063896e+01, 6.6524416525615749e+01
  )
  weights <- c(
    1.6874680185113369e-01, 2.9125436200606059e-01, 2.6668610286699662e-01,
    1.6600245326950186e-01, 7.4826064668790743e-02, 2.4964417309282588e-02,
    6.2025508445721095e-03, 1.1449623864768774e-03, 1.5574177302780828e-04,
    1.5401440865224536e-05, 1.0864863665179549e-06, 5.3301209095566101e-08,
    1.7579811790505475e-09, 3.7255024025121629e-11, 4.7675292515780503e-13,
    3.3728442433624615e-15, 1.1550143395003684e-17, 1.5395221405823035e-20,
    5.2864427255689281e-24, 1.6564566124990854e-28
  )
  list(nodes = nodes, weights = weights)
}

# ---------------------------------------------------------------------------
# Recover (mu, Sigma) from a fitted nicher object.
#
# Mirrors the canonical recipe used by predict.nicher() at
# R/predict_nicher.R:99-124. Lifting it here lets the plotting code recover
# the same niche geometry that predict.nicher() uses, by construction.
#
# Returns a list with: mu (length p), Sigma (p x p), p (integer),
# var_names (character or NULL), eta (numeric scalar).
# Errors if `theta` is malformed.
# ---------------------------------------------------------------------------
.recover_mu_sigma <- function(object) {
  theta <- object$best$theta
  if (!is.numeric(theta) || length(theta) < 2L) {
    stop("`object$best$theta` is missing or malformed.", call. = FALSE)
  }
  k <- length(theta)
  lik <- if (!is.null(object$likelihood)) object$likelihood else ""
  is_skew   <- lik %in% c("skew_normal", "skew_normal_weighted",
                          "skew_t", "skew_t_weighted")
  is_skew_t <- lik %in% c("skew_t", "skew_t_weighted")

  # Gaussian layout : k = 2p + p(p-1)/2          =>  p = (-3 + sqrt(9 + 8k))/2
  # Skew-normal     : k = 3p + p(p-1)/2          =>  p = (-5 + sqrt(25 + 8k))/2
  # Skew-t          : k = 3p + p(p-1)/2 + 1      =>  same as skew-normal but
  #                                                  on (k-1) instead of k
  if (is_skew_t) {
    p_dbl <- (-5 + sqrt(25 + 8 * (k - 1L))) / 2
  } else if (is_skew) {
    p_dbl <- (-5 + sqrt(25 + 8 * k)) / 2
  } else {
    p_dbl <- (-3 + sqrt(9 + 8 * k)) / 2
  }
  p <- as.integer(round(p_dbl))
  expected_k <- if (is_skew_t)
                  3L * p + p * (p - 1L) / 2L + 1L
                else if (is_skew)
                  3L * p + p * (p - 1L) / 2L
                else
                  2L * p + p * (p - 1L) / 2L
  if (abs(p_dbl - p) > 1e-8 || p < 1L || k != expected_k) {
    stop("Cannot infer p from length(object$best$theta) = ", k,
         " (likelihood = '", lik, "').", call. = FALSE)
  }

  n_v <- p * (p - 1L) / 2L
  mu    <- theta[seq_len(p)]
  sigma <- exp(theta[(p + 1L):(2L * p)])
  v     <- if (n_v > 0L) theta[(2L * p + 1L):(2L * p + n_v)] else numeric(0)
  alpha <- if (is_skew) theta[(2L * p + n_v + 1L):(3L * p + n_v)]
           else         rep(0.0, p)
  log_r <- if (is_skew_t) theta[3L * p + n_v + 1L]
           else           NA_real_

  # Same eta-fallback as predict.nicher() for legacy fits without `eta`.
  eta_fit <- if (is.null(object$eta)) 1.0 else object$eta
  L_corr  <- cvine_cholesky(v, d = p, eta = eta_fit)
  L_cov   <- diag(sigma, p) %*% L_corr
  Sigma   <- tcrossprod(L_cov)

  list(
    mu        = mu,
    Sigma     = Sigma,
    alpha     = alpha,
    log_r     = log_r,
    is_skew   = is_skew,
    is_skew_t = is_skew_t,
    p         = p,
    var_names = object$var_names,
    eta       = eta_fit
  )
}

# ---------------------------------------------------------------------------
# Hard 2-D guard. var_names is the source of truth when present; for legacy
# fits without var_names we fall back to recovered `p`.
# ---------------------------------------------------------------------------
.assert_nicher_2d <- function(object, name = "object") {
  if (!inherits(object, "nicher")) {
    stop(sprintf("`%s` must be a nicher object returned by optimize_niche().",
                 name), call. = FALSE)
  }
  p <- if (!is.null(object$var_names)) length(object$var_names)
       else .recover_mu_sigma(object)$p
  if (!identical(p, 2L)) {
    stop(sprintf(paste0(
      "`%s` must be a 2-D nicher fit (got p = %d). Plotting in E-space ",
      "is restricted to two environmental variables; no projections, ",
      "marginalizations, or automatic reductions are supported."),
      name, p), call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Compatibility check across a list of nicher models.
# ref_var_names defaults to models[[1]]$var_names; errors out the offender.
# ---------------------------------------------------------------------------
.assert_var_names_compatible <- function(models, ref_var_names = NULL) {
  if (!is.list(models) || !length(models)) {
    stop("`models` must be a non-empty list of nicher objects.",
         call. = FALSE)
  }
  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "nicher")) {
      stop(sprintf(
        "models[[%d]] (%s) is not a nicher object.",
        i, names(models)[i] %||% "<unnamed>"), call. = FALSE)
    }
  }
  if (is.null(ref_var_names)) ref_var_names <- models[[1L]]$var_names
  if (is.null(ref_var_names)) {
    stop("Reference `var_names` is missing. Pass `ref_var_names` ",
         "explicitly when models do not carry var_names.", call. = FALSE)
  }
  for (i in seq_along(models)) {
    vn <- models[[i]]$var_names
    if (!identical(vn, ref_var_names)) {
      stop(sprintf(paste0(
        "Models are not compatible: models[[%d]] (%s) has var_names = %s, ",
        "expected %s. No automatic reordering is performed."),
        i, names(models)[i] %||% "<unnamed>",
        deparse(vn), deparse(ref_var_names)),
        call. = FALSE)
    }
  }
  invisible(ref_var_names)
}

# ---------------------------------------------------------------------------
# Coerce env data (matrix or data.frame) to a 2-col data.frame using
# positional indexing. Column names default to var_names when supplied so
# axis labels and aesthetics line up with the model.
# ---------------------------------------------------------------------------
.coerce_xy <- function(x, var_names = NULL, name = "env") {
  if (is.matrix(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a matrix or data frame.", name),
         call. = FALSE)
  }
  if (ncol(x) < 2L) {
    stop(sprintf(paste0("`%s` must have at least 2 columns ",
                         "(got ncol = %d). Plotting is positional: ",
                         "the first two columns are used as the ",
                         "two environmental variables."),
                  name, ncol(x)), call. = FALSE)
  }
  out <- data.frame(
    x[[1L]], x[[2L]],
    stringsAsFactors = FALSE
  )
  vn <- if (!is.null(var_names) && length(var_names) >= 2L) {
    as.character(var_names[1:2])
  } else {
    c("x1", "x2")
  }
  names(out) <- vn
  out
}

# ---------------------------------------------------------------------------
# Single-level ellipse path in E-space.
#   level_type = "suitability"  =>  S(x) = level   =>  d^2 = -2 log(level)
#   level_type = "chisq"        =>  d^2 = qchisq(level, df = 2)
# Returns a closed (n+1) x 2 path: last row == first row.
# ---------------------------------------------------------------------------
.ellipse_path <- function(mu, Sigma, level,
                          level_type = c("suitability", "chisq"),
                          n = 200L) {
  level_type <- match.arg(level_type)
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level > 1) {
    stop("`level` must be a single number in (0, 1].", call. = FALSE)
  }
  if (level_type == "suitability" && level == 1) {
    # S = 1 is a single point at mu; degenerate ellipse with r = 0.
    r2 <- 0
  } else {
    r2 <- switch(level_type,
                  suitability = -2 * log(level),
                  chisq       = stats::qchisq(level, df = 2L))
  }

  # Σ = L Lᵀ with L lower-triangular. chol() returns U upper-triangular
  # so we transpose. tryCatch() so we degrade gracefully on a degenerate
  # Σ (e.g., a fit that pinned σ_k → 0) instead of erroring out the whole
  # ggplot build.
  L <- tryCatch(t(chol(Sigma)),
                error = function(e) NULL)
  if (is.null(L)) {
    warning("Sigma is rank-deficient; returning empty ellipse path.",
            call. = FALSE)
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  t  <- seq(0, 2 * pi, length.out = n + 1L)
  uv <- cbind(cos(t), sin(t)) * sqrt(r2)
  pts <- t(L %*% t(uv))
  data.frame(
    x = mu[1] + pts[, 1L],
    y = mu[2] + pts[, 2L]
  )
}

# ---------------------------------------------------------------------------
# Stack of ellipse paths for a vector `level`. Adds a `level` column for
# downstream group/colour/linetype mapping. var_names are used as the
# `x` / `y` column names when supplied.
# ---------------------------------------------------------------------------
.ellipse_paths <- function(mu, Sigma, level,
                            level_type = c("suitability", "chisq"),
                            n = 200L,
                            var_names = NULL) {
  level_type <- match.arg(level_type)
  if (!is.numeric(level) || !length(level) ||
      any(!is.finite(level)) || any(level <= 0) || any(level > 1)) {
    stop("`level` must be a numeric vector with values in (0, 1].",
         call. = FALSE)
  }
  vn <- if (!is.null(var_names) && length(var_names) >= 2L) {
    as.character(var_names[1:2])
  } else {
    c("x1", "x2")
  }

  parts <- lapply(level, function(lv) {
    pth <- .ellipse_path(mu, Sigma, lv, level_type, n)
    if (!nrow(pth)) return(NULL)
    out <- data.frame(pth$x, pth$y, level = lv)
    names(out)[1:2] <- vn
    out
  })
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (!length(parts)) {
    out <- data.frame(numeric(0), numeric(0), level = numeric(0))
    names(out)[1:2] <- vn
    return(out)
  }
  do.call(rbind, parts)
}

# ---------------------------------------------------------------------------
# Tiny null-coalesce. Avoids an rlang dependency for one operator.
# ---------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# ggplot2 renamed `size` -> `linewidth` for line geoms in 3.4.0. Returns
# the right name based on the *installed* ggplot2 so we work cleanly on
# both old (3.3.x) and new (3.4+) installs.
# ---------------------------------------------------------------------------
.linewidth_param <- function(value) {
  has_linewidth <- tryCatch(
    utils::packageVersion("ggplot2") >= "3.4.0",
    error = function(e) FALSE
  )
  out <- list()
  if (has_linewidth) out$linewidth <- value else out$size <- value
  out
}

.require_ggplot2 <- function(fn) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(sprintf("`%s()` requires the 'ggplot2' package. ",
                  fn),
         "Install it with install.packages(\"ggplot2\").",
         call. = FALSE)
  }
  invisible(TRUE)
}
