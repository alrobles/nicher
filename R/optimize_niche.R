#' Optimize niche model log-likelihood with multi-start Sobol design
#'
#' Runs multi-start optimization over a Sobol low-discrepancy sequence
#' (via \pkg{pomp}) of starting points covering the parameter space implied
#' by \code{env_occ} and \code{breadth}.
#'
#' The supported likelihood models are:
#' \itemize{
#'   \item \code{"weighted"} (default): paper-faithful weighted-normal model.
#'         Implements Eq. 5 of Jiménez & Soberón (2022, Ecological
#'         Modelling 438:109982) -- pure ML estimation of a weighted normal
#'         density on M -- plus a weakly-informative ridge prior on
#'         \eqn{\log \sigma} that prevents the well-known Patil & Ord
#'         (1976) \eqn{\sigma \to \infty} drift in the un-regularised
#'         MLE. The ridge strength is controlled by
#'         \code{prior_log_sigma_lambda} (default \code{1.0}, weak); set
#'         it to \code{0} for pure paper Eq. 5 (not recommended on
#'         multimodal M).
#'   \item \code{"kde_bias_corrected"}: legacy KDE-bias-corrected
#'         weighted-normal formula. The KDE of M down-weights
#'         occurrences sitting in densely sampled background (a
#'         sampling-bias correction). Empirically stable, no
#'         \eqn{\sigma \to \infty} drift, but does not match the
#'         published Eq. 5.
#'   \item \code{"presence_only"}: model using only presence points,
#'         no background correction. Unimodal likelihood, useful as a
#'         sanity check or to seed the weighted multistart (see
#'         \code{warm_start}).
#'   \item \code{"skew_normal"}: presence-only fit of a multivariate
#'         skew-normal niche
#'         (Azzalini & Capitanio 1999, J. R. Stat. Soc. Ser. B 61(3):
#'         579-602):
#'         \deqn{S(x) \propto \phi_p(x - \mu; \Sigma) \, \Phi\left(
#'           \sum_k \alpha_k \, (x_k - \mu_k) / \sigma_k \right).}
#'         Adds a length-\code{p} skewness vector \eqn{\alpha} to the
#'         existing \eqn{(\mu, \Sigma)} parameters. \eqn{\alpha = 0}
#'         recovers the symmetric Gaussian niche exactly.
#'         Sobol-start machinery samples \eqn{\alpha_k} in
#'         \eqn{[-3, 3]} (Azzalini & Capitanio 1999, Sec. 5).
#'   \item \code{"skew_normal_weighted"}: paper Eq. 5 with the SN
#'         density in place of the Gaussian, plus the same ridge
#'         prior on \eqn{\log \sigma} used by \code{"weighted"}.
#'         Defaults inherit from \code{"weighted"}.
#'   \item \code{"skew_t"}: presence-only fit of the multivariate
#'         non-central skew-t (NCST) density (Branco & Dey 2001, J.
#'         Multivariate Analysis 79(1):99-113):
#'         \deqn{T = \mu + X \sqrt{r/Y}, \quad X \sim SN_p(0, \Sigma, \alpha),
#'                                                Y \sim \chi^2_r.}
#'         Adds a single degrees-of-freedom parameter \code{log_r}
#'         on top of the skew-normal layout. Smaller \eqn{r} = heavier
#'         tails; \eqn{r \to \infty} recovers \code{"skew_normal"}.
#'         The marginal density has no closed form, so it is computed
#'         by 32-node Gauss-Laguerre quadrature on the \eqn{\chi^2_r}
#'         mixing variable. Sobol-start machinery samples \code{log_r}
#'         in \eqn{[\log 2, \log 100]}.
#'   \item \code{"skew_t_weighted"}: paper Eq. 5 with the NCST
#'         density in place of the Gaussian, plus the same ridge prior
#'         on \eqn{\log \sigma} used by \code{"weighted"}. Defaults
#'         inherit from \code{"weighted"}.
#' }
#'
#' For multimodal or long-tailed \code{env_m}, consider z-scoring
#' \code{env_occ} and \code{env_m} so all variables have comparable
#' spread (e.g. z-score or quantile-rank); the ridge prior is then
#' uniform across axes.
#'
#' @section Optimization backend:
#'
#' \code{optimize_niche()} optimizes via \code{ucminfcpp::ucminf_xptr()} with
#' a compiled C++ objective built by \code{create_niche_obj_ptr()}. The full
#' \code{(theta -> mu, log sigma, v)} unpacking, \code{cvine_cholesky},
#' log-likelihood, and gradient are evaluated in pure C++ with no R-callback
#' overhead. For the weighted likelihood, gradient mode \code{"analytic"} uses
#' closed-form derivatives over \code{mu} and \code{log sigma} and central
#' finite differences over the C-vine partial-correlation block.
#'
#' @section KDE sampling (weighted model):
#'
#' The KDE that weights the M-background depends only on the environmental
#' data (not on \code{theta}), so the weights are computed once before
#' optimization. By default, the KDE reference (\code{m_kde_subsample}) and
#' the denominator subset (\code{m_subsample}) are both capped at
#' 10\,000 distinct background combinations -- the maximum the model can
#' usefully exploit even for very large rasters. A floor of
#' \code{max(500, 50 * 2^p)} rows is used as a heuristic minimum
#' representative sample (Silverman 1986; Wand & Jones 1995); below this
#' floor a \code{warning()} is emitted.
#'
#' @param env_occ Data frame of environmental values at presence points.
#' @param env_m Data frame of background environmental values. Required for
#'   \code{likelihood} in \code{"kde_bias_corrected"}, \code{"weighted"},
#'   \code{"skew_normal_weighted"}, or \code{"skew_t_weighted"};
#'   ignored for \code{"presence_only"}, \code{"skew_normal"}, and
#'   \code{"skew_t"}.
#' @param num_starts Integer. Number of Sobol starting points.
#' @param breadth Numeric in (0, 0.5). Controls the quantile range used to
#'   define starting bounds for \code{mu} parameters. Default \code{0.1}.
#' @param likelihood One of \code{"weighted"} (default; paper Eq. 5 + ridge),
#'   \code{"kde_bias_corrected"} (legacy KDE-bias-corrected formula),
#'   \code{"presence_only"}, \code{"skew_normal"} (presence-only
#'   multivariate skew-normal), \code{"skew_normal_weighted"} (paper
#'   Eq. 5 with skew-normal density + ridge prior on log sigma),
#'   \code{"skew_t"} (presence-only multivariate non-central skew-t
#'   via 32-node Gauss-Laguerre quadrature), or
#'   \code{"skew_t_weighted"} (paper Eq. 5 with NCST density + ridge
#'   prior).
#' @param grad Gradient strategy: \code{"auto"} (default) selects
#'   \code{"analytic"} for the Gaussian weighted models and \code{"central"}
#'   otherwise. Force one of \code{c("analytic", "central", "forward")} to
#'   override.
#' @param m_subsample,m_kde_subsample Optional integer or fraction in (0, 1].
#'   Resolved to \code{min(nrow(env_m), 10000)} when \code{NULL} (default).
#' @param seed Optional integer to make subsampling deterministic.
#' @param warm_start Logical. When \code{likelihood} is
#'   \code{"kde_bias_corrected"}, \code{"weighted"}, or
#'   \code{"skew_normal_weighted"}, run a quick presence-only fit first
#'   and prepend its \code{theta} (padded with \eqn{\alpha = 0} for the
#'   skew variants) as one extra starting point for the weighted
#'   multi-start. Cheap insurance against bad multistart luck on rough or
#'   multimodal weighted likelihoods (e.g. when \code{env_m} contains
#'   regions far from the occurrence cloud); does not replace the Sobol
#'   starts. Ignored for \code{likelihood = "presence_only"} and
#'   \code{"skew_normal"}. Default \code{TRUE}.
#' @param prior_log_sigma_lambda Numeric scalar (\code{>= 0}), strength of
#'   the ridge penalty on \eqn{\log \sigma} used by
#'   \code{likelihood = "weighted"}, \code{"skew_normal_weighted"}, and
#'   \code{"skew_t_weighted"}. The penalty is
#'   \eqn{\lambda \sum_k (\log \sigma_k - \log \hat\sigma_k)^2}. The centre
#'   \eqn{\log \hat\sigma_k} defaults to the presence-only fit's
#'   \eqn{\log \sigma_k} (shrinks the weighted fit toward the PO fit); see
#'   \code{prior_log_sigma_center}. Default \code{1.0} (weak). Larger
#'   values keep \eqn{\sigma} closer to the PO scale; \code{0} reduces the
#'   model to pure ML weighted normal (Eq. 5) and is not recommended
#'   (subject to Patil & Ord 1976 drift).
#' @param prior_log_sigma_center Optional numeric vector of length
#'   \code{ncol(env_occ)}. \code{NULL} (default) anchors the centre at the
#'   presence-only fit's \eqn{\log \sigma} when \code{warm_start = TRUE}
#'   (the PO fit is already run to seed the weighted multi-start); falls
#'   back to \code{log(apply(env_occ, 2, sd))} when \code{warm_start = FALSE}.
#' @param prior_mu_lambda Numeric scalar (\code{>= 0}), strength of a unit-free
#'   ridge penalty on \eqn{\mu} used by the weighted families
#'   (\code{"weighted"}, \code{"skew_normal_weighted"},
#'   \code{"skew_t_weighted"}). The penalty is
#'   \eqn{\lambda_\mu \sum_k \left((\mu_k - \hat\mu_k) / \hat\sigma_k\right)^2}
#'   with anchors \eqn{\hat\mu_k} and \eqn{\hat\sigma_k} from the PO fit
#'   (see \code{prior_mu_center}). Because the penalty is divided by the PO
#'   \eqn{\sigma_k}, it is scale-free across heterogeneous variables (e.g.
#'   bio1 in \eqn{^\circ C} vs bio12 in mm): \code{prior_mu_lambda = 1}
#'   allows \eqn{\mu} to drift ~1 PO standard deviation before the
#'   penalty pushes back. Default \code{1.0}. Set to \code{0} for
#'   unpenalized maximum likelihood.
#' @param prior_alpha_lambda Numeric scalar (\code{>= 0}), strength of a
#'   ridge penalty \eqn{\lambda_\alpha \sum_k \alpha_k^2} on the skew vector
#'   \eqn{\alpha}, shrinking toward the Gaussian sub-model. Applies to
#'   both presence-only and weighted skew families (\code{"skew_normal"},
#'   \code{"skew_normal_weighted"}, \code{"skew_t"},
#'   \code{"skew_t_weighted"}). Fixes the well-known unbounded-MLE
#'   pathology of the skew-normal direct parameterization (Azzalini 1985,
#'   Pewsey 2000). Default \code{0.1} (mild). Set to \code{0} for
#'   unpenalized maximum likelihood.
#' @param prior_mu_center Optional numeric vector of length
#'   \code{ncol(env_occ)}. \code{NULL} (default) anchors the centre at the
#'   presence-only fit's \eqn{\mu} when \code{warm_start = TRUE}; if
#'   \code{warm_start = FALSE} and \code{prior_mu_lambda > 0}, it must be
#'   supplied explicitly.
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
                           likelihood = c("weighted",
                                          "kde_bias_corrected",
                                          "presence_only",
                                          "skew_normal",
                                          "skew_normal_weighted",
                                          "skew_t",
                                          "skew_t_weighted"),
                           grad       = c("auto", "analytic",
                                          "central", "forward"),
                           m_subsample     = NULL,
                           m_kde_subsample = NULL,
                           seed            = NULL,
                           warm_start      = TRUE,
                           prior_log_sigma_lambda = 1.0,
                           prior_log_sigma_center = NULL,
                           prior_mu_lambda        = 1.0,
                           prior_mu_center        = NULL,
                           prior_alpha_lambda     = 0.1,
                           control = list(),
                           verbose = FALSE,
                           ...) {
  likelihood <- match.arg(likelihood)
  grad       <- match.arg(grad)

  # Convenience: treat the weighted family uniformly where logic is shared.
  is_weighted_family <- likelihood %in% c("kde_bias_corrected", "weighted",
                                          "skew_normal_weighted",
                                          "skew_t_weighted")
  # Convenience: skew likelihoods carry an extra alpha block (length p).
  is_skew_family <- likelihood %in% c("skew_normal", "skew_normal_weighted",
                                      "skew_t", "skew_t_weighted")
  # Convenience: skew-t likelihoods carry an additional log_r scalar.
  is_skew_t_family <- likelihood %in% c("skew_t", "skew_t_weighted")

  # Resolve `eta` from `...` so we can both forward it to the objective
  # functions (already done downstream) and persist it on the returned
  # object for `predict.nicher()`. Validate upfront — fails fast before
  # the optimizer wastes any compute on a malformed value.
  .dots <- list(...)
  unknown_dots <- setdiff(names(.dots), "eta")
  if (length(unknown_dots) > 0L) {
    stop("Unused argument(s): ", paste(unknown_dots, collapse = ", "))
  }
  eta <- if (!is.null(.dots$eta)) .dots$eta else 1.0
  if (!is.numeric(eta) || length(eta) != 1L ||
      !is.finite(eta) || eta <= 0) {
    stop("`eta` must be a single positive finite number.")
  }

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  if (!(likelihood %in% c("presence_only", "skew_normal", "skew_t"))) {
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

  # Validate ridge-prior strengths
  .validate_lambda <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop(sprintf("`%s` must be a single non-negative finite number.", name))
    }
  }
  .validate_lambda(prior_log_sigma_lambda, "prior_log_sigma_lambda")
  .validate_lambda(prior_mu_lambda,        "prior_mu_lambda")
  .validate_lambda(prior_alpha_lambda,     "prior_alpha_lambda")

  p_occ <- ncol(env_occ)
  # Validate user-supplied centres upfront. Defaults (NULL) are resolved
  # later -- after the warm-start PO fit, which is the source of the
  # anchor for the weighted families. For non-weighted families we
  # materialise a zero vector so downstream code always sees numeric(p).
  .validate_center <- function(x, name) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.numeric(x) || length(x) != p_occ || any(!is.finite(x))) {
      stop(sprintf(
        "`%s` must be a finite numeric vector of length ncol(env_occ).",
        name))
    }
  }
  .validate_center(prior_log_sigma_center, "prior_log_sigma_center")
  .validate_center(prior_mu_center,        "prior_mu_center")

  # Resolve grad="auto". Analytic gradients exist for the Gaussian weighted
  # family (kde_bias_corrected, weighted) only; the skew-normal and skew-t
  # kernels ship with finite-difference gradients (closed-form Azzalini
  # gradients deferred to a follow-up; for skew-t the integral over the
  # chi^2_r mixing variable would also need a quadrature-aware derivative).
  resolved_grad <- if (grad == "auto") {
    if (likelihood %in% c("kde_bias_corrected", "weighted"))
      "analytic"
    else
      "central"
  } else grad

  # ------------------------------------------------------------------
  # Sobol starts
  # ------------------------------------------------------------------
  quant_vec <- c(breadth, 0.5, 1.0 - breadth)
  starts_df <- start_theta_multiple(
    env_data   = env_occ,
    num_starts = num_starts,
    quant_vec  = quant_vec,
    method     = "sobol",
    skew       = is_skew_family,
    skew_t     = is_skew_t_family
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
  # Resolve KDE subsamples / weights ONCE per fit (weighted family only)
  # ------------------------------------------------------------------
  weighted_inputs <- NULL
  if (is_weighted_family) {
    weighted_inputs <- .resolve_weighted_inputs(
      env_occ         = env_occ,
      env_m           = env_m,
      m_subsample     = m_subsample,
      m_kde_subsample = m_kde_subsample,
      seed            = seed
    )
  }

  # ------------------------------------------------------------------
  # Warm-start (weighted only): prepend the presence-only optimum
  # ------------------------------------------------------------------
  # The presence-only likelihood is unimodal and converges cleanly even on
  # multimodal env_m (e.g. example_vicugna where bio12 in M is bimodal but
  # in env_occ is unimodal). Its theta lives in the same parameter space as
  # the weighted model (same mu, log_sigma, v layout) so we can hand it
  # straight to the weighted optimizer as one additional starting point.
  # Cheap insurance against bad multistart luck on rough weighted surfaces;
  # does NOT replace the Sobol starts.
  po_anchor_mu        <- NULL
  po_anchor_log_sigma <- NULL
  if (warm_start && is_weighted_family) {
    if (verbose) message("Warm-start: running presence-only fit ...")
    po_fit <- tryCatch(
      optimize_niche(
        env_occ    = env_occ,
        env_m      = NULL,
        num_starts = min(20L, as.integer(num_starts)),
        breadth    = breadth,
        likelihood = "presence_only",
        grad       = if (grad == "auto") "central" else grad,
        seed       = seed,
        warm_start = FALSE,
        prior_log_sigma_lambda = 0.0,
        prior_mu_lambda        = 0.0,
        prior_alpha_lambda     = 0.0,
        control    = control,
        verbose    = FALSE,
        ...
      ),
      error = function(e) {
        warning("Warm-start presence-only fit failed: ", conditionMessage(e),
                "; continuing with Sobol starts only.")
        NULL
      }
    )
    if (!is.null(po_fit) && po_fit$best$convergence %in% c(1L, 2L)) {
      warm_theta <- po_fit$best$theta
      # Extract the PO fit's (mu, log_sigma) blocks as anchors for the
      # penalised-MLE centres. Gaussian PO theta layout is
      # [mu(p), log_sigma(p), v(p(p-1)/2)].
      po_anchor_mu        <- as.numeric(warm_theta[seq_len(p_occ)])
      po_anchor_log_sigma <- as.numeric(warm_theta[p_occ + seq_len(p_occ)])
      # For skew_normal_weighted / skew_t_weighted, the PO fit has no alpha
      # (or log_r) block; pad so the warm-start lands at the symmetric
      # Gaussian point. SN_k(μ, Σ, α=0) ≡ N_k(μ, Σ); for skew-t we
      # additionally start at log_r = log(10) (Sobol centre).
      target_len <- length(starts_list[[1L]])
      n_pad <- target_len - length(warm_theta)
      if (is_skew_family && n_pad > 0L) {
        pad_alpha <- rep(0.0, p_occ)
        pad_log_r <- if (is_skew_t_family) log(10) else numeric(0)
        pad <- c(pad_alpha, pad_log_r)[seq_len(n_pad)]
        warm_theta <- c(warm_theta, pad)
      }
      if (length(warm_theta) == length(starts_list[[1L]])) {
        names(warm_theta) <- names(starts_list[[1L]])
        starts_list <- c(list(warm_theta), starts_list)
        if (verbose) {
          message(sprintf(
            "Warm-start theta added (PO loglik = %.6f).",
            po_fit$best$loglik
          ))
        }
      }
    }
  }

  # Resolve default centres for the ridge priors. Priority:
  #   1. User-supplied centre (validated above).
  #   2. PO anchor (from warm-start, if it converged).
  #   3. Fallback: log(sd(env_occ)) for log_sigma; mu has no fallback --
  #      if prior_mu_lambda > 0 the user MUST supply either warm_start
  #      or an explicit prior_mu_center (else error).
  if (is.null(prior_log_sigma_center)) {
    if (!is.null(po_anchor_log_sigma)) {
      prior_log_sigma_center <- po_anchor_log_sigma
      names(prior_log_sigma_center) <- colnames(env_occ)
    } else if (is_weighted_family && prior_log_sigma_lambda > 0) {
      sds <- apply(as.matrix(env_occ), 2L, stats::sd, na.rm = TRUE)
      if (any(!is.finite(sds)) || any(sds <= 0)) {
        stop("Cannot derive default `prior_log_sigma_center`: some env_occ ",
             "columns have zero or non-finite sd. Pass it explicitly or ",
             "enable warm_start.")
      }
      prior_log_sigma_center <- log(sds)
      names(prior_log_sigma_center) <- colnames(env_occ)
    } else {
      prior_log_sigma_center <- rep(0.0, p_occ)
    }
  }
  if (is.null(prior_mu_center)) {
    if (!is.null(po_anchor_mu)) {
      prior_mu_center <- po_anchor_mu
      names(prior_mu_center) <- colnames(env_occ)
    } else if (is_weighted_family && prior_mu_lambda > 0) {
      stop("prior_mu_center is NULL and no PO warm-start is available. ",
           "Either enable warm_start = TRUE, pass prior_mu_center ",
           "explicitly, or set prior_mu_lambda = 0.")
    } else {
      prior_mu_center <- rep(0.0, p_occ)
    }
  }

  # ------------------------------------------------------------------
  # Run all starts
  # ------------------------------------------------------------------
  helper <- .optimize_niche_helper_cpp

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
      prior_log_sigma_center = prior_log_sigma_center,
      prior_log_sigma_lambda = prior_log_sigma_lambda,
      prior_mu_center        = prior_mu_center,
      prior_mu_lambda        = prior_mu_lambda,
      prior_alpha_lambda     = prior_alpha_lambda,
      eta = eta
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
  # Internal safeguard: validate compiled optimizer result against ucminf
  # ------------------------------------------------------------------
  if (best$convergence %in% c(1L, 2L)) {
    .validate_xptr_result(
      best            = best,
      env_occ         = env_occ,
      env_m           = env_m,
      likelihood      = likelihood,
      weighted_inputs = weighted_inputs,
      ctrl            = ctrl,
      prior_log_sigma_center = prior_log_sigma_center,
      prior_log_sigma_lambda = prior_log_sigma_lambda,
      prior_mu_center        = prior_mu_center,
      prior_mu_lambda        = prior_mu_lambda,
      prior_alpha_lambda     = prior_alpha_lambda,
      eta = eta
    )
  }

  # ------------------------------------------------------------------
  # Information-criteria payload: un-penalised log-likelihood evaluated
  # at the converged theta, plus nobs and cheap fingerprints of
  # env_occ / env_m so `compare_nicher()` can refuse non-comparable
  # inputs. `loglik_unpenalised` strips the ridge so AIC/BIC see the
  # bare data likelihood; for `presence_only` and `kde_bias_corrected`
  # there is no penalty so it equals `best$loglik`. See `?logLik.nicher`.
  # ------------------------------------------------------------------
  best$loglik_unpenalised <- tryCatch(
    .compute_unpenalised_loglik(
      theta           = best$theta,
      likelihood      = likelihood,
      env_occ         = env_occ,
      env_m           = env_m,
      weighted_inputs = weighted_inputs,
      eta             = eta
    ),
    error = function(e) {
      warning("Could not compute un-penalised log-likelihood for ",
              "logLik()/AIC()/BIC(): ", conditionMessage(e),
              call. = FALSE)
      NA_real_
    }
  )
  best$nobs                <- nrow(env_occ)
  best$env_occ_fingerprint <- .env_fingerprint(env_occ)
  best$env_m_fingerprint   <- if (likelihood %in%
                                  c("presence_only", "skew_normal",
                                    "skew_t")) {
    NULL
  } else {
    .env_fingerprint(env_m)
  }

  # ------------------------------------------------------------------
  # Persist fit recipe so `cv_nicher()` can replay the same likelihood
  # / knob configuration on a subset of `env_occ`. We deliberately
  # exclude the bulky data (env_occ, env_m); those are passed back to
  # `cv_nicher()` and validated against the fingerprints above.
  # ------------------------------------------------------------------
  fit_args <- list(
    num_starts             = num_starts,
    breadth                = breadth,
    likelihood             = likelihood,
    grad                   = resolved_grad,
    m_subsample            = m_subsample,
    m_kde_subsample        = m_kde_subsample,
    seed                   = seed,
    warm_start             = warm_start,
    prior_log_sigma_lambda = prior_log_sigma_lambda,
    prior_log_sigma_center = prior_log_sigma_center,
    prior_mu_lambda        = prior_mu_lambda,
    prior_mu_center        = prior_mu_center,
    prior_alpha_lambda     = prior_alpha_lambda,
    control                = control,
    eta                    = eta
  )

  if (verbose) {
    message(sprintf(
      "Best log-likelihood: %.6f (convergence = %d)",
      best$loglik, best$convergence
    ))
  }

  out <- new_nicher(
    solutions  = solutions,
    best       = best,
    likelihood = likelihood,
    n_starts   = num_starts,
    var_names  = colnames(env_occ),
    eta        = eta
  )
  out$fit_args <- fit_args
  out
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
                                  weighted_inputs = NULL,
                                  prior_log_sigma_center = NULL,
                                  prior_log_sigma_lambda = 0.0,
                                  prior_mu_center        = NULL,
                                  prior_mu_lambda        = 0.0,
                                  prior_alpha_lambda     = 0.0,
                                  ...) {
  # Mirror the per-start helper: for non-weighted families the mu/log_sigma
  # penalties are disabled so the validator evaluates the *same* objective
  # the xptr kernel minimised. (The alpha ridge still applies.)
  if (!(likelihood %in% c("weighted", "skew_normal_weighted",
                          "skew_t_weighted"))) {
    prior_mu_lambda        <- 0.0
    prior_log_sigma_lambda <- 0.0
  }
  pmc <- if (!is.null(prior_mu_center)) as.numeric(prior_mu_center) else NULL
  plsc <- if (!is.null(prior_log_sigma_center))
            as.numeric(prior_log_sigma_center) else NULL
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
    kde_bias_corrected = function(theta) {
      loglik_niche_math_kde_bias_corrected(
        theta, env_occ = env_occ, env_m = env_m, neg = TRUE,
        den_idx       = weighted_inputs$den_idx,
        kde_idx       = weighted_inputs$kde_idx,
        precomp_w_den = weighted_inputs$w_den,
        ...
      )
    },
    weighted = {
      env_m_mat <- as.matrix(env_m)
      M_den <- env_m_mat[weighted_inputs$den_idx, , drop = FALSE]
      function(theta) {
        loglik_niche_math_weighted_cpp(
          theta = theta,
          env_occ = as.matrix(env_occ),
          M_den   = M_den,
          w_occ   = weighted_inputs$w_occ,
          w_den   = weighted_inputs$w_den,
          prior_log_sigma_center = as.numeric(prior_log_sigma_center),
          prior_log_sigma_lambda = prior_log_sigma_lambda,
          eta = if (!is.null(list(...)$eta)) list(...)$eta else 1.0,
          prior_mu_center    = pmc,
          prior_mu_lambda    = prior_mu_lambda,
          prior_alpha_lambda = prior_alpha_lambda
        )
      }
    },
    skew_normal = function(theta) {
      loglik_niche_math_skew_normal_cpp(
        theta = theta,
        env_occ = as.matrix(env_occ),
        eta = if (!is.null(list(...)$eta)) list(...)$eta else 1.0,
        prior_mu_center        = pmc,
        prior_mu_lambda        = prior_mu_lambda,
        prior_log_sigma_center = plsc,
        prior_log_sigma_lambda = prior_log_sigma_lambda,
        prior_alpha_lambda     = prior_alpha_lambda
      )
    },
    skew_normal_weighted = {
      env_m_mat <- as.matrix(env_m)
      M_den <- env_m_mat[weighted_inputs$den_idx, , drop = FALSE]
      function(theta) {
        loglik_niche_math_skew_normal_weighted_cpp(
          theta = theta,
          env_occ = as.matrix(env_occ),
          M_den   = M_den,
          w_occ   = weighted_inputs$w_occ,
          w_den   = weighted_inputs$w_den,
          prior_log_sigma_center = as.numeric(prior_log_sigma_center),
          prior_log_sigma_lambda = prior_log_sigma_lambda,
          eta = if (!is.null(list(...)$eta)) list(...)$eta else 1.0,
          prior_mu_center    = pmc,
          prior_mu_lambda    = prior_mu_lambda,
          prior_alpha_lambda = prior_alpha_lambda
        )
      }
    },
    skew_t = function(theta) {
      loglik_niche_math_skew_t_cpp(
        theta = theta,
        env_occ = as.matrix(env_occ),
        eta = if (!is.null(list(...)$eta)) list(...)$eta else 1.0,
        prior_mu_center        = pmc,
        prior_mu_lambda        = prior_mu_lambda,
        prior_log_sigma_center = plsc,
        prior_log_sigma_lambda = prior_log_sigma_lambda,
        prior_alpha_lambda     = prior_alpha_lambda
      )
    },
    skew_t_weighted = {
      env_m_mat <- as.matrix(env_m)
      M_den <- env_m_mat[weighted_inputs$den_idx, , drop = FALSE]
      function(theta) {
        loglik_niche_math_skew_t_weighted_cpp(
          theta = theta,
          env_occ = as.matrix(env_occ),
          M_den   = M_den,
          w_occ   = weighted_inputs$w_occ,
          w_den   = weighted_inputs$w_den,
          prior_log_sigma_center = as.numeric(prior_log_sigma_center),
          prior_log_sigma_lambda = prior_log_sigma_lambda,
          eta = if (!is.null(list(...)$eta)) list(...)$eta else 1.0,
          prior_mu_center    = pmc,
          prior_mu_lambda    = prior_mu_lambda,
          prior_alpha_lambda = prior_alpha_lambda
        )
      }
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
                                       prior_log_sigma_center = NULL,
                                       prior_log_sigma_lambda = 0.0,
                                       prior_mu_center        = NULL,
                                       prior_mu_lambda        = 0.0,
                                       prior_alpha_lambda     = 0.0,
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
  if (likelihood %in% c("kde_bias_corrected", "weighted",
                        "skew_normal_weighted", "skew_t_weighted") &&
      !is.null(weighted_inputs)) {
    den_idx       <- weighted_inputs$den_idx
    kde_idx       <- weighted_inputs$kde_idx
    precomp_w_occ <- weighted_inputs$w_occ
    precomp_w_den <- weighted_inputs$w_den
  }

  # Only ridge-penalised variants (weighted, skew_normal_weighted,
  # skew_t_weighted) use the mu/log_sigma penalty; alpha penalty applies
  # to any skew family (PO or weighted). For simplicity we let the C++
  # side decide what to use: we pass everything and rely on lambdas = 0
  # or missing alpha blocks to no-op.
  plsl <- as.numeric(prior_log_sigma_lambda)
  pml  <- as.numeric(prior_mu_lambda)
  pal  <- as.numeric(prior_alpha_lambda)
  # Disable mu/log_sigma penalties for non-weighted families (the kernels
  # that don't carry them would ignore the terms anyway, but passing
  # nonzero lambdas with NULL centres would trigger validation errors).
  if (!(likelihood %in% c("weighted", "skew_normal_weighted",
                          "skew_t_weighted"))) {
    plsl <- 0.0
    pml  <- 0.0
  }
  plsc <- if (plsl > 0 && !is.null(prior_log_sigma_center))
            as.numeric(prior_log_sigma_center) else NULL
  pmc  <- if (pml  > 0 && !is.null(prior_mu_center))
            as.numeric(prior_mu_center) else NULL

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
    gradstep      = gs,
    prior_log_sigma_center = plsc,
    prior_log_sigma_lambda = plsl,
    prior_mu_center        = pmc,
    prior_mu_lambda        = pml,
    prior_alpha_lambda     = pal
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
