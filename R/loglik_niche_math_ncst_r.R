#' Reference (pure-R) negative log-likelihood for the NCST density
#'
#' Evaluates the Non-Central Skew t (NCST) negative log-likelihood using
#' pure base-R arithmetic (no C++ calls except for the C-vine Cholesky).
#' This is the Stage 1 reference implementation following the SSDLC:
#' correctness is everything, speed is irrelevant.
#'
#' The NCST stochastic form (Hasan & Chen 2025, Definition 1) is:
#'   T = X / sqrt(Y / r),
#'   X ~ SN_k(xi, Omega, alpha),
#'   Y ~ chi^2_r,
#' where the location xi enters INSIDE the skew-normal BEFORE chi-squared
#' scaling. This differs from Azzalini/Branco & Dey's skew-t where
#' location is added AFTER scaling.
#'
#' @param theta Numeric vector of length 3*p + p*(p-1)/2 + 1:
#'   \code{[xi(p), log_sigma(p), v(n_v), alpha(p), log_r(1)]}.
#' @param env_occ Numeric matrix (n x p) of occurrence data.
#' @param eta Scalar shape parameter for the C-vine prior (default 1).
#'
#' @return Scalar negative log-likelihood (same constant-dropping
#'   convention as the C++ kernel).
#'
#' @references
#' Hasan, M. R. & Chen, M.-H. (2025). Flexible Modeling of Multivariate
#' Skewed and Heavy-Tailed Data via a Non-Central Skew t Distribution.
#' arXiv:2507.10465v1.
#'
#' @keywords internal
loglik_niche_math_ncst_r <- function(theta, env_occ, eta = 1.0) {
  env_occ <- as.matrix(env_occ)
  p <- ncol(env_occ)
  n_v <- p * (p - 1) / 2
  n_occ <- nrow(env_occ)

  # Parse theta
  xi        <- theta[seq_len(p)]
  log_sigma <- theta[p + seq_len(p)]
  sigma     <- exp(log_sigma)
  v         <- theta[2 * p + seq_len(n_v)]
  alpha     <- theta[2 * p + n_v + seq_len(p)]
  log_r     <- theta[3 * p + n_v + 1]
  r         <- exp(log_r)

  # Build covariance Cholesky: L_cov = diag(sigma) %*% L_corr
  L_corr <- cvine_cholesky(v, d = p, eta = eta)
  L_cov  <- diag(sigma, nrow = p) %*% L_corr

  log_det <- 2 * sum(log(diag(L_cov)))

  # Precompute: L_cov^{-1} xi and per-observation L_cov^{-1} t_i
  y_xi <- forwardsolve(L_cov, xi)
  qxi  <- sum(y_xi^2)
  a_over_s <- alpha / sigma
  zxi <- sum(a_over_s * xi)

  # GL nodes and weights (same 32-node table as C++)
  gl_nodes <- c(
    4.448936583326738858e-02, 2.345261095196182477e-01,
    5.768846293018863314e-01, 1.072448753817818012e+00,
    1.722408776444645628e+00, 2.528336706425794222e+00,
    3.492213273021993913e+00, 4.616456769749767375e+00,
    5.903958504174243949e+00, 7.358126733186240997e+00,
    8.982940924212597267e+00, 1.078301863253997084e+01,
    1.276369798674272538e+01, 1.493113975552255823e+01,
    1.729245433671531629e+01, 1.985586094033605420e+01,
    2.263088901319677504e+01, 2.562863602245924710e+01,
    2.886210181632347371e+01, 3.234662915396473437e+01,
    3.610049480575197123e+01, 4.014571977153944005e+01,
    4.450920799575494158e+01, 4.922439498730864216e+01,
    5.433372133339690890e+01, 5.989250916213401865e+01,
    6.597537728793504641e+01, 7.268762809066271302e+01,
    8.018744697791352394e+01, 8.873534041789240234e+01,
    9.882954286828397983e+01, 1.117513980979376953e+02
  )
  gl_weights <- c(
    1.092183419524251492e-01, 2.104431079387923953e-01,
    2.352132296698429825e-01, 1.959033359728735435e-01,
    1.299837862860684090e-01, 7.057862386571556179e-02,
    3.176091250917396913e-02, 1.191821483483820367e-02,
    3.738816294611417134e-03, 9.808033066149245659e-04,
    2.148649188013577994e-04, 3.920341967987801200e-05,
    5.934541612868447911e-06, 7.416404578667340129e-07,
    7.604567879120552326e-08, 6.350602226625617791e-09,
    4.281382971040738227e-10, 2.305899491891237043e-11,
    9.799379288726772448e-13, 3.237801657729156650e-14,
    8.171823443420548711e-16, 1.542133833393811132e-17,
    2.119792290163546927e-19, 2.054429673788003569e-21,
    1.346982586637361709e-23, 5.661294130397462134e-26,
    1.418560545462968432e-28, 1.913375494453994255e-31,
    1.192248760098201249e-34, 2.671511219239662518e-38,
    1.338616942106424308e-42, 4.510536193898778441e-48
  )

  n_gl <- length(gl_nodes)
  a_exp <- 0.5 * (p + r) - 1  # exponent for u in GL integrand

  sum_log_I <- 0
  for (i in seq_len(n_occ)) {
    t_i <- env_occ[i, ]
    y_t_i <- forwardsolve(L_cov, t_i)
    qt_i  <- sum(y_t_i^2)
    cross_i <- sum(y_t_i * y_xi)
    zt_i <- sum(a_over_s * t_i)

    log_terms <- numeric(n_gl)
    for (q in seq_len(n_gl)) {
      u_q <- gl_nodes[q]
      w_q <- gl_weights[q]
      if (w_q <= 0) {
        log_terms[q] <- -Inf
        next
      }
      s_q <- sqrt(2 * u_q / r)

      # q_iq = ||s_q * y_t_i - y_xi||^2
      q_iq <- s_q^2 * qt_i - 2 * s_q * cross_i + qxi
      z_iq <- s_q * zt_i - zxi

      log_terms[q] <- log(w_q) + a_exp * log(u_q) -
        0.5 * q_iq + stats::pnorm(z_iq, log.p = TRUE)
    }

    # log-sum-exp
    max_t <- max(log_terms)
    if (!is.finite(max_t)) {
      sum_log_I <- sum_log_I + (-1e300)
    } else {
      sum_log_I <- sum_log_I + max_t + log(sum(exp(log_terms - max_t)))
    }
  }

  # Negative log-likelihood (dropping theta-independent constants:
  # -(k/2+1) n log(2) and (k/2) n log(2*pi), matching existing skew-t)
  n <- as.double(n_occ)
  neg_loglik <- 0.5 * n * log_det +
    0.5 * n * p * log_r +
    n * lgamma(0.5 * r) -
    sum_log_I

  if (!is.finite(neg_loglik)) neg_loglik <- 1e300
  neg_loglik
}
