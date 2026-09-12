#' @title Discrete Beta Model
#'
#' @description
#' The Discrete Beta (DBT) distribution models ordinal rating data on a fixed
#' integer scale \eqn{R \in \{1, \dots, k\}} by discretizing an underlying
#' continuous Beta distribution at \eqn{k - 1} evenly-spaced thresholds
#' \eqn{\gamma_j = j / k}. Unlike proportional-odds style models, which fix
#' the underlying distribution and estimate the thresholds, the Discrete Beta
#' fixes the thresholds and estimates the two shape parameters of the
#' underlying Beta distribution instead. This keeps the model parsimonious
#' (only 2 parameters) while remaining flexible enough to reproduce "U" and "J"
#' (non-monotonic convex) shapes that are common in rating data and that
#' proportional-odds models cannot capture (Sciandra et al., 2024).
#'
#' @param n Number of simulated values.
#' @param x,q Vector of quantiles (integer ratings between 1 and `k`, or 0 if
#'   `pzero > 0`).
#' @param p Vector of probabilities.
#' @inheritParams rcogmod_betagate
#' @param k Number of rating categories (a positive integer, `k >= 1`), i.e.
#'   the response scale runs from 1 to `k`.
#' @param pzero Probability of an additional "hurdle" point mass at 0, on top of
#'   the `1:k` rating scale. Defaults to `0`, in which case the distribution
#'   reduces to the pure Discrete Beta model. Useful for rating scales that
#'   include an extra "zero" category (e.g., "not applicable" or a genuine
#'   zero response) that is not part of the underlying `1:k` continuum.
#' @param log,log.p Logical; if `TRUE`, probabilities/densities are returned
#'   on the log scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(R \le q)}, otherwise \eqn{P(R > q)}.
#'
#' @details
#' Writing \eqn{\alpha = \mu \phi} and \eqn{\beta = (1 - \mu)\phi} for the
#' shape parameters of the underlying Beta distribution, the probability mass
#' function is (Sciandra et al., 2024, eq. 2)
#' \deqn{P(R = j) = F_B(j/k; \alpha, \beta) - F_B((j-1)/k; \alpha, \beta), \quad j = 1, \dots, k}
#' where \eqn{F_B} is the Beta CDF.
#'
#' `rcogmod_betadiscrete()` uses the equivalent, faster generative representation:
#' draw a continuous \eqn{X \sim Beta(\alpha, \beta)} and set
#' \eqn{R = \lceil k X \rceil}, clipped to `[1, k]`.
#'
#' When `pzero > 0`, a hurdle is added at 0: with probability `pzero` the
#' response is 0, and with probability `1 - pzero` it is generated from the
#' Discrete Beta distribution described above, i.e.
#' \deqn{P(R = 0) = \code{pzero}, \quad P(R = j) = (1 - \code{pzero}) \times [F_B(j/k) - F_B((j-1)/k)], \quad j = 1, \dots, k}
#'
#' **Special cases:**
#' - `mu = 0.5`, `phi = 1` (i.e. `alpha = beta = 1`): reduces to the discrete
#'   Uniform distribution on `1:k`.
#' - `alpha, beta < 1`: "U"/"J"-shaped, with mass concentrated in the tails.
#' - `alpha, beta > 1`: concave, with mass concentrated around the middle
#'   category.
#' - `phi -> Inf` (with `mu` fixed): mass concentrates on a single category.
#' - `pzero = 0`: reduces to the pure Discrete Beta model (no hurdle).
#'
#' @return `dcogmod_betadiscrete()` returns the probability mass;
#'   `pcogmod_betadiscrete()` returns the cumulative probability;
#'   `qcogmod_betadiscrete()` returns the quantile (an integer between 0 and
#'   `k`); `rcogmod_betadiscrete()` returns simulated ratings. All are numeric
#'   vectors, vectorized over `x`/`q`/`p`, `mu`, `phi`, `pzero` and `k`.
#'   `cogmod_betadiscrete()` returns a `brms::custom_family` object, to put on
#'   a `brms::bf()` formula. `cogmod_betadiscrete_stanvars()` returns a
#'   `brms::stanvars` object holding the family's Stan `functions` block, to
#'   pass to `brms::brm()`, and `cogmod_betadiscrete_lpmf_expose()` compiles
#'   that Stan code and returns it as an R function, for checking the mass
#'   function outside of a model. The remaining functions are `brms`
#'   post-processing methods, called by `brms` rather than directly:
#'   `log_lik_cogmod_betadiscrete()` returns a numeric vector holding one
#'   log-likelihood value per posterior draw for observation `i`,
#'   `posterior_predict_cogmod_betadiscrete()` a draws x 1 matrix of ratings
#'   simulated for observation `i`, and
#'   `posterior_epred_cogmod_betadiscrete()` a draws x observations matrix of
#'   expected ratings.
#'
#' @references
#' - Sciandra, M., Fasola, S., Albano, A., Di Maria, C., & Plaia, A. (2024).
#'   Discrete Beta and Shifted Beta-Binomial models for rating and ranking
#'   data. Environmental and Ecological Statistics, 31, 317-338.
#'   \doi{10.1007/s10651-023-00592-5}
#'
#' @examples
#' x <- 1:10
#' probs <- dcogmod_betadiscrete(x, mu = 0.66, phi = 3.51, k = 10)
#' barplot(probs, names.arg = x)
#'
#' y <- rcogmod_betadiscrete(1000, mu = 0.66, phi = 3.51, k = 10)
#' hist(y, breaks = 0:10)
#'
#' # discrete Uniform special case
#' dcogmod_betadiscrete(1:5, mu = 0.5, phi = 1, k = 5)
#'
#' # hurdle at zero: 20% chance of a 0, otherwise pure Discrete Beta
#' dcogmod_betadiscrete(0:5, mu = 0.66, phi = 3.51, k = 5, pzero = 0.2)
#'
#' @export
rcogmod_betadiscrete <- function(n, mu = 0.5, phi = 3, k = 5, pzero = 0) {
  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) {
    stop("n must be a positive integer.")
  }
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n_out <- max(n, length(mu), length(phi), length(k), length(pzero))
  if (n_out > 1 || n > 1) {
    mu <- rep(mu, length.out = n_out)
    phi <- rep(phi, length.out = n_out)
    k <- rep(k, length.out = n_out)
    pzero <- rep(pzero, length.out = n_out)
  } else {
    n_out <- n
  }

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Simulation ---
  # Discretize a continuous Beta draw at the fixed thresholds gamma_j = j / k
  x <- stats::rbeta(n_out, shape1 = alpha, shape2 = beta)
  out <- ceiling(x * k)
  out <- pmax(1L, pmin(out, k)) # guard against floating-point edge cases at 0 and 1

  # Hurdle: overwrite with 0 wherever the hurdle is crossed
  is_zero <- stats::runif(n_out) < pzero
  out[is_zero] <- 0L

  if (n_out == n && n > 0) {
    return(out[seq_len(n)])
  }
  out
}


#' @rdname rcogmod_betadiscrete
#' @export
dcogmod_betadiscrete <- function(x, mu = 0.5, phi = 3, k = 5, pzero = 0, log = FALSE) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(x)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Probability Mass ---
  # P(R = j) = F_B(j/k) - F_B((j-1)/k); handled at the boundaries explicitly
  # so that floating-point pbeta(0, ...) / pbeta(1, ...) calls are avoided.
  upper <- ifelse(x >= k, 1, stats::pbeta(x / k, alpha, beta))
  lower <- ifelse(x <= 1, 0, stats::pbeta((x - 1) / k, alpha, beta))

  prob <- (1 - pzero) * (upper - lower)
  # Hurdle point mass at 0
  is_zero <- x == 0
  prob[is_zero] <- pzero[is_zero]

  invalid <- (x < 0 | x > k | x != floor(x))
  prob[invalid] <- 0
  prob <- pmax(0, prob) # guard against tiny negative values from cancellation

  if (log) {
    prob <- ifelse(prob > 0, log(prob), -Inf)
  }
  prob
}


#' @rdname rcogmod_betadiscrete
#' @export
pcogmod_betadiscrete <- function(
  q,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(q)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Cumulative Probability ---
  qq <- pmin(pmax(floor(q), 0), k)
  p_orig <- ifelse(
    qq <= 0,
    0,
    ifelse(qq >= k, 1, stats::pbeta(qq / k, alpha, beta))
  )
  # P(R <= q) = pzero + (1 - pzero) * F_B(q/k), with P(R <= 0) = pzero
  p <- pzero + (1 - pzero) * p_orig
  p[q < 0] <- 0

  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}


#' @rdname rcogmod_betadiscrete
#' @export
qcogmod_betadiscrete <- function(
  p,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (any(p < 0 | p > 1)) {
    stop("p must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(p)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Quantile ---
  # Because the thresholds gamma_j = j/k are fixed, the discretized quantile
  # is simply the continuous Beta quantile mapped through the same thresholds.
  # Below the hurdle mass, the quantile is 0; above it, rescale p onto the
  # [0, 1] scale of the underlying (non-hurdle) Discrete Beta distribution.
  p_resc <- pmax(0, (p - pzero) / (1 - pzero))
  xp <- stats::qbeta(p_resc, alpha, beta)
  out <- pmax(1L, pmin(k, ceiling(k * xp)))
  out[p <= pzero] <- 0L
  out
}


# Stan Functions and Custom Family for brms (Discrete Beta)

# Stanvars ----------------------------------------------------------------

#' @keywords internal
.cogmod_betadiscrete_lpmf <- function() {
  "
// Log probability mass function for the (hurdle) Discrete Beta distribution
// (Sciandra et al., 2024, Sect. 3.1)
//   y     : observed rating, integer in {0, 1, ..., k}. y = 0 is only valid
//           when pzero > 0 (hurdle point mass below the 1..k rating scale)
//   mu    : mean of the underlying Beta distribution (0 < mu < 1); the
//           'liking' indicator on the logit scale
//   phi   : precision of the underlying Beta distribution (alpha + beta > 0);
//           the 'agreement' indicator on the log scale
//   pzero : probability of the hurdle point mass at 0 (0 <= pzero < 1)
//   k     : number of rating categories (fixed, passed in as data)
real cogmod_betadiscrete_lpmf(int y, real mu, real phi, real pzero, int k) {
  real alpha;
  real beta_par;
  real upper_lcdf;
  real lower_lcdf;

  if (y < 0 || y > k) {
    reject(\"cogmod_betadiscrete_lpmf: y must be an integer between 0 and k; found y = \", y);
  }

  if (y == 0) {
    return log(pzero);
  }

  alpha = mu * phi * 2;
  beta_par = (1 - mu) * phi * 2;

  // P(R = y) = F_B(y/k) - F_B((y-1)/k), computed on the log scale for
  // numerical stability via log_diff_exp(log(upper), log(lower)).
  upper_lcdf = (y == k) ? 0.0 : beta_lcdf(y * 1.0 / k | alpha, beta_par);
  lower_lcdf = (y == 1) ? negative_infinity() : beta_lcdf((y - 1) * 1.0 / k | alpha, beta_par);

  return log1m(pzero) + log_diff_exp(upper_lcdf, lower_lcdf);
}
"
}


#' @rdname rcogmod_betadiscrete
#' @examples
#' \donttest{
#' # Exposing the Stan function needs cmdstanr and a CmdStan toolchain,
#' # which live outside CRAN - see the package website to install them.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   lpmf <- cogmod_betadiscrete_lpmf_expose()
#'   lpmf(y = 7, mu = 0.66, phi = 3.51, pzero = 0, k = 10)
#' }
#' }
#'
#' @export
cogmod_betadiscrete_lpmf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .cogmod_betadiscrete_lpmf(),
    "\n}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$cogmod_betadiscrete_lpmf
}

#' @rdname rcogmod_betadiscrete
#' @export
cogmod_betadiscrete_stanvars <- function() {
  brms::stanvar(scode = .cogmod_betadiscrete_lpmf(), block = "functions")
}


#' @rdname rcogmod_betadiscrete
#' @param link_mu,link_phi,link_pzero Link functions for the parameters. `pzero`
#'   defaults to a `"logit"` link. By default (i.e., if `pzero` is not included
#'   in the `brms::bf()` formula), it is estimated as a single, intercept-only
#'   value shared across all observations (as is done for `pmid` in
#'   [cogmod_choco()]); it can instead be given predictors to let it vary (`pzero ~ x`),
#'   or fixed to a constant -- e.g., `pzero = 0`, recovering the pure Discrete
#'   Beta model -- directly in `brms::bf()` (as is done for `pmid` in [cogmod_choco()]).
#'
#' @details
#' Note that `y = 0` is always handled by `pzero` alone, and `k` always refers to the
#' number of categories of the *non-zero* `1:k` part of the scale. What
#' *does* require some care is deciding what `k` should be and whether to
#' estimate or fix `pzero`, depending on how the zero in your data arose:
#' - Scale is `0:N` and 0 is *not* a hurdle (just the lowest of `N + 1`
#'   ordinary ordinal categories, e.g., a 0-10 rating scale with no excess of
#'   zeros): recode the data to `1:(N + 1)` (add 1 to every response), use
#'   `vint(N + 1)`, and fix `pzero = 0` as shown below.
#' - Scale is `0:N` and 0 *is* a hurdle (e.g., a mix of a genuine/excess "zero"
#'   response with an ordinal `1:N` scale): keep the data as-is, use `vint(N)`
#'   (i.e., the total number of categories *minus* the hurdle category), and
#'   let `pzero` be estimated (optionally with predictors, `pzero ~ x`).
#' - Scale is `1:N` with excess responses piling up at the low end (e.g., a
#'   floor effect at the lowest category): recode by subtracting 1
#'   (`1:N` becomes `0:(N - 1)`), then proceed as in the previous bullet,
#'   i.e., `vint(N - 1)` and estimate `pzero`.
#'
#' @examples
#' # Fitting with brms. Because `k` is fixed data rather than a distributional
#' # parameter, it is passed through the brms::vint() addition term. Put the
#' # family on the formula, and cogmod_stanvars() supplies the Stan code for it.
#' f <- brms::bf(rating | vint(k) ~ predictor, family = cogmod_betadiscrete())
#' cogmod_stanvars(f)
#'
#' # To also model the hurdle probability (e.g., proportion of zero ratings):
#' brms::bf(rating | vint(k) ~ predictor, pzero ~ predictor,
#'   family = cogmod_betadiscrete()
#' )
#'
#' # To fix pzero at exactly 0, e.g. because your scale has no hurdle:
#' brms::bf(rating | vint(k) ~ predictor, pzero = 0,
#'   family = cogmod_betadiscrete()
#' )
#'
#' @export
cogmod_betadiscrete <- function(
  link_mu = "logit",
  link_phi = "log",
  link_pzero = "logit"
) {
  brms::custom_family(
    name = "cogmod_betadiscrete",
    dpars = c("mu", "phi", "pzero"),
    links = c(link_mu, link_phi, link_pzero),
    lb = c(NA, 0, 0), # phi > 0; 0 <= pzero <= 1
    ub = c(NA, NA, 1),
    type = "int", # discrete outcome
    vars = "vint1[n]" # k, passed via vint()
  )
}

# brms Post-processing Functions -------------------------------------------

# Subsets a dpar that may be a draws x obs matrix (estimated) or a
# draws-length vector (fixed, e.g. via `pzero = 0` in `brms::bf()`, or any
# other intercept-only dpar), into a draws x length(idx) matrix.
#' @keywords internal
.subset_dpar <- function(x, idx, n_draws) {
  if (is.matrix(x)) {
    return(x[, idx, drop = FALSE])
  }
  matrix(x, nrow = n_draws, ncol = length(idx))
}

#' @rdname rcogmod_betadiscrete
#' @export
log_lik_cogmod_betadiscrete <- function(i, prep) {
  y <- prep$data$Y[i]
  k <- prep$data$vint1[i]

  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pzero <- brms::get_dpar(prep, "pzero", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(numeric(0))
  }

  y_vec <- rep(y, length.out = n_draws)

  ll <- dcogmod_betadiscrete(
    x = y_vec,
    mu = mu,
    phi = phi,
    k = k,
    pzero = pzero,
    log = TRUE
  )
  ll[is.nan(ll) | is.na(ll)] <- -Inf
  ll
}


#' @rdname rcogmod_betadiscrete
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_cogmod_betadiscrete <- function(i, prep, ...) {
  k <- prep$data$vint1[i]

  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pzero <- brms::get_dpar(prep, "pzero", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(matrix(numeric(0), ncol = 1))
  }

  final_out <- rcogmod_betadiscrete(
    n = n_draws,
    mu = mu,
    phi = phi,
    k = k,
    pzero = pzero
  )
  as.matrix(final_out)
}


#' @rdname rcogmod_betadiscrete
#' @export
posterior_epred_cogmod_betadiscrete <- function(prep) {
  mu <- brms::get_dpar(prep, "mu") # draws x obs, or draws-length vector if no predictors
  phi <- brms::get_dpar(prep, "phi") # draws x obs, or draws-length vector if no predictors
  pzero <- brms::get_dpar(prep, "pzero") # draws x obs, or draws-length vector if fixed/intercept-only
  k <- prep$data$vint1 # length = obs

  n_draws <- prep$ndraws
  n_obs <- length(k)
  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  # E[R | alpha, beta] = k - sum_{j=1}^{k-1} F_B(j/k; alpha, beta)  (eq. 4)
  # Grouping by unique k avoids redundant computation when (as is typical)
  # all observations share the same number of rating categories. The hurdle
  # simply rescales the non-zero part of the expectation by (1 - pzero).
  for (k_val in unique(k)) {
    idx <- which(k == k_val)
    mu_sub <- .subset_dpar(mu, idx, n_draws)
    phi_sub <- .subset_dpar(phi, idx, n_draws)
    pzero_sub <- .subset_dpar(pzero, idx, n_draws)

    if (k_val <= 1) {
      epred[, idx] <- k_val * (1 - pzero_sub)
      next
    }

    alpha_sub <- mu_sub * phi_sub * 2
    beta_sub <- (1 - mu_sub) * phi_sub * 2

    cum_sum <- matrix(0, nrow = n_draws, ncol = length(idx))
    for (j in seq_len(k_val - 1)) {
      cum_sum <- cum_sum + stats::pbeta(j / k_val, alpha_sub, beta_sub)
    }

    epred[, idx] <- (1 - pzero_sub) * (k_val - cum_sum)
  }

  epred
}
