#' @title Beta-Gate Model
#'
#' @description
#' The Beta-Gate model represents subjective ratings as a mixture of a continuous Beta distribution
#' with additional point masses at the extremes (0 and 1). This structure effectively captures
#' common patterns in subjective rating data where respondents often select extreme values
#' at higher rates than would be expected from a Beta distribution alone.
#'
#' The Beta-Gate model corresponds to a reparametrized ordered beta model
#' (Kubinec, 2023, \doi{10.1017/pan.2022.20}).
#' In the ordered Beta model, the extreme values (0 and 1) arise from censoring an underlying
#' latent process based on cutpoints ("gates"). Values falling past the gates are considered extremes
#' (zeros and ones). The difference from the Ordered Beta is the way the cutpoints are defined,
#' as well as the scale of the precision parameter phi.
#'
#' It differs from the Zero-One-Inflated Beta (ZOIB) model in that the ZOIB model has `zoi`
#' and `coi` parameters, directly controlling the likelihood of extreme values. Instead,
#' Beta-Gate uses `pex` and `bex` to define "cutpoints" after which extreme values become likely.
#' In an ordered beta framework, the boundary probabilities arise through a single underlying
#' ordering process (the location of the cutpoints on the latent scale). In a ZOIB framework,
#' the boundaries are more like additional mass points inserted into a beta distribution.
#' In Beta-gate models, extreme values arise naturally from thresholding a single latent process.
#'
#' @param n Number of simulated values.
#' @param mu Mean of the underlying Beta distribution (`0 < mu < 1`).
#' @param phi Precision parameter of the underlying Beta distribution (must be strictly positive).
#'   Can be conceptualized as an "agreement" indicator: higher `phi` means less
#'   dispersion (more agreement) among ratings, holding `mu` fixed.
#'   Note: In many implementations, `phi` is parametrized differently, and correspond to
#'   the double of our `phi` argument (cogmod's `phi` = standard's `phi` * 2). Our parametrization
#'   Makes it `phi = 1` corresponds to uniform when `mu = 0.5`, which makes setting priors more convenient
#'   (e.g., on the logit scale)
#' @param pex Controls the location of the lower and upper boundary gates (`0 <= pex <= 1`). It defines
#'   the total probability mass allocated to the extremes (0 or 1).  Higher `pex` increases the probability
#'   of extreme values (0 or 1).
#' @param bex Balances the extreme probability mass `pex` between 0 and 1 (`0 <= bex <= 1`). A balance
#'   of `0.5` means that the 'gates' are symmetrically placed around the center of the distribution, and
#'   values higher or lower than `0.5` will shift the relative "ease" of crossing the gates towards 1
#'   or 0, respectively.
#'
#' @details
#' **Special cases:**
#' - When `pex = 0`: Pure Beta distribution with mean `mu` and precision `phi * 2`.
#' - When `pex = 1`: Pure Bernoulli distribution with `P(1) = bex`, `P(0) = 1-bex`.
#' - When `bex = 0` and `pex = 1`: All mass at 0.
#' - When `bex = 1` and `pex = 1`: All mass at 1.
#'
#' **Psychological Interpretation:**
#' - `mu`: Can be interpreted as the underlying average tendency or preference strength,
#'   disregarding extreme "all-or-nothing" responses.
#' - `phi`: Reflects the certainty or consistency of the non-extreme responses. Higher `phi`
#'   indicates responses tightly clustered around `mu` (more certainty), while lower `phi`
#'   (especially `phi = 1`) suggests more uniform or uncertain responses.
#' - `pex`: Represents the overall tendency towards extreme responding (choosing 0 or 1).
#'   This could reflect individual response styles (e.g., acquiescence, yea-saying/nay-saying)
#'   or properties of the item itself (e.g., polarizing questions).
#' - `bex`: Indicates the *direction* of the extreme response bias. `bex > 0.5` suggests a bias
#'   for producing ones more easily, while `bex < 0.5` suggests a bias towards zero.
#'
#' @return `rcogmod_betagate()` returns a numeric vector of `n` simulated
#'   ratings on the unit interval `[0, 1]`, including the exact `0`s and `1`s
#'   produced by the gates. `dcogmod_betagate()` returns the density at each
#'   element of `x` - the log density if `log = TRUE` - recycled to the length
#'   of the longest argument; at `0` and `1` it is the probability mass rather
#'   than a density. `cogmod_betagate()` returns a `brms::custom_family`
#'   object, to put on a `brms::bf()` formula. `cogmod_betagate_stanvars()`
#'   returns a `brms::stanvars` object holding the family's Stan `functions`
#'   block, to pass to `brms::brm()`, and `cogmod_betagate_lpdf_expose()`
#'   compiles that Stan code and returns it as an R function, for checking the
#'   density outside of a model. The remaining functions are `brms`
#'   post-processing methods, called by `brms` rather than directly:
#'   `log_lik_cogmod_betagate()` returns a numeric vector holding one
#'   log-likelihood value per posterior draw for observation `i`,
#'   `posterior_predict_cogmod_betagate()` a draws x 1 matrix of ratings
#'   simulated for observation `i`, and `posterior_epred_cogmod_betagate()` a
#'   draws x observations matrix of expected ratings.
#'
#' @references
#' - Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with
#'     lower and upper bounds. Political Analysis, 31(4), 519-536.
#'
#' @examples
#' # Symmetric gates (c0=0.05, c1=0.95), pex=0.1, bex=0.5
#' x1 <- rcogmod_betagate(10000, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5)
#' hist(x1, breaks=50, main="rcogmod_betagate: Symmetric Cutpoints (pex=0.1)")
#'
#' # Asymmetric gates (c0=0.15, c1=0.95), pex=0.2, bex=0.25
#' x2 <- rcogmod_betagate(10000, mu = 0.5, phi = 3, pex = 0.2, bex = 0.25)
#' hist(x2, breaks=50, main="rcogmod_betagate: Asymmetric Cutpoints (pex=0.2, bex=0.25)")
#'
#' # No gating (pure Beta)
#' x3 <- rcogmod_betagate(10000, mu = 0.7, phi = 5, pex = 0, bex = 0.5)
#' hist(x3, breaks=50, main="rcogmod_betagate: No Extreme Values (pex=0)")
#'
#' @export
rcogmod_betagate <- function(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5) {
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
  if (any(pex < 0 | pex > 1)) {
    stop("pex must be between 0 and 1.")
  }
  if (any(bex < 0 | bex > 1)) {
    stop("bex must be between 0 and 1.")
  }

  # --- Vectorization ---
  n_out <- max(n, length(mu), length(phi), length(pex), length(bex))
  if (n_out > 1 || n > 1) {
    # Ensure vectorization if n>1 even if params are scalar
    mu <- rep(mu, length.out = n_out)
    phi <- rep(phi, length.out = n_out)
    pex <- rep(pex, length.out = n_out)
    bex <- rep(bex, length.out = n_out)
  } else {
    n_out <- n # Case where n=1 and params are scalar
  }

  # --- Parameter Calculation ---
  eps <- .Machine$double.eps # Smallest representable positive number

  # Cutpoints on probability scale
  cutzero <- pex * (1 - bex)
  cutone <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero)
  cutonelog <- stats::qlogis(cutone)

  # Beta distribution parameters
  shape1 <- mu * phi * 2
  shape2 <- (1 - mu) * phi * 2

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # --- Probabilities for Outcome Categories ---
  # P(outcome = 0) = P(eta < cutzerolog) = P(logistic(mu_ql) < cutzero)
  # P(outcome = 1) = P(eta > cutonelog) = P(logistic(mu_ql) > cutone)
  # P(0 < outcome < 1) = P(cutzerolog <= eta <= cutonelog)
  # Using plogis(q, location) = P(X <= q) where X ~ Logistic(location, scale=1)
  # P(eta < cutzerolog) = plogis(cutzerolog, location = mu_ql) = 1 - plogis(mu_ql - cutzerolog)
  # P(eta > cutonelog) = 1 - plogis(cutonelog, location = mu_ql) = plogis(mu_ql - cutonelog)
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE) # P(eta <= cutzerolog)
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE) # P(eta > cutonelog)
  prob_mid <- 1 - prob_0 - prob_1
  # Ensure probabilities sum to 1 and handle potential floating point inaccuracies
  prob_mid <- pmax(0, prob_mid)
  probs_matrix <- cbind(prob_0, prob_mid, prob_1)
  probs_matrix <- probs_matrix / rowSums(probs_matrix) # Normalize row-wise

  # --- Simulation ---
  # Sample outcome category (0, middle, 1) for each observation
  # Uses Gumbel-max trick for efficient multinomial sampling
  gumbel_noise <- -log(
    -log(matrix(stats::runif(n_out * 3), nrow = n_out, ncol = 3))
  )
  outcome_category <- max.col(log(probs_matrix) + gumbel_noise) # 1 for 0, 2 for middle, 3 for 1

  # Generate underlying Beta draws for all (only used for middle category)
  beta_draws <- stats::rbeta(n = n_out, shape1 = shape1, shape2 = shape2)
  # Clamp beta draws slightly away from exact 0/1 for stability if needed downstream
  beta_draws <- pmax(eps, pmin(1 - eps, beta_draws))

  # --- Combine Results ---
  final_out <- numeric(n_out)
  final_out[outcome_category == 1] <- 0
  final_out[outcome_category == 3] <- 1
  mid_indices <- which(outcome_category == 2)
  if (length(mid_indices) > 0) {
    final_out[mid_indices] <- beta_draws[mid_indices]
  }

  # Ensure output length matches original n if n was the max dimension
  if (n_out == n && n > 0) {
    return(final_out[1:n])
  } else if (n_out > n && n == 1) {
    return(final_out) # Return vector if n=1 but params were vectors
  } else {
    return(final_out) # Default case
  }
}


#' @rdname rcogmod_betagate
#' @param x Vector of quantiles (values at which to evaluate the density). Must be between 0 and 1, inclusive.
#' @param log Logical; if TRUE, returns the log-density.
#' @examples
#' x <- seq(0, 1, length.out = 1001)
#' densities <- dcogmod_betagate(x, mu = 0.5, phi = 5, pex = 0.2, bex = 0.5)
#' plot(x, densities, type = "l", main = "Density Function", xlab = "y", ylab = "Density")
#' @export
dcogmod_betagate <- function(x, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5, log = FALSE) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(pex < 0 | pex > 1)) {
    stop("pex must be between 0 and 1.")
  }
  if (any(bex < 0 | bex > 1)) {
    stop("bex must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(x)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  pex <- rep(pex, length.out = n)
  bex <- rep(bex, length.out = n)

  # --- Parameter Calculation ---
  eps <- .Machine$double.eps # Smallest representable positive number

  # Cutpoints on probability scale
  cutzero <- pex * (1 - bex)
  cutone <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero)
  cutonelog <- stats::qlogis(cutone)

  # Beta distribution parameters
  shape1 <- mu * phi * 2
  shape2 <- (1 - mu) * phi * 2

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # --- Probabilities for Outcome Categories ---
  # See rcogmod_betagate for derivation
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE)
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE)
  prob_mid <- 1 - prob_0 - prob_1
  prob_mid <- pmax(0, prob_mid) # Handle potential floating point inaccuracies

  # --- Calculate Density ---
  # Initialize density vector
  density <- numeric(n)

  # Indices for different cases
  idx_zero <- which(x == 0)
  idx_one <- which(x == 1)
  idx_mid <- which(x > 0 & x < 1)
  idx_outside <- which(x < 0 | x > 1) # Should have 0 density

  # Density for x = 0
  if (length(idx_zero) > 0) {
    density[idx_zero] <- prob_0[idx_zero]
  }

  # Density for x = 1
  if (length(idx_one) > 0) {
    density[idx_one] <- prob_1[idx_one]
  }

  # Density for 0 < x < 1
  if (length(idx_mid) > 0) {
    # Calculate Beta density component
    # Need to handle case where prob_mid is zero (e.g., pex=1)
    beta_dens <- ifelse(
      prob_mid[idx_mid] > eps,
      stats::dbeta(
        x[idx_mid],
        shape1 = shape1[idx_mid],
        shape2 = shape2[idx_mid],
        log = FALSE
      ),
      0
    )
    # Total density is P(middle category) * BetaPDF(x | params)
    density[idx_mid] <- prob_mid[idx_mid] * beta_dens
  }

  # Density for x outside [0, 1] is 0 (already initialized)

  # --- Return Log Density if Requested ---
  if (log) {
    # Use log1p for potentially small probabilities near 0? No, log(prob) is fine.
    # Handle density = 0 case -> log(0) = -Inf
    density <- ifelse(density > 0, log(density), -Inf)
  }

  density
}


# Stan Functions and Custom Family for brms (Beta-Gate)

# Stanvars ----------------------------------------------------------------

#' @keywords internal
.cogmod_betagate_lpdf <- function() {
  "
// Log probability density function for the Beta-Gate distribution
real cogmod_betagate_lpdf(real y, real mu, real phi, real pex, real bex) {
  // Tolerance for floating point comparisons near 0 and 1
  real eps = 1e-10;

  // --- Parameter Validation ---
  // Ensure parameters are within valid ranges. Note: brms often handles this
  // via link functions, but explicit checks add robustness.
  if (!(mu > 0.0 && mu < 1.0) || !(phi > 0.0) ||
      !(pex >= 0.0 && pex <= 1.0) || !(bex >= 0.0 && bex <= 1.0) ||
      !(y >= 0.0 && y <= 1.0)) {
    return negative_infinity();
  }

  // --- Calculate Cutpoints ---
  // Calculate probability-scale cutpoints and apply logit scale
  real cutzerolog = logit(pex * (1.0 - bex));
  real cutonelog = logit(1.0 - pex * bex);

  // --- Calculate Log Probability based on y ---
  // Location parameter on logit scale
  real mu_ql = logit(mu);

  if (abs(y - 0.0) < eps) { // Case: y = 0
    // Log probability P(latent <= cutzerolog)
    return log1m_inv_logit(mu_ql - cutzerolog);

  } else if (abs(y - 1.0) < eps) { // Case: y = 1
    // Log probability P(latent > cutonelog)
    return log_inv_logit(mu_ql - cutonelog);

  } else { // Case: 0 < y < 1
    // Log probability P(cutzerolog < latent <= cutonelog)
    real log_prob_middle = log_diff_exp(log_inv_logit(mu_ql - cutzerolog),
                                        log_inv_logit(mu_ql - cutonelog));

    // Beta distribution parameters
    real shape1 = mu * phi * 2.0;
    real shape2 = (1.0 - mu) * phi * 2.0;

    // Log Beta density for y
    real log_beta_dens = beta_lpdf(y | shape1, shape2);

    // Total log probability: log(P(middle)) + log(BetaPDF(y))
    return log_prob_middle + log_beta_dens;
  }
}
"
}


#' @rdname rcogmod_betagate
#' @examples
#' \donttest{
#' # Exposing the Stan function needs cmdstanr and a CmdStan toolchain,
#' # which live outside CRAN - see the package website to install them.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   lpdf <- cogmod_betagate_lpdf_expose()
#'   lpdf(y = 0.5, mu = 0.6, phi = 10, pex = 0.2, bex = 0.5)
#' }
#' }
#'
#' @export
cogmod_betagate_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .cogmod_betagate_lpdf(),
    "\n}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$cogmod_betagate_lpdf
}

#' @rdname rcogmod_betagate
#' @export
cogmod_betagate_stanvars <- function() {
  brms::stanvar(scode = .cogmod_betagate_lpdf(), block = "functions")
}


#' @rdname rcogmod_betagate
#' @param link_mu,link_phi,link_pex,link_bex Link functions for the parameters.
#' @export
cogmod_betagate <- function(link_mu = "logit", link_phi = "softplus", link_pex = "logit", link_bex = "logit") {
  brms::custom_family(
    name = "cogmod_betagate",
    dpars = c("mu", "phi", "pex", "bex"),
    links = c(link_mu, link_phi, link_pex, link_bex),
    lb = c(NA, 0, 0, 0), # Lower bounds: phi>0, pex>=0, bex>=0
    ub = c(NA, NA, 1, 1), # Upper bounds: pex<=1, bex<=1
    type = "real" # Outcome variable type
  )
}

# brms Post-processing Functions -------------------------------------------

#' @rdname rcogmod_betagate
#' @export
log_lik_cogmod_betagate <- function(i, prep) {
  # Extract observed value for the i-th observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y_scalar <- prep$data$Y[i] # This is a single value

  # Extract model draws (vectors) for the i-th observation
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  # Determine number of draws
  n_draws <- length(mu)
  if (n_draws == 0) {
    return(numeric(0))
  } # Handle case with no draws

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y_scalar, length.out = n_draws)

  # Calculate log-likelihood using the vectorized dcogmod_betagate function
  ll <- dcogmod_betagate(x = y_vec, mu = mu, phi = phi, pex = pex, bex = bex, log = TRUE)

  # Ensure no NaN/NA values (dcogmod_betagate should return -Inf for zero density)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll # Return the vector of log-likelihoods for all draws
}


#' @rdname rcogmod_betagate
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_cogmod_betagate <- function(i, prep, ...) {
  # Extract draws for each parameter for the i-th observation
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  n_draws <- length(mu) # Number of posterior draws
  if (n_draws == 0) {
    return(matrix(numeric(0), ncol = 1))
  } # Handle case with no draws

  # Use rcogmod_betagate() to generate predictions for each draw
  final_out <- rcogmod_betagate(n = n_draws, mu = mu, phi = phi, pex = pex, bex = bex)

  # Return as a matrix with one column (standard for posterior_predict)
  as.matrix(final_out)
}


#' @rdname rcogmod_betagate
#' @export
posterior_epred_cogmod_betagate <- function(prep) {
  # Extract draws for the necessary parameters (draws x observations matrices)
  mu <- brms::get_dpar(prep, "mu")
  phi <- brms::get_dpar(prep, "phi") # phi is not directly needed for epred, but kept for consistency
  pex <- brms::get_dpar(prep, "pex")
  bex <- brms::get_dpar(prep, "bex")

  # Calculate the expectation (mean) based on the Ordered Beta logic
  # E[Y] = P(Y=1) * 1 + P(0 < Y < 1) * E[Y | 0 < Y < 1]
  # E[Y | 0 < Y < 1] is simply the mean of the Beta distribution, which is mu.
  # E[Y] = P(Y=1) + P(0 < Y < 1) * mu

  # --- Calculate Probabilities (Vectorized) ---
  # Cutpoints on probability scale
  cutzero_p <- pex * (1 - bex)
  cutone_p <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero_p)
  cutonelog <- stats::qlogis(cutone_p)

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # Probabilities for the three categories
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE) # P(eta <= cutzerolog) -> Y=0
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE) # P(eta > cutonelog) -> Y=1
  prob_mid <- 1 - prob_0 - prob_1 # P(cutzerolog < eta <= cutonelog) -> Y ~ Beta(mu, phi)
  # Ensure probabilities sum to 1 and handle potential floating point inaccuracies
  prob_mid <- pmax(0, prob_mid)

  # --- Calculate Expectation ---
  # E[Y] = P(Y=1) * 1 + P(0 < Y < 1) * E[Y | 0 < Y < 1]
  # E[Y] = prob_1 * 1 + prob_mid * mu
  epred <- prob_1 + prob_mid * mu

  # Handle edge cases where parameters might lead to NaN (e.g., invalid mu)
  epred[is.nan(epred) | is.na(epred)] <- NA # Or a default value if preferred

  epred # Return the matrix of posterior expectations (draws x observations)
}
