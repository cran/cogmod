#' @title Ex-Gaussian Model (Classical Parameterization)
#'
#' @description
#' Density, random generation, and `brms` custom family for the Ex-Gaussian
#' distribution, using the "classical" parameterization familiar to
#' experimental psychologists, in which `mu` and `sigma` are the mean and SD
#' of the Gaussian component alone, and `tau` is the mean of the exponential
#' component (the tail). This is unlike `brms`'s built-in `exgaussian()`
#' family, in which `mu` indexes the mean of the *entire* distribution
#' (Gaussian + exponential components combined).
#'
#' Functions:
#' - `rcogmod_exgaussian()`: Simulates random draws from the Ex-Gaussian distribution.
#' - `dcogmod_exgaussian()`: Computes the density (likelihood) of the Ex-Gaussian distribution.
#' - `pcogmod_exgaussian()`: Computes the cumulative distribution function (CDF) or survival.
#' - `cogmod_exgaussian()`: Creates a `brms::custom_family()` for use in `brms` models.
#'
#' @details
#' The Ex-Gaussian distribution is the sum of an independent Normal (Gaussian)
#' random variable with mean `mu` and SD `sigma`, and an Exponential random
#' variable with mean `tau` (rate `1 / tau`). Unlike `brms`'s built-in
#' `exgaussian()` family - in which `mu` indexes the mean of the *entire*
#' distribution (Gaussian + exponential components combined) - here `mu` is
#' the mean of the Gaussian component alone, so that the mean of the full
#' distribution is `mu + tau`.
#'
#' This distinction matters because changes in the Gaussian location (`mu`)
#' and changes in the exponential tail (`tau`) can offset one another at the
#' level of the overall mean, so effects estimated on `brms`'s default `mu`
#' can lead to different (and potentially incorrect) inferences than effects
#' estimated on this classical `mu`.
#'
#' In the `brms` custom family (`cogmod_exgaussian()`), `sigma` and `tau` use a
#' `"softplus"` link (`log(1 + exp(x))`) rather than `"log"`. Both are scales
#' and must be strictly positive for the density to exist at all. A `"log"` link
#' would enforce that too, but its curvature explodes as the linear predictor
#' departs from zero, producing extreme gradients and making priors and sampling
#' harder to calibrate - and `tau` is on the RT scale (seconds), where it can
#' take comparatively large values. `"softplus"` is positive-constrained like
#' `"log"` but behaves almost linearly (`softplus(x) ~ x`) away from zero, so
#' weakly-informative priors can be stated directly on the RT scale.
#'
#' `mu` is different, and uses `"identity"`. It is the **location** of the
#' Gaussian component, not a scale: the convolution is well defined for any real
#' value, and the density integrates to one at `mu = 0` or below just as it does
#' above (the Stan `lpdf` has always accepted a non-positive `mu` - it checks
#' only `sigma` and `tau`). Nothing is gained by constraining it, and two things
#' are lost. First, interpretability, which is most of the reason to prefer the
#' ex-Gaussian in the first place: behind a `softplus` link a coefficient is not
#' in seconds, and the conversion factor moves with the intercept - the local
#' slope is 0.33 at `mu = 0.4` s, 0.39 at 0.5 s and 0.63 at 1 s, so the same
#' effect reads as a different number depending on where the intercept sits. On
#' `"identity"` a coefficient is seconds, full stop. Second, fidelity: for fast,
#' heavily-tailed data the Gaussian component genuinely belongs near or below
#' zero with `tau` carrying the mass, and forcing `mu > 0` distorts the `mu`/`tau`
#' split in exactly the cases where that decomposition is the thing being
#' estimated. This also matches every other implementation - `brms`'s own
#' `exgaussian()`, `retimes`, and the estimates reported in the literature - so
#' fitted values are directly comparable.
#'
#' The cost is that `brms`'s default intercept prior is no longer sensible for
#' `mu` (`student_t(3, 0, 2.5)` centred at zero seconds), which is why
#' [cogmod_priors()] now supplies one.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mu Mean of the Gaussian component. Unbounded - it is a location, not a
#'   scale - though for RT data it is normally positive. Range: (-Inf, Inf).
#' @param sigma SD of the Gaussian component. Must be positive. Range: (0, Inf).
#' @param tau Mean of the exponential component (the tail). Must be positive. Range: (0, Inf).
#'
#' @references
#' - Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian
#'     and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*,
#'     *16*(5), 798-817. \doi{10.3758/PBR.16.5.798}
#'
#' @return `rcogmod_exgaussian()` returns a numeric vector of `n` simulated
#'   reaction times, in seconds. `dcogmod_exgaussian()` returns the density at
#'   each element of `x` - the log density if `log = TRUE` - recycled to the
#'   length of the longest argument. `pcogmod_exgaussian()` returns the
#'   cumulative probability at each element of `q`, honouring `lower.tail` and
#'   `log.p`. `cogmod_exgaussian()` returns a
#'   `brms::custom_family` object, to put on a `brms::bf()` formula.
#'   `cogmod_exgaussian_stanvars()` returns a `brms::stanvars` object holding
#'   the family's Stan `functions` block, to pass to `brms::brm()`, and
#'   `cogmod_exgaussian_lpdf_expose()` compiles that Stan code and returns it
#'   as an R function, for checking the density outside of a model. The
#'   remaining functions are `brms` post-processing methods, called by `brms`
#'   rather than directly: `log_lik_cogmod_exgaussian()` returns a numeric
#'   vector holding one log-likelihood value per posterior draw for
#'   observation `i`, `posterior_predict_cogmod_exgaussian()` a draws x 1
#'   matrix of reaction times simulated for observation `i`, and
#'   `posterior_epred_cogmod_exgaussian()` a draws x observations matrix of
#'   expected reaction times.
#'
#' @examples
#' # Simulate 1000 RTs
#' rts <- rcogmod_exgaussian(1000, mu = 0.5, sigma = 0.1, tau = 0.2)
#' hist(rts, breaks = 50, main = "Simulated Ex-Gaussian RTs", xlab = "Reaction Time")
#'
#' @export
rcogmod_exgaussian <- function(n, mu = 0.5, sigma = 0.1, tau = 0.2) {
    # Prepare and validate all inputs for RNG
    params <- .prepare_exgaussian(
        x = NULL,
        n = n,
        mu = mu,
        sigma = sigma,
        tau = tau
    )

    stats::rnorm(params$ndraws, mean = params$mu, sd = params$sigma) +
        stats::rexp(params$ndraws, rate = 1 / params$tau)
}


#' @rdname rcogmod_exgaussian
#' @param x Vector of quantiles (observed reaction times).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @export
dcogmod_exgaussian <- function(x, mu = 0.5, sigma = 0.1, tau = 0.2, log = FALSE) {
    # Prepare and validate inputs for density
    params <- .prepare_exgaussian(
        x = x,
        n = NULL,
        mu = mu,
        sigma = sigma,
        tau = tau
    )

    # Log-density of the classical Ex-Gaussian (Gaussian convolved with exponential),
    # computed on the log scale for numerical stability. Equivalent to Stan's
    # exp_mod_normal_lpdf(x | mu, sigma, 1 / tau).
    log_density <- -log(params$tau) +
        (params$mu / params$tau) +
        (params$sigma^2) / (2 * params$tau^2) -
        (params$x / params$tau) +
        stats::pnorm(
            params$x,
            mean = params$mu + (params$sigma^2) / params$tau,
            sd = params$sigma,
            log.p = TRUE
        )

    # Guard against numerical over/underflow
    log_density[is.nan(log_density)] <- -Inf

    if (log) {
        log_density
    } else {
        exp(log_density)
    }
}


#' @rdname rcogmod_exgaussian
#' @param q Vector of quantiles (reaction times, in seconds).
#' @param lower.tail Logical; if TRUE (default), probabilities are `P[X <= q]`,
#'   otherwise `P[X > q]` - the survival, which is what a right-censored
#'   response contributes to the likelihood under `brms::bf(rt | cens(x) ~ ...)`.
#'   See the *Censoring* section of [rcogmod_invgaussian()].
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @export
pcogmod_exgaussian <- function(q, mu = 0.5, sigma = 0.1, tau = 0.2,
                               lower.tail = TRUE, log.p = FALSE) {
    params <- .prepare_exgaussian(x = q, n = NULL, mu = mu, sigma = sigma,
                                  tau = tau)
    lp <- if (lower.tail) {
        .lcdf_exgaussian(params$x, params$mu, params$sigma, params$tau)
    } else {
        .lsurv_exgaussian(params$x, params$mu, params$sigma, params$tau)
    }
    if (log.p) lp else exp(lp)
}


# Internals ---------------------------------------------------------------

# log survival of the ex-Gaussian, as the sum of its two positive terms
#
#   S(x) = Phi(-z) + exp(sigma^2 / (2 tau^2) - (x - mu) / tau) * Phi(z - sigma / tau)
#
# rather than log(1 - F). The counterpart of .lcdf_exgaussian() (model_geg.R),
# which takes the left tail through the same two terms subtracted; each one is
# written in the tail where its terms do not cancel. Mirrors the Stan
# cogmod_exgaussian_logsurv() (.EXGAUSSIAN_STAN_PRELUDE below).
#' @keywords internal
.lsurv_exgaussian <- function(x, mu, sigma, tau) {
    z <- (x - mu) / sigma
    la <- stats::pnorm(-z, log.p = TRUE)
    lb <- sigma^2 / (2 * tau^2) - (x - mu) / tau +
        stats::pnorm(z - sigma / tau, log.p = TRUE)
    m <- pmax(la, lb)
    out <- m + log(exp(la - m) + exp(lb - m))
    out[!is.finite(m)] <- -Inf
    pmin(out, 0)
}

#' @keywords internal
.prepare_exgaussian <- function(x = NULL, n = NULL, mu, sigma, tau) {
    # Validate parameters once. `mu` is deliberately NOT checked: it is the
    # location of the Gaussian component, the convolution is defined for any
    # real value, and the Stan lpdf rejects only the two scales. Rejecting
    # `mu <= 0` here would refuse inputs the sampler is perfectly happy to fit.
    # See ?rcogmod_exgaussian.
    if (any(sigma <= 0)) {
        stop("`sigma` must be positive.")
    }
    if (any(tau <= 0)) {
        stop("`tau` must be positive.")
    }

    # Determine target length:
    if (!is.null(x)) {
        # For PDF: based on x length
        m <- length(x)
    } else if (!is.null(n)) {
        # For RNG: based on n draw count
        if (length(n) != 1 || n < 0 || n != floor(n)) {
            stop("n must be a single non-negative integer.")
        }
        m <- n
    } else {
        stop("Either 'x' or 'n' must be provided.")
    }

    # Recycle vectors to length m
    params <- list(
        mu = rep(mu, length.out = m),
        sigma = rep(sigma, length.out = m),
        tau = rep(tau, length.out = m)
    )

    # Incorporate x if provided
    if (!is.null(x)) {
        params$x <- rep(x, length.out = m)
    }
    params$ndraws <- m
    params
}


#' @rdname rcogmod_exgaussian
#' @param link_mu,link_sigma,link_tau Character of the type of link used to
#'   model the ex-Gaussian parameters. Defaults to `"identity"` for `mu` and
#'   `"softplus"` for `sigma` and `tau` (see Details).
#' @export
cogmod_exgaussian <- function(
    link_mu = "identity",
    link_sigma = "softplus",
    link_tau = "softplus"
) {
    brms::custom_family(
        name = "cogmod_exgaussian",
        dpars = c("mu", "sigma", "tau"), # mu/sigma = Gaussian location/SD, tau = exponential mean
        # mu is a LOCATION and is unbounded; only the two scales are positive.
        # See the Details section: the convolution is defined for any real mu,
        # and the Stan lpdf below has always accepted one.
        lb = c(NA, 0, 0),
        links = c(link_mu, link_sigma, link_tau),
        ub = c(NA, NA, NA),
        type = "real" # Continuous outcome variable (RT)
    )
}

# The ex-Gaussian's density and its two tails in Stan, the counterparts of
# dcogmod_exgaussian() and .lsurv_exgaussian() above and of .lcdf_exgaussian()
# (model_geg.R), line for line. Its own object because cogmod_geg() needs the
# CDF inside its *density* - the alpha-power construction is shape * log F -
# and so pastes this in front of its own functions too (see
# .cogmod_geg_lpdf()); helper-stan.R strips the second copy when it
# concatenates every family into one test program.
#
# All of it goes through cogmod_log_Phi() rather than Stan's
# exp_mod_normal_lcdf(). The built-in's value is fine - it agreed with the R
# side to 1e-8 over the grid test-model_geg.R walks - but its partials are not,
# and the GEG differentiates its log-likelihood through them. Measured
# 2026-09-18 at sigma / tau = 9.3, where the two terms of F cancel hardest:
# Stan's own central differences put d/d sigma at 28.5 and its autodiff at
# 77.5, with d/d tau 113.6 against 53.0. Subtracting the two terms ourselves,
# in log space, is both the accurate value and the right derivative, and it is
# the formula the R side has always used.
#
# The three are separate `real` functions rather than one returning the pair
# the GEG wants. That was tried: a Stan function returning `vector[2]`
# allocates one on the autodiff stack per call, and at one call per observation
# that outweighed everything the sharing saved - half again the cost of writing
# the pair out inline, measured against the built-ins at 1.47 with the vector
# and 1.20 without over the same seven blocks. cogmod_geg_lpdf() therefore
# writes the shared arithmetic out itself.
#' @keywords internal
.EXGAUSSIAN_STAN_PRELUDE <- paste0(.LOG_PHI_STAN_PRELUDE, "
// log f(y) for the ex-Gaussian. With z = (y - mu) / sigma this is
//
//   lb - log(tau),   lb = sigma^2 / (2 tau^2) - (y - mu) / tau
//                           + log Phi(z - sigma / tau),
//
// the same expression dcogmod_exgaussian() evaluates in R and the same one
// Stan's exp_mod_normal_lpdf() does, over cogmod_log_Phi() rather than the
// built-in's bare erfc: identical in the body of the distribution, and still
// finite about 38 standardized units into the left tail, where erfc is not.
real cogmod_exgaussian_ldens(real y, real mu, real sigma, real tau) {
  return square(sigma) / (2 * square(tau)) - (y - mu) / tau - log(tau)
           + cogmod_log_Phi((y - mu) / sigma - sigma / tau);
}

// log F(y) = log(Phi(z) - exp(lb)), with `lb` as above - note that the term F
// subtracts IS the density, up to the 1 / tau, which is what cogmod_geg_lpdf()
// exploits (see the note there).
//
// The difference is taken with log1m_exp() rather than by subtracting the two
// terms, because in the left tail they are individually tiny and very close
// together - exactly where cogmod_geg() needs the CDF, since shape < 1
// multiplies log F by a negative number and any error there is amplified.
// F > 0 everywhere, so the difference is positive; the floor on `d` is against
// rounding alone, it is the same one .lcdf_exgaussian() uses so that the two
// sides agree bit for bit rather than one of them returning -inf, and it is
// what leaves the sum negative without a closing fmin().
real cogmod_exgaussian_logcdf(real y, real mu, real sigma, real tau) {
  real z = (y - mu) / sigma;
  real la = cogmod_log_Phi(z);
  real lb = square(sigma) / (2 * square(tau)) - (y - mu) / tau
              + cogmod_log_Phi(z - sigma / tau);
  real d = fmin(lb - la, -2.220446049250313e-16);
  return la + log1m_exp(d);
}

// log S(y), as the SUM of the two positive terms
//   S(y) = Phi(-z) + exp(sigma^2 / (2 tau^2) - (y - mu) / tau) Phi(z - sigma / tau),
// rather than log1m_exp() of the above: a right-censored slow response sits
// exactly where 1 - F has no digits left. Each tail is written where its terms
// do not cancel.
real cogmod_exgaussian_logsurv(real y, real mu, real sigma, real tau) {
  real z = (y - mu) / sigma;
  return fmin(log_sum_exp(
    cogmod_log_Phi(-z),
    square(sigma) / (2 * square(tau)) - (y - mu) / tau
      + cogmod_log_Phi(z - sigma / tau)
  ), 0);
}
")


#' @keywords internal
.cogmod_exgaussian_lpdf <- function() {
    paste0(.EXGAUSSIAN_STAN_PRELUDE, "
// Log-likelihood for a single observation from the classical Ex-Gaussian distribution.
// Y: observed reaction time.
// mu: mean of the Gaussian component. A LOCATION, so unbounded - the
//     convolution is defined for any real value, which is why the check below
//     leaves it alone.
// sigma: SD of the Gaussian component (> 0).
// tau: mean of the exponential component, i.e., the tail (> 0).
real cogmod_exgaussian_lpdf(real Y, real mu, real sigma, real tau) {
    // Parameter checks
    if (sigma <= 0 || tau <= 0) return negative_infinity();

    // The same expression as Stan's exp_mod_normal_lpdf(Y | mu, sigma, 1 / tau)
    // and as dcogmod_exgaussian() in R, written through cogmod_log_Phi() rather
    // than the built-in's erfc: the value agrees to the last bit in the body of
    // the distribution and keeps going where erfc underflows, about 38
    // standardized units into the left tail. It also shares its one normal tail
    // with the CDF, which is what cogmod_geg() is built on.
    return cogmod_exgaussian_ldens(Y, mu, sigma, tau);
}

// CDF and survival of the same distribution, for brms's cens() addition term
// (see ?rcogmod_invgaussian for what censoring a reaction time means). Both are
// the helpers above; these wrappers exist only because brms writes
// `<family>_lcdf(y | ...)` and Stan reserves those two suffixes for the `|`
// call syntax.
real cogmod_exgaussian_lcdf(real Y, real mu, real sigma, real tau) {
    if (sigma <= 0 || tau <= 0) return negative_infinity();
    return cogmod_exgaussian_logcdf(Y, mu, sigma, tau);
}

real cogmod_exgaussian_lccdf(real Y, real mu, real sigma, real tau) {
    if (sigma <= 0 || tau <= 0) return negative_infinity();
    return cogmod_exgaussian_logsurv(Y, mu, sigma, tau);
}
")
}

#' @rdname rcogmod_exgaussian
#' @export
cogmod_exgaussian_lpdf_expose <- function() {
    insight::check_if_installed("cmdstanr")

    # Wrap the function Stan block
    stancode <- paste0(
        "functions {
",
        .cogmod_exgaussian_lpdf(),
        "
}"
    )

    mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
    mod$expose_functions()
    mod$functions$cogmod_exgaussian_lpdf
}

#' @rdname rcogmod_exgaussian
#' @export
cogmod_exgaussian_stanvars <- function() {
    brms::stanvar(scode = .cogmod_exgaussian_lpdf(), block = "functions")
}


# brms methods ------------------------------------------------------------

#' @rdname rcogmod_exgaussian
#' @inheritParams cogmod_lnr
#' @export
log_lik_cogmod_exgaussian <- function(i, prep) {
    # Extract observation
    if (!"Y" %in% names(prep$data)) {
        stop("Outcome variable 'Y' not found in prep$data.")
    }
    y <- prep$data$Y[i]
    if (is.na(y)) {
        return(NA_real_)
    }

    # Get parameters for observation i across all draws
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)

    # Determine number of draws
    n_draws <- length(mu)
    if (n_draws == 0) {
        return(numeric(0))
    }

    # Replicate the scalar y to match the number of draws
    y_vec <- rep(y, length.out = n_draws)

    # Log-density of the classical Ex-Gaussian (Gaussian convolved with exponential),
    # computed on the log scale for numerical stability. Equivalent to Stan's
    # exp_mod_normal_lpdf(y | mu, sigma, 1 / tau). Wrapped in .censor_ll() so
    # that a `cens()` term on the formula is honoured here the way it is in the
    # Stan program - brms leaves that to the custom family's own method.
    ll <- .censor_ll(
        prep, i, y,
        ldens = function(y) {
            -log(tau) +
                (mu / tau) +
                (sigma^2) / (2 * tau^2) -
                (y_vec / tau) +
                stats::pnorm(
                    y_vec,
                    mean = mu + (sigma^2) / tau,
                    sd = sigma,
                    log.p = TRUE
                )
        },
        lcdf = function(y, lower.tail) {
            yv <- rep(y, length.out = n_draws)
            if (lower.tail) {
                .lcdf_exgaussian(yv, mu, sigma, tau)
            } else {
                .lsurv_exgaussian(yv, mu, sigma, tau)
            }
        }
    )

    # Ensure correct handling for invalid parameters
    ll[sigma <= 0 | tau <= 0] <- -Inf

    # Ensure no other NaN/NA values
    ll[is.nan(ll) | is.na(ll)] <- -Inf

    ll
}

#' @rdname rcogmod_exgaussian
#' @inheritParams cogmod_lnr
#' @export
posterior_predict_cogmod_exgaussian <- function(i, prep, ...) {
    # Get parameters for observation i across all draws
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)

    # Number of posterior draws
    n_draws <- length(mu)

    # Simulate as Gaussian + Exponential (vectorized)
    stats::rnorm(n = n_draws, mean = mu, sd = sigma) +
        stats::rexp(n = n_draws, rate = 1 / tau)
}

#' @rdname rcogmod_exgaussian
#' @export
posterior_epred_cogmod_exgaussian <- function(prep) {
    # Extract draws for the necessary parameters (matrices: draws x observations)
    mu <- brms::get_dpar(prep, "mu")
    tau <- brms::get_dpar(prep, "tau")

    # E[ExGaussian] = E[Gaussian] + E[Exponential] = mu + tau
    mu + tau
}
