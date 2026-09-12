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
#'   length of the longest argument. `cogmod_exgaussian()` returns a
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


# Internals ---------------------------------------------------------------

#' @keywords internal
.prepare_exgaussian <- function(x = NULL, n = NULL, mu, sigma, tau) {
    # Validate parameters once. `mu` is deliberately NOT checked: it is the
    # location of the Gaussian component, the convolution is defined for any
    # real value, and the Stan lpdf has only ever rejected the two scales. The
    # R side used to reject `mu <= 0`, which meant it refused inputs the sampler
    # was perfectly happy to fit. See ?rcogmod_exgaussian.
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

#' @keywords internal
.cogmod_exgaussian_lpdf <- function() {
    "
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

    // Stan's built-in exp_mod_normal is parameterized with the rate of the
    // exponential component (beta = 1 / tau)
    return exp_mod_normal_lpdf(Y | mu, sigma, inv(tau));
}
"
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
    # exp_mod_normal_lpdf(y | mu, sigma, 1 / tau).
    ll <- -log(tau) +
        (mu / tau) +
        (sigma^2) / (2 * tau^2) -
        (y_vec / tau) +
        stats::pnorm(
            y_vec,
            mean = mu + (sigma^2) / tau,
            sd = sigma,
            log.p = TRUE
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
