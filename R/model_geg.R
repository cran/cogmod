#' @title Generalised Ex-Gaussian (GEG) Distribution
#'
#' @description
#' The Generalised Ex-Gaussian of Marmolejo-Ramos et al. (2023): the ex-Gaussian
#' with one extra `shape` parameter, obtained by raising its CDF to a power.
#'
#' \deqn{F_{GEG}(x) = \left[F_{EG}(x)\right]^{shape}}
#'
#' so that the density is
#'
#' \deqn{f_{GEG}(x) = shape \cdot \left[F_{EG}(x)\right]^{shape-1} f_{EG}(x)}
#'
#' At `shape = 1` this is the ex-Gaussian exactly - not approximately - so
#' [cogmod_exgaussian()] is nested inside it and `loo_compare()` between the two
#' is like-for-like.
#'
#' @details
#' # What the shape parameter buys
#'
#' A wider range of shapes than the ex-Gaussian can reach. Sweeping `sigma` and
#' `tau` across the values RT data occupy, the ex-Gaussian spans skewness 0 to 2
#' and excess kurtosis 0 to 6; freeing `shape` widens that to roughly -0.4 to 4.8
#' and 0 to 35. In particular the GEG can be **negatively skewed**, which the
#' ex-Gaussian cannot be at any parameter value.
#'
#' # What it costs
#'
#' Interpretability, and specifically the one property that makes the
#' ex-Gaussian worth using as a *descriptive* model.
#'
#' - **The mean is no longer `mu + tau`.** With `mu = 0.4` and `tau = 0.2`, the
#'   mean runs 0.31 at `shape = 0.2`, 0.60 at `shape = 1` and 1.15 at
#'   `shape = 20`. There is no closed form for it - see
#'   [posterior_epred()][brms::posterior_epred] below - and the bulk/tail
#'   decomposition that the ex-Gaussian is normally reported for does not
#'   survive.
#' - **`shape` is strongly confounded with `mu`.** Fitted by maximum likelihood
#'   to the lexical-decision data used in the package vignettes, the correlation
#'   between the two at the optimum is about -0.98, and the other estimates move
#'   with it: on one condition `mu` goes 0.429 to 0.508, `sigma` 0.051 to 0.037
#'   and `tau` 0.119 to 0.162 once `shape` is freed. `shape` does not add an
#'   independent axis so much as re-slice the same bulk-and-tail split.
#'
#' The practical consequence is that `cogmod_priors()` gives `shape` a
#' deliberately informative prior centred on the ex-Gaussian, and that `mu`,
#' `sigma` and `tau` should not be read as the Gaussian centre, the Gaussian SD
#' and the mean of the tail once `shape` is free. If those quantities are the
#' point of the analysis, fit [cogmod_exgaussian()] instead. If a better-fitting
#' descriptive family is the point, [cogmod_logstudent()] and
#' [cogmod_loggamma()] decouple skew from tail weight with parameters that stay
#' interpretable.
#'
#' # Construction
#'
#' The power transform is Durrans' alpha-power (or "exponentiated") family, and
#' for integer `shape` it is the distribution of the maximum of `shape`
#' independent ex-Gaussian draws. That is a mathematical device rather than an
#' account of a process, so unlike [cogmod_invgaussian()]'s `sigmadrift` there is
#' no mechanism attached to it.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to be
#'   the number required.
#' @param mu Location of the Gaussian component. Unbounded. Range: (-Inf, Inf).
#' @param sigma SD of the Gaussian component. Must be positive. Range: (0, Inf).
#' @param tau Mean of the exponential component. Must be positive. Range: (0, Inf).
#' @param shape Power applied to the ex-Gaussian CDF. Must be positive.
#'   `shape = 1` gives the ex-Gaussian back exactly. Range: (0, Inf).
#'
#' @references
#' - Marmolejo-Ramos, F., Barrera-Causil, C., Kuang, S., Fazlali, Z., Wegener, D.,
#'     Kneib, T., De Bastiani, F., & Martinez-Florez, G. (2023). Generalised
#'     exponential-Gaussian distribution: A method for neural reaction time
#'     analysis. *Cognitive Neurodynamics*, *17*(1), 221-237.
#'     \doi{10.1007/s11571-022-09813-2}
#' - Durrans, S. R. (1992). Distributions of fractional order statistics in
#'     hydrology. *Water Resources Research*, *28*(6), 1649-1655.
#'     \doi{10.1029/92WR00554}
#'
#' @return `rcogmod_geg()` returns a numeric vector of `n` simulated reaction
#'   times, in seconds. `dcogmod_geg()` returns the density at each element of
#'   `x` - the log density if `log = TRUE` - and `pcogmod_geg()` the cumulative
#'   probability at each element of `q`, honouring `lower.tail` and `log.p`;
#'   both are recycled to the length of the longest argument. `cogmod_geg()`
#'   returns a `brms::custom_family` object, to put on a `brms::bf()` formula.
#'   `cogmod_geg_stanvars()` returns a `brms::stanvars` object holding the
#'   family's Stan `functions` block, to pass to `brms::brm()`, and
#'   `cogmod_geg_lpdf_expose()` compiles that Stan code and returns it as an R
#'   function, for checking the density outside of a model. The remaining
#'   functions are `brms` post-processing methods, called by `brms` rather than
#'   directly: `log_lik_cogmod_geg()` returns a numeric vector holding one
#'   log-likelihood value per posterior draw for observation `i`,
#'   `posterior_predict_cogmod_geg()` a draws x 1 matrix of reaction times
#'   simulated for observation `i`, and `posterior_epred_cogmod_geg()` a
#'   draws x observations matrix of expected reaction times, obtained by
#'   numerical integration.
#'
#' @examples
#' # shape = 1 is the ex-Gaussian, to machine precision
#' x <- seq(0.2, 2, length.out = 5)
#' dcogmod_geg(x, 0.4, 0.1, 0.2, shape = 1)
#' dcogmod_exgaussian(x, 0.4, 0.1, 0.2)
#'
#' rts <- rcogmod_geg(1000, mu = 0.4, sigma = 0.1, tau = 0.2, shape = 2)
#' hist(rts, breaks = 50, xlab = "RT (s)")
#'
#' @export
rcogmod_geg <- function(n, mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1) {
    params <- .prepare_geg(x = NULL, n = n, mu = mu, sigma = sigma,
                           tau = tau, shape = shape)

    # F_GEG(x) = F_EG(x)^shape, so inverting it means inverting the ex-Gaussian
    # at u^(1 / shape). Equivalently u^(1 / shape) ~ Beta(shape, 1).
    u <- stats::runif(params$ndraws)
    .qcogmod_exgaussian(u^(1 / params$shape), params$mu, params$sigma, params$tau)
}


#' @rdname rcogmod_geg
#' @param x Vector of quantiles (observed reaction times).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @export
dcogmod_geg <- function(x, mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1,
                        log = FALSE) {
    params <- .prepare_geg(x = x, n = NULL, mu = mu, sigma = sigma,
                           tau = tau, shape = shape)

    log_density <- base::log(params$shape) +
        (params$shape - 1) *
            .lcdf_exgaussian(params$x, params$mu, params$sigma, params$tau) +
        dcogmod_exgaussian(params$x, params$mu, params$sigma, params$tau,
                           log = TRUE)

    log_density[is.nan(log_density)] <- -Inf

    if (log) log_density else exp(log_density)
}


#' @rdname rcogmod_geg
#' @param q Vector of quantiles.
#' @param lower.tail Logical; if TRUE (default), probabilities are `P(X <= q)`.
#' @param log.p Logical; if TRUE, probabilities are returned on the log scale.
#' @export
pcogmod_geg <- function(q, mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1,
                        lower.tail = TRUE, log.p = FALSE) {
    params <- .prepare_geg(x = q, n = NULL, mu = mu, sigma = sigma,
                           tau = tau, shape = shape)

    # log F_GEG = shape * log F_EG - the whole point of the construction
    lp <- params$shape *
        .lcdf_exgaussian(params$x, params$mu, params$sigma, params$tau)

    if (!lower.tail) {
        # log(1 - exp(lp)), stable for lp near 0 and near -Inf alike
        lp <- .log1mexp(lp)
    }
    if (log.p) lp else exp(lp)
}


# Internals ---------------------------------------------------------------

#' @keywords internal
.prepare_geg <- function(x = NULL, n = NULL, mu, sigma, tau, shape) {
    # `mu` is unchecked for the same reason as in .prepare_exgaussian(): it is a
    # location, and the convolution is defined for any real value.
    if (any(sigma <= 0)) {
        stop("`sigma` must be positive.")
    }
    if (any(tau <= 0)) {
        stop("`tau` must be positive.")
    }
    if (any(shape <= 0)) {
        stop("`shape` must be positive.")
    }

    if (!is.null(x)) {
        m <- length(x)
    } else if (!is.null(n)) {
        if (length(n) != 1 || n < 0 || n != floor(n)) {
            stop("n must be a single non-negative integer.")
        }
        m <- n
    } else {
        stop("Either 'x' or 'n' must be provided.")
    }

    params <- list(
        mu = rep(mu, length.out = m),
        sigma = rep(sigma, length.out = m),
        tau = rep(tau, length.out = m),
        shape = rep(shape, length.out = m)
    )
    if (!is.null(x)) {
        params$x <- rep(x, length.out = m)
    }
    params$ndraws <- m
    params
}


# log(1 - exp(a)) for a <= 0, split at -log(2) the usual way.
#' @keywords internal
.log1mexp <- function(a) {
    out <- numeric(length(a))
    big <- a > -0.693147180559945
    out[big] <- base::log(-expm1(a[big]))
    out[!big] <- log1p(-exp(a[!big]))
    out[a >= 0] <- -Inf
    out
}


# log CDF of the ex-Gaussian.
#
#   F_EG(x) = Phi(z) - exp(sigma^2 / (2 tau^2) - (x - mu) / tau) * Phi(z - sigma/tau)
#
# Both terms are computed on the log scale and subtracted with log1mexp, because
# in the left tail they are individually tiny and very close together - which is
# exactly where the GEG needs the CDF, since `shape < 1` multiplies log F_EG by a
# negative number and any error there is amplified.
#' @keywords internal
.lcdf_exgaussian <- function(x, mu, sigma, tau) {
    z <- (x - mu) / sigma
    la <- stats::pnorm(z, log.p = TRUE)
    lb <- (sigma^2 / (2 * tau^2) - (x - mu) / tau) +
        stats::pnorm(z - sigma / tau, log.p = TRUE)

    d <- lb - la
    # F_EG > 0 everywhere, so lb < la always; clamp only against rounding.
    d[!is.finite(d)] <- -Inf
    d <- pmin(d, -.Machine$double.eps)

    out <- la + .log1mexp(d)
    out[is.nan(out)] <- -Inf
    pmin(out, 0)
}


# Vectorised quantile function of the ex-Gaussian, by bisection on the closed-form
# CDF. There is no closed form, and pulling in gamlss.dist for `qexGAUS` would be
# a dependency for one monotone root-find.
#' @keywords internal
.qcogmod_exgaussian <- function(p, mu, sigma, tau) {
    n <- length(p)
    lo <- mu - 10 * sigma
    hi <- mu + 10 * sigma + tau * 40

    # Widen until the target is bracketed. The ex-Gaussian is unbounded both
    # ways, so neither end can be assumed.
    for (k in 1:60) {
        need <- exp(.lcdf_exgaussian(lo, mu, sigma, tau)) > p
        if (!any(need)) break
        lo[need] <- lo[need] - 10 * sigma[need] - tau[need]
    }
    for (k in 1:60) {
        need <- exp(.lcdf_exgaussian(hi, mu, sigma, tau)) < p
        if (!any(need)) break
        hi[need] <- hi[need] + 10 * sigma[need] + 40 * tau[need]
    }

    for (k in 1:100) {
        mid <- 0.5 * (lo + hi)
        left <- exp(.lcdf_exgaussian(mid, mu, sigma, tau)) < p
        lo[left] <- mid[left]
        hi[!left] <- mid[!left]
    }
    out <- 0.5 * (lo + hi)
    out[p <= 0] <- -Inf
    out[p >= 1] <- Inf
    out
}


# Mean of the GEG, by quadrature. There is no closed form: for integer `shape`
# this is the expected maximum of `shape` ex-Gaussian draws, and the paper says
# outright that the moments have to be found numerically.
#
# Uses E[X] = lo + int_lo^hi (1 - F(x)) dx over a bracket taken from the GEG's
# own quantiles, so the truncation error is bounded by the tail mass left
# outside it (1e-9 either side) times the range. Simpson's rule on 128 intervals
# is then far more accurate than that bound.
#
# Cost is ~130 CDF evaluations per element, so this is the one generic that is
# materially slower here than for a closed-form family. It is fine on the
# prediction grids `modelbased::estimate_means()` builds; on a full data set of
# many thousands of rows, summarise `posterior_predict()` draws instead.
#' @keywords internal
.mean_geg <- function(mu, sigma, tau, shape, nodes = 128L) {
    dims <- dim(mu)
    mu <- as.vector(mu); sigma <- as.vector(sigma)
    tau <- as.vector(tau); shape <- as.vector(shape)

    eps <- 1e-9
    lo <- .qcogmod_exgaussian(eps^(1 / shape), mu, sigma, tau)
    hi <- .qcogmod_exgaussian((1 - eps)^(1 / shape), mu, sigma, tau)

    h <- (hi - lo) / nodes
    acc <- numeric(length(mu))
    for (k in 0:nodes) {
        # Simpson weights: 1, 4, 2, 4, ..., 4, 1
        w <- if (k == 0 || k == nodes) 1 else if (k %% 2 == 1) 4 else 2
        x <- lo + k * h
        surv <- -expm1(shape * .lcdf_exgaussian(x, mu, sigma, tau))
        acc <- acc + w * surv
    }
    out <- lo + acc * h / 3

    if (!is.null(dims)) dim(out) <- dims
    out
}


# Family ------------------------------------------------------------------

#' @rdname rcogmod_geg
#' @param link_mu,link_sigma,link_tau,link_shape Character of the type of link
#'   used to model the GEG parameters. Defaults to `"identity"` for `mu`,
#'   `"softplus"` for `sigma` and `tau`, and `"log"` for `shape`.
#'
#'   `shape` is on a `log` link so that zero on the link scale is `shape = 1`,
#'   the ex-Gaussian. A prior centred at zero is then a prior centred on the
#'   nested model, which is what [cogmod_priors()] supplies.
#' @export
cogmod_geg <- function(
    link_mu = "identity",
    link_sigma = "softplus",
    link_tau = "softplus",
    link_shape = "log"
) {
    brms::custom_family(
        name = "cogmod_geg",
        dpars = c("mu", "sigma", "tau", "shape"),
        lb = c(NA, 0, 0, 0),
        links = c(link_mu, link_sigma, link_tau, link_shape),
        ub = c(NA, NA, NA, NA),
        type = "real"
    )
}


#' @keywords internal
.cogmod_geg_lpdf <- function() {
    "
// Log-likelihood for a single observation from the Generalised Ex-Gaussian.
// Y: observed reaction time.
// mu: location of the Gaussian component. Unbounded, as for the ex-Gaussian.
// sigma: SD of the Gaussian component (> 0).
// tau: mean of the exponential component (> 0).
// shape: power applied to the ex-Gaussian CDF (> 0). shape = 1 is the
//        ex-Gaussian exactly, and the two lines below then reduce to
//        exp_mod_normal_lpdf alone.
real cogmod_geg_lpdf(real Y, real mu, real sigma, real tau, real shape) {
    // Parameter checks
    if (sigma <= 0 || tau <= 0 || shape <= 0) return negative_infinity();

    // Stan's exp_mod_normal is parameterized by the exponential RATE (1 / tau),
    // and ships both the lpdf and the lcdf, so the alpha-power construction
    //     f(x) = shape * F(x)^(shape - 1) * f(x)
    // needs no numerical work of its own.
    real lambda = inv(tau);
    return log(shape)
         + (shape - 1) * exp_mod_normal_lcdf(Y | mu, sigma, lambda)
         + exp_mod_normal_lpdf(Y | mu, sigma, lambda);
}
"
}


#' @rdname rcogmod_geg
#' @export
cogmod_geg_lpdf_expose <- function() {
    insight::check_if_installed("cmdstanr")

    stancode <- paste0(
        "functions {
",
        .cogmod_geg_lpdf(),
        "
}"
    )

    mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
    mod$expose_functions()
    mod$functions$cogmod_geg_lpdf
}


#' @rdname rcogmod_geg
#' @export
cogmod_geg_stanvars <- function() {
    brms::stanvar(scode = .cogmod_geg_lpdf(), block = "functions")
}


# brms methods ------------------------------------------------------------

#' @rdname rcogmod_geg
#' @inheritParams cogmod_lnr
#' @export
log_lik_cogmod_geg <- function(i, prep) {
    if (!"Y" %in% names(prep$data)) {
        stop("Outcome variable 'Y' not found in prep$data.")
    }
    y <- prep$data$Y[i]
    if (is.na(y)) {
        return(NA_real_)
    }

    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)

    n_draws <- length(mu)
    if (n_draws == 0) {
        return(numeric(0))
    }
    y_vec <- rep(y, length.out = n_draws)

    # A draw with an out-of-support parameter is -Inf, not an error - one bad
    # draw must not abort the whole call. The mask is applied *before* the
    # arithmetic rather than after, because `log(shape)` on a negative `shape`
    # would otherwise warn on its way to the NaN we were going to discard.
    bad <- sigma <= 0 | tau <= 0 | shape <= 0
    ll <- rep(-Inf, n_draws)
    if (!all(bad)) {
        k <- !bad
        ll[k] <- base::log(shape[k]) +
            (shape[k] - 1) *
                .lcdf_exgaussian(y_vec[k], mu[k], sigma[k], tau[k]) +
            dcogmod_exgaussian(y_vec[k], mu[k], sigma[k], tau[k], log = TRUE)
    }

    ll[is.nan(ll) | is.na(ll)] <- -Inf
    ll
}


#' @rdname rcogmod_geg
#' @inheritParams cogmod_lnr
#' @export
posterior_predict_cogmod_geg <- function(i, prep, ...) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)

    rcogmod_geg(length(mu), mu = mu, sigma = sigma, tau = tau, shape = shape)
}


#' @rdname rcogmod_geg
#' @export
posterior_epred_cogmod_geg <- function(prep) {
    mu <- brms::get_dpar(prep, "mu")
    sigma <- brms::get_dpar(prep, "sigma")
    tau <- brms::get_dpar(prep, "tau")
    shape <- brms::get_dpar(prep, "shape")

    # No closed form - see .mean_geg(). This is the one generic that is
    # materially slower here than for the ex-Gaussian, whose mean is mu + tau.
    .mean_geg(mu, sigma, tau, shape)
}
