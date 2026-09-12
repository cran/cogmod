# cogmod_geg(): the ex-Gaussian with its CDF raised to `shape`.
#
# The single most important property is the nesting: at shape = 1 this has to be
# cogmod_exgaussian() *exactly*, in R and in Stan alike, or loo_compare() between
# the two is comparing different likelihoods rather than nested models.

test_that("shape = 1 is cogmod_exgaussian() exactly", {
    x <- seq(0.15, 3, length.out = 200)
    expect_equal(
        dcogmod_geg(x, 0.4, 0.1, 0.2, shape = 1, log = TRUE),
        dcogmod_exgaussian(x, 0.4, 0.1, 0.2, log = TRUE)
    )
    # and away from 1 it is genuinely a different density
    expect_false(isTRUE(all.equal(
        dcogmod_geg(x, 0.4, 0.1, 0.2, shape = 2, log = TRUE),
        dcogmod_exgaussian(x, 0.4, 0.1, 0.2, log = TRUE)
    )))
})


test_that("dcogmod_geg integrates to one", {
    for (s in c(0.2, 0.5, 1, 3, 10)) {
        I <- stats::integrate(
            function(z) dcogmod_geg(z, 0.4, 0.1, 0.2, shape = s),
            -Inf, Inf, rel.tol = 1e-10
        )$value
        expect_equal(I, 1, tolerance = 1e-6, label = paste("shape =", s))
    }
})


test_that("pcogmod_geg is the integral of dcogmod_geg", {
    q <- c(0.3, 0.5, 0.8, 1.5)
    for (s in c(0.3, 1, 4)) {
        num <- vapply(q, function(u) {
            stats::integrate(
                function(z) dcogmod_geg(z, 0.4, 0.1, 0.2, shape = s),
                -Inf, u, rel.tol = 1e-10
            )$value
        }, numeric(1))
        expect_equal(pcogmod_geg(q, 0.4, 0.1, 0.2, shape = s), num,
                     tolerance = 1e-7, label = paste("shape =", s))
    }

    # the two flags
    p <- pcogmod_geg(0.7, shape = 2)
    expect_equal(pcogmod_geg(0.7, shape = 2, lower.tail = FALSE), 1 - p)
    expect_equal(pcogmod_geg(0.7, shape = 2, log.p = TRUE), log(p))
})


test_that("rcogmod_geg draws from pcogmod_geg", {
    skip_on_cran()
    set.seed(1)
    for (s in c(0.4, 1, 3)) {
        r <- rcogmod_geg(20000, 0.4, 0.1, 0.2, shape = s)
        u <- pcogmod_geg(r, 0.4, 0.1, 0.2, shape = s)
        expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.001)
    }
})


test_that("the mean has no closed form but the quadrature finds it", {
    for (s in c(0.2, 1, 5)) {
        num <- stats::integrate(
            function(z) z * dcogmod_geg(z, 0.4, 0.1, 0.2, shape = s),
            -Inf, Inf, rel.tol = 1e-10
        )$value
        expect_equal(cogmod:::.mean_geg(0.4, 0.1, 0.2, s), num, tolerance = 1e-5,
                     label = paste("shape =", s))
    }
    # at shape = 1 it must land on the ex-Gaussian's mu + tau
    expect_equal(cogmod:::.mean_geg(0.4, 0.1, 0.2, 1), 0.6, tolerance = 1e-6)

    # and the mean genuinely moves with shape - this is the interpretability
    # cost the family is documented as carrying, so it is pinned rather than
    # left as prose
    expect_lt(cogmod:::.mean_geg(0.4, 0.1, 0.2, 0.2), 0.6)
    expect_gt(cogmod:::.mean_geg(0.4, 0.1, 0.2, 20), 0.6)
})


test_that("shape reaches shapes the ex-Gaussian cannot", {
    # negative skew is the headline claim: the ex-Gaussian is a Normal plus a
    # positive exponential, so its skew is bounded below by zero at every
    # parameter value, while the power transform can pull it negative.
    skew <- function(s, sigma = 0.2, tau = 0.05) {
        f <- function(z) dcogmod_geg(z, 0.4, sigma, tau, shape = s)
        m <- stats::integrate(function(z) z * f(z), -Inf, Inf, rel.tol = 1e-9)$value
        v <- stats::integrate(function(z) (z - m)^2 * f(z), -Inf, Inf, rel.tol = 1e-9)$value
        stats::integrate(function(z) (z - m)^3 * f(z), -Inf, Inf, rel.tol = 1e-9)$value / v^1.5
    }
    expect_lt(skew(0.2), 0)
    expect_gt(skew(1), -1e-6)
})


test_that("parameter checks reject exactly what Stan rejects", {
    expect_error(rcogmod_geg(5, sigma = 0), "sigma")
    expect_error(rcogmod_geg(5, tau = 0), "tau")
    expect_error(rcogmod_geg(5, shape = 0), "shape")
    expect_error(rcogmod_geg(5, shape = -1), "shape")
    # mu is a location and is deliberately unchecked, as for the ex-Gaussian
    expect_silent(rcogmod_geg(5, mu = -1))
})


test_that("cogmod_geg() builds a valid brms custom family", {
    fam <- cogmod_geg()

    expect_s3_class(fam, "customfamily")
    expect_identical(fam$dpars, c("mu", "sigma", "tau", "shape"))
    expect_identical(
        unname(c(fam$link, fam$link_sigma, fam$link_tau, fam$link_shape)),
        c("identity", "softplus", "softplus", "log")
    )
    # `shape` is on `log` so that zero on the link scale is the ex-Gaussian
    expect_equal(fam$lb, list(mu = NA_character_, sigma = "0", tau = "0",
                              shape = "0"))
})


test_that("log_lik_cogmod_geg matches dcogmod_geg and rejects bad parameters", {
    mu <- c(0.4, 0.5)
    prep <- structure(
        list(
            data = list(Y = 0.8),
            dpars = list(
                mu = matrix(mu, nrow = 2, ncol = 1),
                sigma = matrix(c(0.1, 0.08), nrow = 2, ncol = 1),
                tau = matrix(c(0.2, 0.25), nrow = 2, ncol = 1),
                shape = matrix(c(1, 2.5), nrow = 2, ncol = 1)
            ),
            ndraws = 2
        ),
        class = "brmsprep"
    )
    expect_equal(
        log_lik_cogmod_geg(1, prep),
        dcogmod_geg(rep(0.8, 2), mu, c(0.1, 0.08), c(0.2, 0.25), c(1, 2.5),
                    log = TRUE)
    )

    prep$dpars$shape[1, 1] <- -1
    expect_equal(log_lik_cogmod_geg(1, prep)[1], -Inf)
})


test_that("cogmod_priors and cogmod_inits know the family", {
    set.seed(3)
    d <- data.frame(
        RT = rcogmod_geg(200, mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1.5),
        Condition = factor(rep(c("a", "b"), length.out = 200))
    )
    f <- brms::bf(RT ~ Condition, sigma ~ Condition, tau ~ Condition,
                  shape ~ Condition, family = cogmod_geg())

    p <- expect_silent(cogmod_priors(f, d))
    pick <- function(cls, dpar, coef = "") {
        p$prior[p$class == cls & p$dpar == dpar & p$coef == coef]
    }
    # centred on the ex-Gaussian: log(1) = 0. The prior exists because `shape`
    # is ~0.98 correlated with `mu` and will otherwise wander up that ridge.
    expect_equal(pick("Intercept", "shape"), "normal(0, 0.5)")
    expect_equal(pick("Intercept", ""), "normal(0.4, 0.25)")

    vals <- cogmod_inits(f, d, jitter = 0)(1)
    # log link, so the ex-Gaussian start of shape = 1 is zero on the link scale
    expect_equal(vals$Intercept_shape, 0)

    # omitted from bf(), shape lands on the natural scale instead
    g <- brms::bf(RT ~ 1, family = cogmod_geg())
    expect_equal(cogmod_priors(g, d)$prior[cogmod_priors(g, d)$class == "shape"],
                 "lognormal(0, 0.5)")
    expect_equal(cogmod_inits(g, d, jitter = 0)(1)$shape, 1)
})


test_that("cogmod_stanvars() finds the family", {
    f <- brms::bf(RT ~ 1, family = cogmod_geg())
    sv <- expect_silent(cogmod_stanvars(f))
    expect_true(grepl("cogmod_geg_lpdf", sv[[1]]$scode, fixed = TRUE))
})


test_that("the Stan lpdf agrees with the R density", {
    skip_on_cran()
    skip_if_not_installed("cmdstanr")

    lpdf <- stan_fun("cogmod_geg")

    g <- covering_grid(
        y = c(0.25, 0.5, 0.9, 1.8),
        mu = c(-0.1, 0.4, 0.7),
        sigma = c(0.03, 0.1, 0.3),
        tau = c(0.05, 0.2, 0.6),
        shape = c(0.25, 1, 2.5, 8)
    )
    for (i in seq_len(nrow(g))) {
        stan <- lpdf(g$y[i], g$mu[i], g$sigma[i], g$tau[i], g$shape[i])
        r <- dcogmod_geg(g$y[i], g$mu[i], g$sigma[i], g$tau[i], g$shape[i],
                         log = TRUE)
        expect_lt(abs(stan - r) / max(1, abs(r)), 1e-8)
    }

    # shape = 1 has to reduce to the ex-Gaussian on the Stan side too, or the
    # nesting the family is sold on holds only in R
    eg <- stan_fun("cogmod_exgaussian")
    expect_equal(lpdf(0.8, 0.4, 0.1, 0.2, 1), eg(0.8, 0.4, 0.1, 0.2),
                 tolerance = 1e-12)

    # invalid arguments are rejected on both sides
    expect_equal(lpdf(0.5, 0.4, 0, 0.2, 1), -Inf)
    expect_equal(lpdf(0.5, 0.4, 0.1, 0, 1), -Inf)
    expect_equal(lpdf(0.5, 0.4, 0.1, 0.2, 0), -Inf)
})
