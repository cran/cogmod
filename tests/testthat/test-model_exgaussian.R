context("Ex-Gaussian (Classical Parameterization) - brms")

test_that("log_lik_cogmod_exgaussian matches brms::dexgaussian under reparameterization", {
    # brms::dexgaussian() uses mu = mean of the *whole* distribution and
    # beta = mean of the exponential component. Our classical mu (Gaussian
    # location) relates to brms's mu via: mu_brms = mu_classical + tau.
    mu <- 0.5
    sigma <- 0.12
    tau <- 0.25
    y_vals <- c(0.1, 0.3, 0.6, 0.9, 1.5, 2.5)

    n_draws <- 10
    prep <- structure(
        list(
            data = list(Y = y_vals),
            dpars = list(
                mu = rep(mu, n_draws),
                sigma = rep(sigma, n_draws),
                tau = rep(tau, n_draws)
            )
        ),
        class = "brmsprep"
    )

    for (i in seq_along(y_vals)) {
        ll <- log_lik_cogmod_exgaussian(i, prep)
        expected <- brms::dexgaussian(
            y_vals[i],
            mu = mu + tau,
            sigma = sigma,
            beta = tau,
            log = TRUE
        )

        expect_equal(
            ll,
            rep(expected, n_draws),
            tolerance = 1e-6,
            label = sprintf("log_lik at y = %.2f", y_vals[i])
        )
    }
})

test_that("log_lik_cogmod_exgaussian returns -Inf for invalid parameters", {
    prep <- structure(
        list(
            data = list(Y = c(0.5)),
            dpars = list(
                mu = c(0.5, 0.5),
                sigma = c(0.1, -0.1),
                tau = c(-0.2, 0.2)
            )
        ),
        class = "brmsprep"
    )

    ll <- suppressWarnings(log_lik_cogmod_exgaussian(1, prep))
    expect_true(all(ll == -Inf))
})

test_that("posterior_predict_cogmod_exgaussian recovers theoretical mean and SD", {
    set.seed(123)
    n_draws <- 50000
    mu <- 0.4
    sigma <- 0.08
    tau <- 0.2

    prep <- structure(
        list(
            dpars = list(
                mu = rep(mu, n_draws),
                sigma = rep(sigma, n_draws),
                tau = rep(tau, n_draws)
            )
        ),
        class = "brmsprep"
    )

    rts <- posterior_predict_cogmod_exgaussian(1, prep)

    theo_mean <- mu + tau
    theo_sd <- sqrt(sigma^2 + tau^2)

    expect_equal(
        mean(rts),
        theo_mean,
        tolerance = 0.02,
        label = "Mean recovery"
    )
    expect_equal(sd(rts), theo_sd, tolerance = 0.02, label = "SD recovery")
    expect_true(all(is.finite(rts)), label = "All simulated RTs are finite")
})

test_that("posterior_epred_cogmod_exgaussian equals mu + tau", {
    mu <- matrix(c(0.3, 0.5, 0.7), nrow = 3, ncol = 2)
    tau <- matrix(c(0.1, 0.2, 0.3), nrow = 3, ncol = 2)
    prep <- structure(
        list(dpars = list(mu = mu, tau = tau)),
        class = "brmsprep"
    )

    epred <- posterior_epred_cogmod_exgaussian(prep)
    expect_equal(epred, mu + tau)
})

test_that("cogmod_exgaussian() builds a valid brms custom family", {
    fam <- cogmod_exgaussian()

    expect_s3_class(fam, "customfamily")
    expect_identical(fam$dpars, c("mu", "sigma", "tau"))
    expect_identical(
        unname(c(fam$link, fam$link_sigma, fam$link_tau)),
        c("identity", "softplus", "softplus")
    )
    # `mu` is a location, not a scale, so it is unbounded - the convolution is
    # defined for any real value and the Stan lpdf only ever checked the two SDs
    expect_equal(fam$lb, list(mu = NA_character_, sigma = "0", tau = "0"))
})
