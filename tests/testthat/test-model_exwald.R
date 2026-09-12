context("ex-Wald")

# Brute-force convolution of the Wald with the exponential, independent of
# everything the package does. Bounded parameter ranges only: at an extreme
# threshold the integrand is near-singular at zero and integrate() gives up, so
# the normalisation test below is the one that covers the whole grid.
ref_conv <- function(t, drift, boundary, tau) {
  g <- 1 / tau
  dw <- function(u) {
    boundary / sqrt(2 * pi * u^3) * exp(-(boundary - drift * u)^2 / (2 * u))
  }
  vapply(t, function(x) {
    tryCatch(
      stats::integrate(function(u) dw(u) * g * exp(-g * (x - u)), 0, x,
                       rel.tol = 1e-12, subdivisions = 2000)$value,
      error = function(e) NA_real_
    )
  }, numeric(1))
}

test_that("both branches match a brute-force convolution", {
  tt <- c(0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3)
  grid <- covering_grid(
    drift = c(2, 3, 5, 8),
    boundary = c(0.3, 0.5, 1, 2),
    tau = c(0.05, 0.1, 0.2, 0.3),
    # tau > 2 / drift^2 is the closed form; below it the Faddeeva branch. Both
    # must appear, so one of each is pinned in.
    always = function(g) {
      (g$drift == 3 & g$boundary == 0.5 & g$tau == 0.3) |
        (g$drift == 3 & g$boundary == 0.5 & g$tau == 0.1)
    }
  )
  seen <- c(closed = FALSE, faddeeva = FALSE)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    want <- ref_conv(tt, g$drift, g$boundary, g$tau)
    keep <- !is.na(want) & want > 1e-8
    if (!any(keep)) next
    seen[[if (g$drift^2 >= 2 / g$tau) "closed" else "faddeeva"]] <- TRUE
    got <- dcogmod_exwald(tt, mu = g$drift, boundary = g$boundary, tau = g$tau,
                          ndt = 0)
    expect_equal(got[keep], want[keep], tolerance = 1e-4,
                 label = sprintf("mu=%s boundary=%s tau=%s",
                                 g$drift, g$boundary, g$tau))
  }
  expect_true(all(seen))
})


# The sharper test, and the one that spans the whole parameter region: a density
# that is wrong anywhere does not integrate to 1, and the ex-Wald mean is known
# exactly. Neither depends on a numerical reference.
test_that("the density normalises and has the right mean", {
  grid <- covering_grid(
    drift = c(1.5, 2, 3, 5, 8, 12),
    boundary = c(0.2, 0.5, 1, 3),
    tau = c(0.03, 0.1, 0.2, 0.4),
    ndt = c(0, 0.2),
    always = function(g) g$drift == 1.5 & g$tau == 0.03
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    f <- function(x) {
      dcogmod_exwald(x, mu = g$drift, boundary = g$boundary, tau = g$tau,
                     ndt = g$ndt)
    }
    lab <- sprintf("mu=%s boundary=%s tau=%s ndt=%s",
                   g$drift, g$boundary, g$tau, g$ndt)
    tot <- stats::integrate(f, g$ndt, Inf, rel.tol = 1e-9,
                            subdivisions = 2000)$value
    expect_equal(tot, 1, tolerance = 1e-6, label = lab)
    m <- stats::integrate(function(x) x * f(x), g$ndt, Inf, rel.tol = 1e-9,
                          subdivisions = 2000)$value
    expect_equal(m, g$ndt + g$boundary / g$drift + g$tau, tolerance = 1e-5,
                 label = lab)
  }
})


test_that("the two branches meet at mu^2 = 2 / tau", {
  tt <- c(0.05, 0.2, 0.5, 1.2, 3)
  for (p in list(c(3, 0.5), c(5, 2), c(2, 0.3), c(8, 3))) {
    drift <- p[1]; boundary <- p[2]
    seam <- 2 / drift^2
    hi <- dcogmod_exwald(tt, mu = drift, boundary = boundary,
                         tau = seam * (1 + 1e-9), ndt = 0, log = TRUE)
    lo <- dcogmod_exwald(tt, mu = drift, boundary = boundary,
                         tau = seam * (1 - 1e-9), ndt = 0, log = TRUE)
    expect_equal(lo, hi, tolerance = 1e-6,
                 label = sprintf("mu=%s boundary=%s", drift, boundary))
  }
})


test_that("Re w matches the Faddeeva function where it is known in closed form", {
  # w(0) = 1, and w(i y) = exp(y^2) erfc(y) for real y > 0.
  expect_equal(cogmod:::.re_faddeeva(0, 0), 1, tolerance = 1e-9)
  y <- c(0.05, 0.2, 1, 2, 4)
  expect_equal(cogmod:::.re_faddeeva(0, y),
               exp(y^2) * 2 * stats::pnorm(-y * sqrt(2)), tolerance = 1e-8)
  # On the real axis Re w(x) = exp(-x^2).
  x <- c(0.05, 0.5, 1.5, 3)
  expect_equal(cogmod:::.re_faddeeva(x, 0), exp(-x^2), tolerance = 1e-8)
  # Positive throughout the upper half plane - it is the Voigt function.
  gr <- expand.grid(x = c(0, 0.1, 1, 5, 20), y = c(1e-4, 0.01, 1, 5, 20))
  expect_true(all(cogmod:::.re_faddeeva(gr$x, gr$y) > 0))
})


test_that("the mixture and the outlier component behave", {
  x <- c(0.25, 0.4, 0.8, 2)
  got <- dcogmod_exwald(x, ndt = 0.2, poutlier = 0.03, log = TRUE)
  dec <- cogmod:::.dexwald_raw(x - 0.2, 3, 0.5, 0.15)
  out <- log(2) + stats::dnorm(x, 0, 0.2, log = TRUE)
  expect_equal(got, log(0.03 * exp(out) + 0.97 * exp(dec)), tolerance = 1e-12)

  expect_equal(dcogmod_exwald(0.1, ndt = 0.2, poutlier = 0.03, log = TRUE),
               log(0.03) + log(2) + stats::dnorm(0.1, 0, 0.2, log = TRUE),
               tolerance = 1e-12)
  expect_equal(dcogmod_exwald(0.1, ndt = 0.2, poutlier = 0), 0)
})


test_that("rcogmod_exwald matches its density", {
  set.seed(9)
  x <- rcogmod_exwald(3e5, mu = 3, boundary = 0.5, tau = 0.15, ndt = 0.2)
  expect_equal(mean(x), 0.2 + 0.5 / 3 + 0.15, tolerance = 0.01)
  cdf <- function(q) {
    stats::integrate(function(u) exp(cogmod:::.dexwald_raw(u, 3, 0.5, 0.15)),
                     0, q - 0.2, rel.tol = 1e-9)$value
  }
  for (p in c(0.1, 0.5, 0.9)) {
    expect_equal(cdf(unname(stats::quantile(x, p))), p, tolerance = 0.01)
  }
})


test_that("parameter bounds are enforced", {
  expect_warning(expect_equal(dcogmod_exwald(0.5, mu = 0), 0), "mu")
  expect_warning(expect_equal(dcogmod_exwald(0.5, boundary = 0), 0), "boundary")
  expect_warning(expect_equal(dcogmod_exwald(0.5, tau = 0), 0), "tau")
  expect_error(rcogmod_exwald(5, mu = 0), "greater than 0")
  expect_error(rcogmod_exwald(5, tau = 0), "greater than 0")
  expect_error(rcogmod_exwald(5, ndt = -0.1), "ndt")
})


test_that("the family, its priors and its epred are wired up", {
  fam <- cogmod_exwald()
  expect_equal(fam$dpars, c("mu", "boundary", "tau", "ndt", "poutlier"))
  expect_equal(fam$name, "cogmod_exwald")

  set.seed(5)
  d <- data.frame(
    RT = rcogmod_exwald(200, ndt = 0.2, poutlier = 0.02),
    Condition = factor(rep(c("a", "b"), 100))
  )
  f <- brms::bf(RT ~ Condition, tau ~ Condition, ndt ~ 1, poutlier ~ 1,
                family = cogmod_exwald())
  # brms leaves `tau` flat, so this is a fill rather than an override.
  raw <- brms::get_prior(f, d, family = f$family)
  expect_true(all(!nzchar(raw$prior[raw$dpar == "tau"])))

  p <- expect_silent(cogmod_priors(f, d))
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "tau"],
               "normal(-1.5, 0.7)")
  expect_equal(p$prior[p$class == "b" & p$dpar == "tau" & p$coef == ""],
               "normal(0, 0.5)")

  expect_silent(
    brms::make_stancode(f, data = d, family = f$family, prior = p,
                        stanvars = cogmod_stanvars(f))
  )
})


context("ex-Wald - brms")

test_that("Stan cogmod_exwald_lpdf matches R dcogmod_exwald", {
  fun <- stan_fun("cogmod_exwald")
  grid <- covering_grid(
    drift = c(1.5, 3, 5, 8),
    boundary = c(0.2, 0.5, 1, 3),
    tau = c(0.03, 0.1, 0.2, 0.4),
    ndt = c(0, 0.2),
    poutlier = c(0, 0.02, 0.3),
    always = function(g) {
      (g$drift == 3 & g$tau == 0.3 & g$ndt == 0.2) |
        (g$drift == 3 & g$tau == 0.03 & g$ndt == 0.2)
    }
  )
  x <- c(0.05, 0.19, 0.21, 0.5, 1.5, 6)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    r <- dcogmod_exwald(x, mu = g$drift, boundary = g$boundary, tau = g$tau,
                        ndt = g$ndt, poutlier = g$poutlier, log = TRUE)
    s <- vapply(x, function(xx) {
      fun(xx, g$drift, g$boundary, g$tau, g$ndt, g$poutlier)
    }, numeric(1))
    expect_equal(s, r, tolerance = 1e-8,
                 label = sprintf("mu=%s boundary=%s tau=%s ndt=%s p=%s",
                                 g$drift, g$boundary, g$tau, g$ndt, g$poutlier))
  }
  expect_equal(fun(0.5, 0, 0.5, 0.15, 0.2, 0.02), -Inf)
  expect_equal(fun(0.5, 3, 0.5, 0, 0.2, 0.02), -Inf)
  expect_equal(fun(-1, 3, 0.5, 0.15, 0.2, 0.02), -Inf)
})


test_that("ex-Wald recovers its parameters with brms", {
  skip_if_not_slow()
  skip_if_not_installed("cmdstanr")

  set.seed(31)
  d <- data.frame(
    RT = rcogmod_exwald(2000, mu = 3, boundary = 0.5, tau = 0.15, ndt = 0.25,
                        poutlier = 0.02)
  )
  f <- brms::bf(RT ~ 1, boundary ~ 1, tau ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_exwald())
  m <- brms::brm(
    f, data = d, prior = cogmod_priors(f, d), init = cogmod_inits(f, d),
    stanvars = cogmod_stanvars(f), chains = 2, iter = 1000,
    backend = "cmdstanr", refresh = 0
  )
  sp <- function(z) log1p(exp(z))
  fx <- brms::fixef(m)
  expect_equal(unname(sp(fx["Intercept", "Estimate"])), 3, tolerance = 1.2)
  expect_equal(unname(sp(fx["tau_Intercept", "Estimate"])), 0.15,
               tolerance = 0.1)
  expect_equal(unname(exp(fx["ndt_Intercept", "Estimate"])), 0.25,
               tolerance = 0.1)
  expect_true(all(brms::rhat(m)[!is.na(brms::rhat(m))] < 1.05))
})
