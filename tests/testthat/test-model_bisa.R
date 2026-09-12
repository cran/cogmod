context("Birnbaum-Saunders")

# Everything about this family is closed form, so nothing here needs a numerical
# reference. Three independent statements pin the density: it integrates to one
# with the stated moments, it equals an equal mixture of the Wald and that
# Wald's length-biased version, and the transform (mu t - boundary) / sqrt(t) is
# exactly standard normal - which is the sharpest of the three, and the one that
# checks the sampler rather than the density.

# The textbook (a, b) parameterization, written from scratch, so a mistake in
# the mechanistic mapping cannot hide.
ref_ab <- function(t, a, b) {
  z <- (sqrt(t / b) - sqrt(b / t)) / a
  stats::dnorm(z) * (sqrt(t / b) + sqrt(b / t)) / (2 * a * t)
}

test_that("the mechanistic parameters match the textbook (a, b) form", {
  grid <- covering_grid(
    drift = c(0.5, 1.5, 3, 5, 12),
    boundary = c(0.05, 0.2, 0.5, 1, 3)
  )
  tt <- c(0.02, 0.05, 0.1, 0.25, 0.6, 1.5, 4)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    # b is the scale, a the shape. See ?rcogmod_bisa.
    want <- ref_ab(tt, a = 1 / sqrt(g$drift * g$boundary), b = g$boundary / g$drift)
    got <- dcogmod_bisa(tt, mu = g$drift, boundary = g$boundary, ndt = 0)
    expect_equal(got, want, tolerance = 1e-10,
                 label = sprintf("mu=%s boundary=%s", g$drift, g$boundary))
  }
})


test_that("the density is an equal mixture of a Wald and its length-biased twin", {
  grid <- covering_grid(
    drift = c(0.5, 2, 3, 8),
    boundary = c(0.1, 0.5, 1, 2.5)
  )
  tt <- c(0.03, 0.1, 0.3, 0.9, 2.5)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    wald <- exp(cogmod:::.dwald_raw(tt, g$drift, g$boundary))
    # The length-biased Wald is the Wald reweighted by t / E[T], and E[T] here
    # is the Wald's own mean boundary / drift.
    biased <- wald * tt / (g$boundary / g$drift)
    got <- dcogmod_bisa(tt, mu = g$drift, boundary = g$boundary, ndt = 0)
    expect_equal(got, 0.5 * wald + 0.5 * biased, tolerance = 1e-10,
                 label = sprintf("mu=%s boundary=%s", g$drift, g$boundary))
  }
})


test_that("the density normalises and has the stated mean and variance", {
  grid <- covering_grid(
    drift = c(0.8, 1.5, 3, 6, 12),
    boundary = c(0.1, 0.3, 0.5, 1, 3),
    ndt = c(0, 0.2),
    always = function(g) g$drift == 0.8 & g$boundary == 3
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    f <- function(x) dcogmod_bisa(x, mu = g$drift, boundary = g$boundary,
                                  ndt = g$ndt)
    lab <- sprintf("mu=%s boundary=%s ndt=%s", g$drift, g$boundary, g$ndt)
    tot <- stats::integrate(f, g$ndt, Inf, rel.tol = 1e-10,
                            subdivisions = 2000)$value
    expect_equal(tot, 1, tolerance = 1e-7, label = lab)

    m_want <- g$ndt + g$boundary / g$drift + 1 / (2 * g$drift^2)
    m <- stats::integrate(function(x) x * f(x), g$ndt, Inf, rel.tol = 1e-10,
                          subdivisions = 2000)$value
    expect_equal(m, m_want, tolerance = 1e-6, label = lab)

    v_want <- g$boundary / g$drift^3 + 5 / (4 * g$drift^4)
    v <- stats::integrate(function(x) (x - m_want)^2 * f(x), g$ndt, Inf,
                          rel.tol = 1e-10, subdivisions = 2000)$value
    expect_equal(v, v_want, tolerance = 1e-5, label = lab)
  }
})


test_that("the median is exactly ndt + boundary / mu", {
  # The CDF is Phi((mu t - boundary) / sqrt(t)), which is 1/2 at t = boundary /
  # mu whatever the parameters - so the median needs no root finding.
  for (p in list(c(0.5, 2), c(3, 0.5), c(9, 0.15))) {
    drift <- p[1]; boundary <- p[2]
    half <- stats::integrate(
      function(x) dcogmod_bisa(x, mu = drift, boundary = boundary, ndt = 0.2),
      0.2, 0.2 + boundary / drift, rel.tol = 1e-10
    )$value
    expect_equal(half, 0.5, tolerance = 1e-7,
                 label = sprintf("mu=%s boundary=%s", drift, boundary))
  }
})


test_that("the mixture and the outlier component behave", {
  x <- c(0.25, 0.4, 0.8, 2)
  got <- dcogmod_bisa(x, ndt = 0.2, poutlier = 0.03, log = TRUE)
  dec <- cogmod:::.dbisa_raw(x - 0.2, 3, 0.5)
  out <- log(2) + stats::dnorm(x, 0, 0.2, log = TRUE)
  expect_equal(got, log(0.03 * exp(out) + 0.97 * exp(dec)), tolerance = 1e-12)

  expect_equal(dcogmod_bisa(0.1, ndt = 0.2, poutlier = 0.03, log = TRUE),
               log(0.03) + log(2) + stats::dnorm(0.1, 0, 0.2, log = TRUE),
               tolerance = 1e-12)
  expect_equal(dcogmod_bisa(0.1, ndt = 0.2, poutlier = 0), 0)
})


test_that("rcogmod_bisa draws exactly from the distribution", {
  # (mu t - boundary) / sqrt(t) is standard normal by construction, so this
  # tests the sampler against the definition rather than against the density,
  # and it is exact - no binning, no tolerance on a moment.
  set.seed(11)
  grid <- covering_grid(
    drift = c(0.7, 2, 3, 9),
    boundary = c(0.1, 0.5, 1, 2)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    x <- rcogmod_bisa(20000, mu = g$drift, boundary = g$boundary, ndt = 0)
    z <- (g$drift * x - g$boundary) / sqrt(x)
    expect_gt(stats::ks.test(z, "pnorm")$p.value, 1e-4)
  }

  x <- rcogmod_bisa(3e5, mu = 3, boundary = 0.5, ndt = 0.2)
  expect_equal(mean(x), 0.2 + 0.5 / 3 + 1 / 18, tolerance = 0.005)
  expect_equal(stats::median(x), 0.2 + 0.5 / 3, tolerance = 0.005)
  expect_equal(stats::var(x), 0.5 / 27 + 5 / 324, tolerance = 0.005)
})


test_that("parameter bounds are enforced", {
  expect_warning(expect_equal(dcogmod_bisa(0.5, mu = 0), 0), "mu")
  expect_warning(expect_equal(dcogmod_bisa(0.5, boundary = 0), 0), "boundary")
  expect_error(rcogmod_bisa(5, mu = 0), "greater than 0")
  expect_error(rcogmod_bisa(5, boundary = 0), "greater than 0")
  expect_error(rcogmod_bisa(5, ndt = -0.1), "ndt")
})


test_that("the family, its priors and its epred are wired up", {
  fam <- cogmod_bisa()
  expect_equal(fam$dpars, c("mu", "boundary", "ndt", "poutlier"))
  expect_equal(fam$name, "cogmod_bisa")

  set.seed(5)
  d <- data.frame(
    RT = rcogmod_bisa(200, ndt = 0.2, poutlier = 0.02),
    Condition = factor(rep(c("a", "b"), 100))
  )
  f <- brms::bf(RT ~ Condition, boundary ~ Condition, ndt ~ 1, poutlier ~ 1,
                family = cogmod_bisa())
  # Two identified parameters and no flat direction, so unlike
  # cogmod_invgaussian()'s sigmadrift neither decision dpar needs a prior of its
  # own; cogmod_priors() fills ndt and poutlier and leaves the rest blanket.
  p <- expect_silent(cogmod_priors(f, d))
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "ndt"],
               "normal(-1.2, 0.2)")
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "poutlier"],
               "normal(-5, 1)")

  expect_silent(
    brms::make_stancode(f, data = d, family = f$family, prior = p,
                        stanvars = cogmod_stanvars(f))
  )
})


context("Birnbaum-Saunders - brms")

test_that("Stan cogmod_bisa_lpdf matches R dcogmod_bisa", {
  fun <- stan_fun("cogmod_bisa")
  grid <- covering_grid(
    drift = c(0.8, 1.5, 3, 5, 12),
    boundary = c(0.1, 0.3, 0.5, 1, 3),
    ndt = c(0, 0.2),
    poutlier = c(0, 0.02, 0.3)
  )
  x <- c(0.05, 0.19, 0.21, 0.5, 1.5, 6)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    r <- dcogmod_bisa(x, mu = g$drift, boundary = g$boundary, ndt = g$ndt,
                      poutlier = g$poutlier, log = TRUE)
    s <- vapply(x, function(xx) {
      fun(xx, g$drift, g$boundary, g$ndt, g$poutlier)
    }, numeric(1))
    expect_equal(s, r, tolerance = 1e-8,
                 label = sprintf("mu=%s boundary=%s ndt=%s p=%s",
                                 g$drift, g$boundary, g$ndt, g$poutlier))
  }
  expect_equal(fun(0.5, 0, 0.5, 0.2, 0.02), -Inf)
  expect_equal(fun(0.5, 3, 0, 0.2, 0.02), -Inf)
  expect_equal(fun(-1, 3, 0.5, 0.2, 0.02), -Inf)
})


test_that("Birnbaum-Saunders recovers its parameters with brms", {
  skip_if_not_slow()
  skip_if_not_installed("cmdstanr")

  set.seed(31)
  d <- data.frame(
    RT = rcogmod_bisa(2000, mu = 3, boundary = 0.5, ndt = 0.25,
                      poutlier = 0.02)
  )
  f <- brms::bf(RT ~ 1, boundary ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_bisa())
  m <- brms::brm(
    f, data = d, prior = cogmod_priors(f, d), init = cogmod_inits(f, d),
    stanvars = cogmod_stanvars(f), chains = 2, iter = 1000,
    backend = "cmdstanr", refresh = 0
  )
  sp <- function(z) log1p(exp(z))
  fx <- brms::fixef(m)
  expect_equal(unname(sp(fx["Intercept", "Estimate"])), 3, tolerance = 0.8)
  expect_equal(unname(sp(fx["boundary_Intercept", "Estimate"])), 0.5,
               tolerance = 0.2)
  expect_equal(unname(exp(fx["ndt_Intercept", "Estimate"])), 0.25,
               tolerance = 0.1)
  expect_true(all(brms::rhat(m)[!is.na(brms::rhat(m))] < 1.05))
})
