context("Log-Student-t")

# The R density, written out independently of the registry entry, so that a
# mistake in one is not copied into the other.
ref_ldec <- function(t, mu, sigma, dof) {
  stats::dt((log(t) - mu) / sigma, df = dof, log = TRUE) - log(sigma) - log(t)
}

test_that("the decision component is a normalised shifted scaled t", {
  # Integrated in u = (log(RT - ndt) - mu) / sigma, where the decision component
  # is a plain Student-t. Two things force that substitution and the finite
  # range that goes with it. On the raw scale the density is unbounded at `ndt`
  # (see below), which adaptive quadrature cannot resolve at all. And the tail
  # is heavy enough that integrate() cannot be trusted to reach it: over
  # (-Inf, Inf) at dof = 1 it returns 0.993 and reports no error, while a finite
  # range on the log scale truncates the Cauchy outright.
  #
  # So the range is finite and the target is the EXACT mass inside it, which
  # makes the check sharp at every dof rather than only at the light-tailed end.
  # `poutlier` is left out because the outlier component's own truncated mass
  # would have to be added back; the mixture weighting is checked exactly in the
  # next block instead.
  U <- 30
  grid <- covering_grid(
    mu = c(-1.2, -0.7, 0),
    sigma = c(0.2, 0.4, 0.8),
    dof = c(1, 2, 5, 30),
    ndt = c(0, 0.2),
    always = function(g) g$dof == 1
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    tot <- stats::integrate(
      function(u) {
        lt <- g$mu + g$sigma * u
        # The Jacobian is folded in on the LOG scale: far out, the density has
        # underflowed to zero while exp(lt) has overflowed to Inf, and their
        # product on the natural scale is NaN rather than 0.
        exp(dcogmod_logstudent(g$ndt + exp(lt), mu = g$mu, sigma = g$sigma,
                               dof = g$dof, ndt = g$ndt, poutlier = 0,
                               log = TRUE) + log(g$sigma) + lt)
      },
      -U, U, rel.tol = 1e-10, subdivisions = 2000
    )$value
    expected <- stats::pt(U, g$dof) - stats::pt(-U, g$dof)
    expect_equal(tot, expected, tolerance = 1e-6,
                 label = sprintf("mu=%s sigma=%s dof=%s ndt=%s",
                                 g$mu, g$sigma, g$dof, g$ndt))
  }
})


test_that("dcogmod_logstudent is the shifted scaled t, mixed with the outlier component", {
  x <- c(0.25, 0.4, 0.8, 2, 6)
  got <- dcogmod_logstudent(x, mu = -0.7, sigma = 0.4, dof = 5, ndt = 0.2,
                            poutlier = 0.03, log = TRUE)
  dec <- ref_ldec(x - 0.2, -0.7, 0.4, 5)
  out <- log(2) + stats::dnorm(x, 0, 0.2, log = TRUE)
  expect_equal(got, log(0.03 * exp(out) + 0.97 * exp(dec)), tolerance = 1e-12)

  # Below ndt only the outlier component has any density, and without it the
  # response is impossible.
  expect_equal(dcogmod_logstudent(0.1, ndt = 0.2, poutlier = 0.03, log = TRUE),
               log(0.03) + log(2) + stats::dnorm(0.1, 0, 0.2, log = TRUE),
               tolerance = 1e-12)
  expect_equal(dcogmod_logstudent(0.1, ndt = 0.2, poutlier = 0), 0)
})


test_that("rcogmod_logstudent matches the density it is drawn from", {
  set.seed(11)
  # Quantiles rather than moments: this family has none. log(RT - ndt) is a
  # scaled t exactly, which is what makes the check sharp.
  x <- rcogmod_logstudent(4e5, mu = -0.7, sigma = 0.4, dof = 5, ndt = 0.2)
  z <- (log(x - 0.2) + 0.7) / 0.4
  p <- c(0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98)
  expect_equal(unname(stats::quantile(z, p)), stats::qt(p, df = 5),
               tolerance = 0.02)

  # With an outlier component, the fast spike is there and matches its weight.
  set.seed(12)
  y <- rcogmod_logstudent(2e5, dof = 5, ndt = 0.3, poutlier = 0.1)
  expect_equal(mean(y < 0.3), 0.1 * (2 * stats::pnorm(0.3 / 0.2) - 1),
               tolerance = 0.01)
})


test_that("dof -> Inf is the LogNormal", {
  x <- c(0.25, 0.5, 1, 3)
  ln <- dcogmod_lognormal(x, mu = -0.7, sigma = 0.4, ndt = 0.2, log = TRUE)
  # The limit is approached slowly - that is the flat direction cogmod_priors()
  # fences - so the tolerance has to loosen as dof falls.
  for (d in c(1e6, 1e4, 1e2)) {
    got <- dcogmod_logstudent(x, mu = -0.7, sigma = 0.4, dof = d, ndt = 0.2,
                              log = TRUE)
    expect_equal(got, ln, tolerance = 50 / d,
                 label = paste("dof =", d))
  }
  # ...and at a dof anyone would actually fit, it is a different distribution.
  # The claim the family is for is about the slow TAIL, which is far bigger than
  # the density ratio at any one point: P(decision time > 5 s) is five orders of
  # magnitude larger than the LogNormal's, where the densities at 5 s differ by
  # a factor of 85.
  tail_t <- stats::pt((log(5 - 0.2) + 0.7) / 0.4, df = 5, lower.tail = FALSE)
  tail_ln <- stats::plnorm(5 - 0.2, -0.7, 0.4, lower.tail = FALSE)
  expect_gt(tail_t / tail_ln, 1e5)
  expect_gt(dcogmod_logstudent(5, dof = 5, ndt = 0.2) /
              dcogmod_lognormal(5, ndt = 0.2), 50)
})


test_that("the density is unbounded at ndt but the spike is integrable", {
  # Taken at ndt = 0, because `ndt + t` is not representable for the small `t`
  # this needs: 0.2 + 1e-30 is exactly 0.2 in double precision. The shift is a
  # translation, so the behaviour at the shift is the behaviour at zero.
  t <- 10^-c(2, 8, 30, 60)
  d <- dcogmod_logstudent(t, dof = 2, ndt = 0)
  # `t` decreases along the vector, and the density grows as it does.
  expect_true(all(diff(d) > 0))
  expect_gt(d[length(d)], 1e40)
  # A LogNormal goes the other way, to zero.
  expect_lt(dcogmod_lognormal(t[2], ndt = 0), 1e-80)

  # Integrable all the same: the mass just above ndt stays finite and small.
  # P(decision time < eps) = pt(log(eps) ... ), taken on the log scale.
  p_spike <- stats::pt((log(1e-3) + 0.7) / 0.4, df = 2)
  expect_lt(p_spike, 0.01)
  expect_gt(p_spike, 0)
})


test_that("parameter bounds are enforced", {
  # The density warns and returns zero where the RNG errors - the shared
  # convention for every family on this parameterization.
  expect_warning(expect_equal(dcogmod_logstudent(0.5, dof = 0), 0), "dof")
  expect_warning(expect_equal(dcogmod_logstudent(0.5, dof = -1), 0), "dof")
  expect_warning(expect_equal(dcogmod_logstudent(0.5, sigma = 0), 0), "sigma")
  expect_warning(dcogmod_logstudent(0.5, poutlier = 1.5), "poutlier")

  expect_error(rcogmod_logstudent(5, dof = 0), "greater than 0")
  expect_error(rcogmod_logstudent(5, sigma = 0), "greater than 0")
  expect_error(rcogmod_logstudent(5, ndt = -0.1), "ndt")
  # dof is strictly positive but otherwise unconstrained: below 1 is legal, and
  # is just a heavier tail than the Cauchy.
  expect_silent(rcogmod_logstudent(5, dof = 0.5))
})


test_that("the family declares dof and posterior_epred refuses to guess", {
  fam <- cogmod_logstudent()
  expect_equal(fam$dpars, c("mu", "sigma", "dof", "ndt", "poutlier"))
  expect_equal(fam$name, "cogmod_logstudent")
  expect_equal(cogmod_logstudent(link_dof = "softplus")$link_dof, "softplus")

  # No expectation exists anywhere in the parameter space.
  expect_error(posterior_epred_cogmod_logstudent(NULL), "no finite mean")
})


test_that("cogmod_priors fills dof, and brms leaves it flat", {
  set.seed(4)
  d <- data.frame(
    RT = rcogmod_logstudent(200, ndt = 0.2, poutlier = 0.02),
    Condition = factor(rep(c("a", "b"), 100))
  )
  f <- brms::bf(RT ~ Condition, dof ~ Condition, ndt ~ 1, poutlier ~ 1,
                family = cogmod_logstudent())

  # brms has no opinion about a dpar called `dof` - which is the reason it is
  # not called `nu` - so this is a fill rather than an override.
  raw <- brms::get_prior(f, d, family = f$family)
  expect_true(all(!nzchar(raw$prior[raw$dpar == "dof"])))

  p <- expect_silent(cogmod_priors(f, d))
  pick <- function(cls, coef = "") {
    p$prior[p$class == cls & p$dpar == "dof" & p$coef == coef]
  }
  expect_equal(pick("Intercept"), "normal(1.8, 0.7)")
  expect_equal(pick("b"), "normal(0, 0.3)")

  # Omitted from bf(), dof lives on the natural scale instead.
  f2 <- brms::bf(RT ~ Condition, family = cogmod_logstudent())
  p2 <- expect_silent(cogmod_priors(f2, d))
  expect_equal(p2$prior[p2$class == "dof"], "lognormal(1.8, 0.7)")

  # inits knows the family, and brms builds the program from the result. The
  # tolerance is for the jitter cogmod_inits() puts on every starting value.
  expect_equal(cogmod_inits(f, d)(1)$Intercept_dof, log(5), tolerance = 0.15)
  expect_silent(
    brms::make_stancode(f, data = d, family = f$family, prior = p,
                        stanvars = cogmod_stanvars(f))
  )
})


context("Log-Student-t - brms")

test_that("Stan cogmod_logstudent_lpdf matches R dcogmod_logstudent", {
  fun <- stan_fun("cogmod_logstudent")
  grid <- covering_grid(
    mu = c(-1.2, -0.7, 0),
    sigma = c(0.2, 0.4, 0.8),
    dof = c(1, 2, 5, 30, 200),
    ndt = c(0, 0.2),
    poutlier = c(0, 0.02, 0.3),
    always = function(g) g$dof == 1 & g$ndt == 0.2 & g$poutlier == 0.02
  )
  x <- c(0.05, 0.19, 0.21, 0.5, 1.5, 8)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    r <- dcogmod_logstudent(x, mu = g$mu, sigma = g$sigma, dof = g$dof,
                            ndt = g$ndt, poutlier = g$poutlier, log = TRUE)
    s <- vapply(x, function(xx) {
      fun(xx, g$mu, g$sigma, g$dof, g$ndt, g$poutlier)
    }, numeric(1))
    expect_equal(s, r, tolerance = 1e-8,
                 label = sprintf("mu=%s sigma=%s dof=%s ndt=%s p=%s",
                                 g$mu, g$sigma, g$dof, g$ndt, g$poutlier))
  }

  # Invalid parameters are rejected on both sides, identically.
  expect_equal(fun(0.5, -0.7, 0.4, 0, 0.2, 0.02), -Inf)
  expect_equal(fun(0.5, -0.7, 0, 5, 0.2, 0.02), -Inf)
  expect_equal(fun(-1, -0.7, 0.4, 5, 0.2, 0.02), -Inf)
})


test_that("Log-Student-t recovers its parameters with brms", {
  skip_if_not_slow()
  skip_if_not_installed("cmdstanr")

  set.seed(21)
  n <- 2000
  d <- data.frame(
    RT = rcogmod_logstudent(n, mu = -0.7, sigma = 0.4, dof = 4, ndt = 0.25,
                            poutlier = 0.02)
  )
  f <- brms::bf(RT ~ 1, sigma ~ 1, dof ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_logstudent())
  m <- brms::brm(
    f, data = d, prior = cogmod_priors(f, d), init = cogmod_inits(f, d),
    stanvars = cogmod_stanvars(f), chains = 2, iter = 1000,
    backend = "cmdstanr", refresh = 0
  )
  fx <- brms::fixef(m)
  expect_equal(unname(fx["Intercept", "Estimate"]), -0.7, tolerance = 0.15)
  expect_equal(unname(exp(fx["dof_Intercept", "Estimate"])), 4, tolerance = 3)
  expect_equal(unname(exp(fx["ndt_Intercept", "Estimate"])), 0.25,
               tolerance = 0.06)
  expect_true(all(brms::rhat(m)[!is.na(brms::rhat(m))] < 1.05))
})
