context("Shifted Log-Gamma - brms")

# Reference implementation of the decision density, written out longhand from
# the generalized-gamma definition rather than from the internals under test:
# T = exp(mu + sigma * w), w = log(shape^2 * g) / shape, g ~ Gamma(1 / shape^2, 1).
ref_dec <- function(t, mu, sigma, shape) {
  if (t <= 0) return(0)
  if (shape == 0) return(dlnorm(t, mu, sigma))
  k <- 1 / shape^2
  w <- (log(t) - mu) / sigma
  exp(-log(sigma) - log(t) + log(abs(shape)) + k * log(k) - lgamma(k) +
    k * (shape * w - exp(shape * w)))
}

ref_dens <- function(y, mu, sigma, shape, ndt, poutlier) {
  d_out <- 2 * dnorm(y, 0, 0.2)
  poutlier * d_out + (1 - poutlier) * ref_dec(y - ndt, mu, sigma, shape)
}


# dcogmod_loggamma ------------------------------------------------------------

test_that("dcogmod_loggamma matches the mixture density", {
  for (shape in c(-0.8, -0.3, 0, 0.3, 1, 2)) {
    for (y in c(0.05, 0.31, 0.5, 0.9, 1.5, 3)) {
      expect_equal(
        dcogmod_loggamma(y, mu = -0.7, sigma = 0.5, shape = shape, ndt = 0.3,
                     poutlier = 0.02),
        ref_dens(y, -0.7, 0.5, shape, 0.3, 0.02),
        tolerance = 1e-10,
        label = sprintf("density at shape = %.1f, y = %.2f", shape, y)
      )
    }
  }
})

test_that("shape = 0 is exactly the shifted LogNormal", {
  y <- c(0.05, 0.4, 0.9, 2)
  expect_equal(
    dcogmod_loggamma(y, mu = -0.7, sigma = 0.5, shape = 0, ndt = 0.3, poutlier = 0.02),
    dcogmod_lognormal(y, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02),
    tolerance = 1e-14
  )
  # and it is reached exactly, not approached: the log-density is continuous
  # through shape = 0 to well below what any sampler would resolve
  for (shape in c(1e-12, 1e-9, 1e-6, 1e-4)) {
    expect_equal(
      dcogmod_loggamma(y, -0.7, 0.5, shape, 0.3, 0, log = TRUE),
      dcogmod_lognormal(y, -0.7, 0.5, 0.3, 0, log = TRUE),
      tolerance = 5 * shape + 1e-12,
      label = sprintf("shape = %g", shape)
    )
  }
})

test_that("the nested special cases come out exactly", {
  t <- c(0.2, 0.5, 1, 2)
  mu <- -0.7
  s <- 0.5

  # shape = 1: Weibull with shape 1 / sigma, scale exp(mu)
  expect_equal(
    dcogmod_loggamma(t, mu, s, shape = 1, ndt = 0),
    dweibull(t, shape = 1 / s, scale = exp(mu)),
    tolerance = 1e-12
  )
  # shape = sigma: Gamma with shape 1 / sigma^2, scale exp(mu) * sigma^2
  expect_equal(
    dcogmod_loggamma(t, mu, s, shape = s, ndt = 0),
    dgamma(t, shape = 1 / s^2, scale = exp(mu) * s^2),
    tolerance = 1e-12
  )
  # shape = -1: inverse Weibull
  expect_equal(
    dcogmod_loggamma(t, mu, s, shape = -1, ndt = 0),
    (1 / s) * (exp(mu) / t)^(1 / s) / t * exp(-(exp(mu) / t)^(1 / s)),
    tolerance = 1e-12
  )
})

test_that("dcogmod_loggamma integrates to one across the shape range", {
  for (shape in c(-0.9, -0.5, 0, 0.5, 1, 2)) {
    total <- integrate(
      function(z) dcogmod_loggamma(z, shape = shape, ndt = 0.3, poutlier = 0.05),
      lower = 0, upper = Inf
    )$value
    expect_equal(total, 1, tolerance = 1e-4, label = sprintf("shape = %.1f", shape))
  }
})

test_that("responses faster than ndt keep positive density", {
  d <- dcogmod_loggamma(0.1, shape = 0.5, ndt = 0.4, poutlier = 0.02)
  expect_gt(d, 0)
  expect_equal(d, 0.02 * 2 * dnorm(0.1, 0, 0.2), tolerance = 1e-12)

  expect_equal(dcogmod_loggamma(0.1, shape = 0.5, ndt = 0.4, poutlier = 0), 0)
  expect_true(
    is.infinite(dcogmod_loggamma(0.1, shape = 0.5, ndt = 0.4, poutlier = 0, log = TRUE))
  )
})

test_that("dcogmod_loggamma is vectorized over every parameter", {
  d <- dcogmod_loggamma(
    c(0.5, 0.9), mu = c(-0.7, -0.5), sigma = 0.5, shape = c(0, 0.5), ndt = 0.3
  )
  expect_length(d, 2)
  expect_equal(
    d,
    c(ref_dec(0.2, -0.7, 0.5, 0), ref_dec(0.6, -0.5, 0.5, 0.5)),
    tolerance = 1e-12
  )
})

test_that("dcogmod_loggamma returns 0 density for invalid parameters", {
  expect_warning(d <- dcogmod_loggamma(0.5, sigma = -0.5))
  expect_equal(d, 0)
  expect_warning(d <- dcogmod_loggamma(0.5, ndt = -0.3))
  expect_equal(d, 0)
  expect_warning(d <- dcogmod_loggamma(0.5, poutlier = 1.5))
  expect_equal(d, 0)
})

# rcogmod_loggamma ------------------------------------------------------------

test_that("rcogmod_loggamma draws follow dcogmod_loggamma", {
  set.seed(123)
  for (shape in c(-0.5, 0, 0.5, 2)) {
    x <- rcogmod_loggamma(3e4, mu = -0.7, sigma = 0.5, shape = shape, ndt = 0.3)
    for (p in c(0.25, 0.5, 0.75)) {
      emp <- integrate(
        function(z) dcogmod_loggamma(z, shape = shape, ndt = 0.3),
        lower = 0, upper = unname(quantile(x, p))
      )$value
      expect_equal(emp, p, tolerance = 0.02,
                   label = sprintf("shape = %.1f, p = %.2f", shape, p))
    }
  }
})

test_that("rcogmod_loggamma recovers the mixture mean", {
  set.seed(1)
  mu <- -0.7
  sigma <- 0.5
  shape <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  rts <- rcogmod_loggamma(2e4, mu, sigma, shape, ndt, poutlier)
  theo <- (1 - poutlier) * (cogmod:::.mloggamma(mu, sigma, shape) + ndt) +
    poutlier * cogmod:::.mcontam()

  expect_equal(mean(rts), theo, tolerance = 0.02)
  expect_true(all(rts > 0))
  expect_length(rts, 2e4)
})

test_that("rcogmod_loggamma produces responses below ndt only via outliers", {
  set.seed(1)
  expect_true(all(rcogmod_loggamma(5000, shape = 0.5, ndt = 0.3, poutlier = 0) > 0.3))
  expect_true(any(rcogmod_loggamma(5000, shape = 0.5, ndt = 0.3, poutlier = 0.5) < 0.3))
})

test_that("rcogmod_loggamma errors on invalid parameters", {
  expect_error(rcogmod_loggamma(10, sigma = -1), "sigma")
  expect_error(rcogmod_loggamma(10, ndt = -1), "ndt")
  expect_error(rcogmod_loggamma(10, poutlier = 2), "poutlier")
  expect_error(rcogmod_loggamma(-5), "non-negative integer")
})


# analytic mean -----------------------------------------------------------

test_that(".mloggamma matches numerical integration, and is Inf where it must be", {
  for (shape in c(-0.9, -0.5, 0, 0.5, 1, 2)) {
    num <- integrate(
      function(z) z * cogmod:::.dloggamma(z, -0.7, 0.5, shape),
      lower = 0, upper = Inf
    )$value
    expect_equal(cogmod:::.mloggamma(-0.7, 0.5, shape), num, tolerance = 1e-5,
                 label = sprintf("shape = %.1f", shape))
  }

  # For shape < 0 the right tail is a power law, with no mean once sigma * |shape| >= 1
  expect_true(is.infinite(cogmod:::.mloggamma(-0.7, 0.5, -2)))
  expect_true(is.infinite(cogmod:::.mloggamma(-0.7, 1, -1)))
  expect_true(is.finite(cogmod:::.mloggamma(-0.7, 0.5, -1.9)))
  # ...but shape > 0 always has one
  expect_true(is.finite(cogmod:::.mloggamma(-0.7, 1, 5)))
})


# brms methods ------------------------------------------------------------

make_prep <- function(y, mu, sigma, shape, ndt, poutlier, n_draws = 10) {
  structure(
    list(
      data = list(Y = y),
      dpars = list(
        mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
        shape = rep(shape, n_draws), ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}

test_that("log_lik_cogmod_loggamma matches dcogmod_loggamma", {
  for (y in c(0.1, 0.35, 0.9, 2)) {
    prep <- make_prep(y, -0.7, 0.5, 0.5, 0.3, 0.02)
    expect_equal(
      log_lik_cogmod_loggamma(1, prep),
      rep(log(ref_dens(y, -0.7, 0.5, 0.5, 0.3, 0.02)), 10),
      tolerance = 1e-10,
      label = sprintf("log_lik at y = %.2f", y)
    )
  }
})

test_that("log_lik_cogmod_loggamma returns -Inf for invalid parameters", {
  prep <- structure(
    list(
      data = list(Y = 0.5),
      dpars = list(
        mu = rep(-0.7, 3), sigma = c(-0.5, 0.5, 0.5), shape = rep(0.5, 3),
        ndt = c(0.3, -0.3, 0.3), poutlier = c(0.02, 0.02, 1.5)
      )
    ),
    class = "brmsprep"
  )
  expect_true(all(suppressWarnings(log_lik_cogmod_loggamma(1, prep)) == -Inf))
})

test_that("posterior_predict_cogmod_loggamma excludes outliers by default", {
  set.seed(123)
  n_draws <- 1e5
  mu <- -0.7
  sigma <- 0.5
  shape <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  prep <- structure(
    list(dpars = list(
      mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
      shape = rep(shape, n_draws), ndt = rep(ndt, n_draws),
      poutlier = rep(poutlier, n_draws)
    )),
    class = "brmsprep"
  )

  rts <- posterior_predict_cogmod_loggamma(1, prep)
  expect_equal(mean(rts), cogmod:::.mloggamma(mu, sigma, shape) + ndt,
               tolerance = 0.02)
  expect_true(all(rts > ndt))

  rts_mix <- posterior_predict_cogmod_loggamma(1, prep, predict_outliers = TRUE)
  expect_true(any(rts_mix < ndt))
  expect_true(all(rts_mix > 0))
})

test_that("posterior_epred_cogmod_loggamma excludes outliers by default", {
  mat <- function(v) matrix(v, nrow = 2, ncol = 2)
  prep <- structure(
    list(dpars = list(
      mu = mat(-0.7), sigma = mat(0.5), shape = mat(0.5), ndt = mat(0.3),
      poutlier = mat(0.02)
    )),
    class = "brmsprep"
  )

  expect_equal(
    posterior_epred_cogmod_loggamma(prep),
    mat(cogmod:::.mloggamma(-0.7, 0.5, 0.5) + 0.3)
  )
  expect_equal(
    posterior_epred_cogmod_loggamma(prep, predict_outliers = TRUE),
    mat(0.98 * (cogmod:::.mloggamma(-0.7, 0.5, 0.5) + 0.3) +
          0.02 * cogmod:::.mcontam())
  )

  # shape = 0 gives back the LogNormal expectation exactly
  prep$dpars$shape <- mat(0)
  expect_equal(
    posterior_epred_cogmod_loggamma(prep),
    mat(exp(-0.7 + 0.5^2 / 2) + 0.3)
  )
})

test_that("posterior_epred_cogmod_loggamma keeps the draws x obs shape", {
  m <- matrix(c(-0.7, -0.5, -0.9, -0.6), nrow = 2)
  prep <- structure(
    list(dpars = list(
      mu = m, sigma = matrix(0.5, 2, 2), shape = matrix(c(0, 0.5, -0.3, 1), 2),
      ndt = matrix(0.3, 2, 2), poutlier = matrix(0.02, 2, 2)
    )),
    class = "brmsprep"
  )
  out <- posterior_epred_cogmod_loggamma(prep)
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(
    out[2, 1],
    cogmod:::.mloggamma(-0.5, 0.5, 0.5) + 0.3
  )
})

test_that("posterior_epred_cogmod_loggamma keeps the shape whichever dpar is constant", {
  # get_dpar() returns a bare scalar for a dpar constant across draws and
  # observations, so the draws x obs shape cannot be read off `mu` alone.
  mk <- function(...) {
    structure(list(dpars = list(...)), class = "brmsprep")
  }
  m <- function(v) matrix(v, 2, 3)
  cases <- list(
    list(mu = m(-0.7), sigma = m(0.5), shape = m(0.5), ndt = m(0.3),
         poutlier = m(0.02)),
    list(mu = -0.7, sigma = m(0.5), shape = 0.5, ndt = 0.3, poutlier = 0.02),
    list(mu = -0.7, sigma = 0.5, shape = 0.5, ndt = m(0.3), poutlier = 0.02),
    list(mu = m(-0.7), sigma = 0.5, shape = 0.5, ndt = 0.3, poutlier = 0.02),
    list(mu = -0.7, sigma = 0.5, shape = m(0.5), ndt = 0.3, poutlier = m(0.02))
  )
  ref <- m(cogmod:::.mloggamma(-0.7, 0.5, 0.5) + 0.3)
  for (i in seq_along(cases)) {
    prep <- do.call(mk, cases[[i]])
    expect_equal(posterior_epred_cogmod_loggamma(prep), ref, label = paste("case", i))
    expect_equal(
      dim(posterior_epred_cogmod_loggamma(prep, predict_outliers = TRUE)),
      c(2L, 3L),
      label = paste("case", i, "mixture")
    )
  }
})

test_that("the family flag drives predictions when no argument is given", {
  set.seed(123)
  n_draws <- 2e4
  dpars <- list(
    mu = rep(-0.7, n_draws), sigma = rep(0.5, n_draws),
    shape = rep(0.5, n_draws), ndt = rep(0.3, n_draws),
    poutlier = rep(0.2, n_draws)
  )
  on <- structure(
    list(dpars = dpars,
         family = list(name = "cogmod_loggamma", predict_outliers = TRUE)),
    class = "brmsprep"
  )
  off <- structure(
    list(dpars = dpars,
         family = list(name = "cogmod_loggamma", predict_outliers = FALSE)),
    class = "brmsprep"
  )

  expect_true(any(posterior_predict_cogmod_loggamma(1, on) < 0.3))
  expect_true(all(posterior_predict_cogmod_loggamma(1, off) > 0.3))
  expect_true(all(
    posterior_predict_cogmod_loggamma(1, on, predict_outliers = FALSE) > 0.3
  ))
})

# family ------------------------------------------------------------------

test_that("cogmod_loggamma() builds a valid brms custom family", {
  fam <- cogmod_loggamma()

  expect_s3_class(fam, "customfamily")
  expect_identical(fam$dpars, c("mu", "sigma", "shape", "ndt", "poutlier"))
  expect_identical(
    unname(c(fam$link, fam$link_sigma, fam$link_shape, fam$link_ndt,
             fam$link_poutlier)),
    c("identity", "softplus", "identity", "log", "logit")
  )
  # shape must be unconstrained: shape = 0 (the LogNormal) has to sit in the interior
  expect_true(is.na(fam$lb[["shape"]]))
  expect_true(is.na(fam$ub[["shape"]]))
  expect_false("minrt" %in% fam$dpars)
  expect_false(fam$predict_outliers)
  expect_true(cogmod_loggamma(predict_outliers = TRUE)$predict_outliers)
})

test_that("the family survives the validation brm() performs", {
  fam <- cogmod_loggamma(predict_outliers = TRUE)
  expect_true(brms:::validate_family(fam)$predict_outliers)

  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  f <- brms::bf(RT ~ 1, sigma ~ 1, shape ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  vf <- brms:::validate_formula(f, data = d)
  expect_true(vf$family$predict_outliers)
})

test_that("stanvars carry the likelihood, and brms builds Stan code with it", {
  sv <- cogmod_loggamma_stanvars()
  code <- paste(vapply(sv, function(x) x$scode, character(1)), collapse = "\n")

  expect_true(grepl("cogmod_loggamma_lpdf", code, fixed = TRUE))
  expect_true(grepl("cogmod_loggamma_lconst", code, fixed = TRUE))
  expect_true(grepl("cogmod_loggamma_lkernel", code, fixed = TRUE))
  expect_true(grepl("log_mix", code, fixed = TRUE))
  # the outlier component is the half Normal, folded into a literal
  expect_true(grepl("12.5 * square(Y)", code, fixed = TRUE))
  expect_true(grepl("half Normal with scale 0.2", code, fixed = TRUE))

  # the dpar order in the Stan signature must match the family's
  expect_true(grepl(
    "cogmod_loggamma_lpdf(real Y, real mu, real sigma, real shape, real ndt, real poutlier)",
    code, fixed = TRUE
  ))

  set.seed(1)
  d <- data.frame(RT = rcogmod_loggamma(50, shape = 0.3, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, sigma ~ 1, shape ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_loggamma())
  expect_silent(
    brms::make_stancode(f, data = d, family = cogmod_loggamma(),
                        prior = cogmod_priors(f, d), stanvars = sv)
  )
})


# shared machinery --------------------------------------------------------

test_that("with_outliers works on cogmod_loggamma fits", {
  fake <- structure(
    list(family = list(name = "cogmod_loggamma"),
         formula = list(family = list(name = "cogmod_loggamma"))),
    class = "brmsfit"
  )
  m <- with_outliers(fake)
  expect_true(m$formula$family$predict_outliers)
  expect_true(m$family$predict_outliers)
  expect_false(without_outliers(m)$formula$family$predict_outliers)
})

test_that("p_outlier uses the log-gamma decision density", {
  y <- c(0.05, 0.4, 0.9, 3.0) # below ndt / just after / bulk / far tail
  nd <- 4
  prep <- list(
    data = list(Y = y),
    family = list(name = "cogmod_loggamma"),
    dpars = list(
      mu = matrix(-0.9, nd, length(y)), sigma = matrix(0.5, nd, length(y)),
      shape = matrix(0.5, nd, length(y)), ndt = matrix(0.25, nd, length(y)),
      poutlier = matrix(0.02, nd, length(y))
    )
  )

  out <- with_mocked_bindings(
    p_outlier(structure(list(), class = "brmsfit")),
    prepare_predictions = function(object, ...) prep,
    get_dpar = function(prep, dpar, ...) prep$dpars[[dpar]],
    .package = "brms"
  )

  ref <- vapply(y, function(v) {
    g <- 2 * dnorm(v, 0, 0.2)
    f <- ref_dec(v - 0.25, -0.9, 0.5, 0.5)
    0.02 * g / (0.02 * g + 0.98 * f)
  }, numeric(1))

  expect_equal(out$p_outlier, ref, tolerance = 1e-10)
  expect_equal(out$p_outlier[1], 1)
  expect_lt(out$p_outlier[3], out$p_outlier[2])
})


# priors ------------------------------------------------------------------

test_that("cogmod_priors replaces every brms default that is wrong here", {
  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  f <- brms::bf(RT ~ 1, sigma ~ 1, shape ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_loggamma())

  raw <- brms::get_prior(f, data = d, family = cogmod_loggamma())
  # ndt and poutlier are left flat; `shape` is worse than flat, because brms
  # knows the name and supplies a positive-support prior for it
  expect_true(all(raw$prior[raw$class == "Intercept" &
                              raw$dpar %in% c("ndt", "poutlier")] == ""))
  expect_equal(raw$prior[raw$class == "Intercept" & raw$dpar == "shape"],
               "gamma(0.01, 0.01)")

  p <- cogmod_priors(f, d)
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "shape"],
               "normal(0, 0.5)")
  expect_true(all(nzchar(
    p$prior[p$class == "Intercept" & p$dpar %in% c("shape", "ndt", "poutlier")]
  )))
})

test_that("dpars omitted from bf() get a natural-scale prior, not a link-scale one", {
  # Omitting a dpar makes brms declare it as a plain auxiliary parameter -
  # class "<name>", dpar "", no link applied - instead of building a linear
  # predictor for it. Matching on dpar alone missed those rows entirely, which
  # left an omitted `poutlier` on brms's flat prior over [0, 1].
  set.seed(1)
  d <- data.frame(RT = rcogmod_loggamma(100, shape = 0.3, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, family = cogmod_loggamma())

  raw <- brms::get_prior(f, data = d, family = cogmod_loggamma())
  # what brms leaves behind, and why each is unacceptable
  expect_equal(raw$prior[raw$class == "shape"], "gamma(0.01, 0.01)") # positive only
  expect_equal(raw$prior[raw$class == "poutlier"], "")               # flat on [0, 1]
  expect_match(raw$prior[raw$class == "ndt"], "min_Y")               # min-RT bound

  p <- cogmod_priors(f, d)
  expect_equal(p$prior[p$class == "shape" & !nzchar(p$dpar)], "normal(0, 0.5)")
  expect_equal(p$prior[p$class == "poutlier" & !nzchar(p$dpar)],
               "exponential(100)")
  # lognormal on the natural scale is the same belief as normal on log(ndt)
  expect_equal(p$prior[p$class == "ndt" & !nzchar(p$dpar)],
               "lognormal(-1.2, 0.2)")

  code <- brms::make_stancode(f, data = d, family = cogmod_loggamma(), prior = p,
                              stanvars = cogmod_loggamma_stanvars())
  expect_true(grepl("normal_lpdf(shape | 0, 0.5)", code, fixed = TRUE))
  expect_true(grepl("exponential_lpdf(poutlier | 100)", code, fixed = TRUE))
  expect_true(grepl("lognormal_lpdf(ndt | -1.2, 0.2)", code, fixed = TRUE))
  # the min-RT bound is gone
  expect_false(grepl("uniform_lpdf(ndt", code, fixed = TRUE))
})

test_that("brms's own gamma(0.01, 0.01) default for `shape` is overridden", {
  # brms recognises the name `shape` from its gamma/weibull families and gives
  # it a prior supported only on the positives - for a parameter whose whole
  # negative half is the inverse-Weibull side of the family. Stan would reject
  # every negative proposal, truncating `shape` at 0 without saying so. The row
  # arrives non-empty, so it must be overridden, not merely filled.
  set.seed(1)
  d <- data.frame(RT = rcogmod_loggamma(100, shape = 0.3, ndt = 0.3, poutlier = 0.02),
                  x = rnorm(100))

  modelled <- brms::bf(RT ~ 1, sigma ~ 1, shape ~ 1, ndt ~ 1, poutlier ~ 1,
                       family = cogmod_loggamma())
  omitted <- brms::bf(RT ~ 1, family = cogmod_loggamma())

  # what brms would have done, in both forms
  raw_m <- brms::get_prior(modelled, data = d, family = cogmod_loggamma())
  raw_o <- brms::get_prior(omitted, data = d, family = cogmod_loggamma())
  expect_equal(raw_m$prior[raw_m$class == "Intercept" & raw_m$dpar == "shape"],
               "gamma(0.01, 0.01)")
  expect_equal(raw_o$prior[raw_o$class == "shape"], "gamma(0.01, 0.01)")

  for (f in list(modelled, omitted)) {
    p <- cogmod_priors(f, d)
    expect_true(all(p$prior[p$class == "shape" | p$dpar == "shape"] ==
                      "normal(0, 0.5)"))
    code <- brms::make_stancode(f, data = d, family = cogmod_loggamma(), prior = p,
                                stanvars = cogmod_loggamma_stanvars())
    expect_false(grepl("gamma_lpdf(shape | 0.01", code, fixed = TRUE))
    expect_false(grepl("gamma_lpdf(Intercept_shape | 0.01", code, fixed = TRUE))
    expect_true(grepl("normal_lpdf(shape | 0, 0.5)", code, fixed = TRUE) ||
                  grepl("normal_lpdf(Intercept_shape | 0, 0.5)", code, fixed = TRUE))
  }
})

test_that("modelled and omitted dpars can be mixed in one formula", {
  set.seed(1)
  d <- data.frame(RT = rcogmod_loggamma(100, shape = 0.3, ndt = 0.3, poutlier = 0.02),
                  x = rnorm(100))
  f <- brms::bf(RT ~ 1, shape ~ x, family = cogmod_loggamma())
  p <- cogmod_priors(f, d)

  # shape went through the regression path; ndt and poutlier did not
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "shape"],
               "normal(0, 0.5)")
  expect_equal(p$prior[p$class == "b" & p$dpar == "shape" & !nzchar(p$coef)],
               "normal(0, 0.2)")
  expect_equal(p$prior[p$class == "ndt"], "lognormal(-1.2, 0.2)")
  expect_equal(p$prior[p$class == "poutlier"], "exponential(100)")

  expect_silent(
    brms::make_stancode(f, data = d, family = cogmod_loggamma(), prior = p,
                        stanvars = cogmod_loggamma_stanvars())
  )
})

test_that("cogmod_priors survives arbitrary cogmod_loggamma formula shapes", {
  set.seed(1)
  d <- data.frame(
    RT = rcogmod_loggamma(200, shape = 0.3, ndt = 0.3, poutlier = 0.02),
    x = rnorm(200), g = factor(sample(letters[1:5], 200, TRUE))
  )
  fam <- cogmod_loggamma()
  forms <- list(
    plain = brms::bf(RT ~ 1, sigma ~ 1, shape ~ 1, ndt ~ 1, poutlier ~ 1,
                     family = fam),
    q_varies = brms::bf(RT ~ x, sigma ~ 1, shape ~ x, ndt ~ 1 + (1 | g),
                        poutlier ~ 1, family = fam),
    q_fixed = brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  )

  for (nm in names(forms)) {
    p <- cogmod_priors(forms[[nm]], d)
    expect_s3_class(p, "brmsprior")
    expect_silent(
      brms::make_stancode(forms[[nm]], data = d, family = fam, prior = p,
                          stanvars = cogmod_loggamma_stanvars())
    )
  }
})
