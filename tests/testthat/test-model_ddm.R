context("DDM - shifted Drift Diffusion Model with outliers")

# Reference implementation of the mixture. The decision component is delegated
# to brms::dwiener(), which is what the migration is *not* about; what is
# written out longhand here is the mixture around it and the 1/K on the outlier
# component.
ref_dens <- function(y, drift, boundary, bias, ndt, response, poutlier) {
  d_out <- 2 * stats::dnorm(y, 0, 0.2) / 2
  d_dec <- if (y - ndt > 0) {
    brms::dwiener(y, alpha = boundary, tau = ndt, beta = bias, delta = drift,
                  resp = response)
  } else {
    0
  }
  poutlier * d_out + (1 - poutlier) * d_dec
}

make_prep <- function(y, dec, drift, boundary, bias, ndt, poutlier,
                      sigmadrift = 0, sigmabias = 0, sigmandt = 0,
                      n_draws = 10, family = cogmod_ddm()) {
  structure(
    list(
      data = list(Y = y, dec = dec),
      family = family,
      dpars = list(
        mu = rep(drift, n_draws),
        boundary = rep(boundary, n_draws),
        bias = rep(bias, n_draws),
        sigmadrift = rep(sigmadrift, n_draws),
        sigmabias = rep(sigmabias, n_draws),
        sigmandt = rep(sigmandt, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}


# dcogmod_ddm -------------------------------------------------------------

test_that("dcogmod_ddm matches the mixture density", {
  pars <- list(drift = 0.5, boundary = 1, bias = 0.4, ndt = 0.2)
  for (poutlier in c(0, 0.02, 0.4)) {
    for (response in 0:1) {
      for (y in c(0.05, 0.15, 0.25, 0.6, 1.5, 4)) {
        expect_equal(
          do.call(dcogmod_ddm, c(list(x = y), pars,
                                 list(response = response,
                                      poutlier = poutlier))),
          do.call(ref_dens, c(list(y = y), pars,
                              list(response = response,
                                   poutlier = poutlier))),
          tolerance = 1e-10,
          info = sprintf("y = %.2f, response = %d, poutlier = %.2f",
                         y, response, poutlier)
        )
      }
    }
  }
})


test_that("the density sums to one over responses and integrates to one", {
  # Including the three between-trial variability parameters one at a time and
  # all together, since each takes a different branch of the density.
  grid <- list(
    c(drift = 0.5, boundary = 1, bias = 0.5, ndt = 0.2, sv = 0, sw = 0, st = 0),
    c(drift = -1.5, boundary = 1.6, bias = 0.35, ndt = 0.25, sv = 0, sw = 0,
      st = 0),
    c(drift = 0.8, boundary = 1, bias = 0.6, ndt = 0.2, sv = 0.8, sw = 0,
      st = 0),
    c(drift = 0.8, boundary = 1, bias = 0.5, ndt = 0.2, sv = 0, sw = 0.5,
      st = 0),
    c(drift = 0.8, boundary = 1, bias = 0.5, ndt = 0.2, sv = 0, sw = 0,
      st = 0.08),
    c(drift = 0.3, boundary = 1.2, bias = 0.45, ndt = 0.15, sv = 0.6, sw = 0.4,
      st = 0.06),
    c(drift = 0, boundary = 1, bias = 0.5, ndt = 0, sv = 0, sw = 0, st = 0)
  )
  for (g in grid) {
    for (poutlier in c(0, 0.03, 0.4)) {
      total <- sum(vapply(0:1, function(k) {
        stats::integrate(
          function(t) {
            dcogmod_ddm(t, g[["drift"]], g[["boundary"]], g[["bias"]],
                        g[["ndt"]], response = k, sigmadrift = g[["sv"]],
                        sigmabias = g[["sw"]], sigmandt = g[["st"]],
                        poutlier = poutlier)
          },
          lower = 0, upper = Inf, subdivisions = 4000
        )$value
      }, numeric(1)))
      expect_equal(total, 1, tolerance = 1e-4,
                   info = sprintf("sv=%.1f sw=%.1f st=%.2f poutlier=%.2f",
                                  g[["sv"]], g[["sw"]], g[["st"]], poutlier))
    }
  }
})


test_that("responses faster than ndt keep positive density (the whole point)", {
  y <- c(0.001, 0.05, 0.1, 0.19)
  for (poutlier in c(0.02, 0.4)) {
    for (response in 0:1) {
      d <- dcogmod_ddm(y, drift = 0.5, boundary = 1, bias = 0.5, ndt = 0.2,
                       response = response, poutlier = poutlier)
      expect_true(all(d > 0))
      # exactly the outlier component, thinned by the two response options
      expect_equal(d, poutlier * 0.5 * 2 * stats::dnorm(y, 0, 0.2),
                   tolerance = 1e-12)
    }
  }
})


test_that("poutlier = 0 recovers the plain shifted diffusion", {
  y <- c(0.05, 0.15, 0.2)
  expect_equal(dcogmod_ddm(y, drift = 0.5, ndt = 0.2, response = 1,
                           poutlier = 0),
               rep(0, length(y)))
  expect_equal(dcogmod_ddm(y, drift = 0.5, ndt = 0.2, response = 1,
                           poutlier = 0, log = TRUE),
               rep(-Inf, length(y)))
})


test_that("the tau0 offset is exactly that: an offset", {
  # brms::dwiener() refuses a zero non-decision time, but the density depends on
  # the time and the non-decision time only through their difference, so the
  # decision component is evaluated at (t + tau0, tau0). This checks that the
  # shift really does cancel: the mixture density at ndt must equal the plain
  # Wiener density with that same ndt.
  y <- c(0.21, 0.3, 0.8, 3)
  for (k in 0:1) {
    expect_equal(
      dcogmod_ddm(y, drift = 0.7, boundary = 1, bias = 0.4, ndt = 0.2,
                  response = k, poutlier = 0, log = TRUE),
      brms::dwiener(y, alpha = 1, tau = 0.2, beta = 0.4, delta = 0.7,
                    resp = k, log = TRUE),
      tolerance = 1e-12
    )
  }
})


test_that("dcogmod_ddm is vectorized over every argument", {
  n <- 7
  d <- dcogmod_ddm(seq(0.3, 1.2, length.out = n),
                   drift = seq(-1, 1, length.out = n),
                   boundary = seq(0.8, 1.6, length.out = n),
                   bias = seq(0.3, 0.7, length.out = n),
                   ndt = seq(0.05, 0.25, length.out = n),
                   response = rep(0:1, length.out = n),
                   poutlier = seq(0, 0.1, length.out = n))
  expect_length(d, n)
  expect_true(all(is.finite(d) & d >= 0))
})


test_that("dcogmod_ddm returns 0 density for invalid parameters", {
  args <- list(x = 0.5, drift = 0.5, boundary = 1, bias = 0.5, ndt = 0.2,
               response = 1)
  for (bad in list(list(boundary = 0), list(bias = 0), list(bias = 1),
                   list(sigmadrift = -1), list(sigmabias = 1),
                   list(ndt = -0.1), list(poutlier = 1.5),
                   list(response = 2))) {
    expect_warning(d <- do.call(dcogmod_ddm, modifyList(args, bad)),
                   info = names(bad))
    expect_equal(d, 0)
  }
})


test_that("the variability parameters are legitimately zero", {
  # That is the classic 4-parameter DDM, so their lower bounds are closed -
  # unlike `boundary` and `bias`, which are not.
  expect_silent(d <- dcogmod_ddm(0.5, drift = 0.5, response = 1,
                                 sigmadrift = 0, sigmabias = 0, sigmandt = 0))
  expect_true(is.finite(d) && d > 0)
  expect_silent(rcogmod_ddm(10, drift = 0.5, sigmadrift = 0, sigmabias = 0,
                            sigmandt = 0))
})


# rcogmod_ddm -------------------------------------------------------------

test_that("rcogmod_ddm returns rt and response", {
  set.seed(1)
  sim <- rcogmod_ddm(500, drift = 0.5, boundary = 1, bias = 0.5, ndt = 0.2)
  expect_true(is.data.frame(sim))
  expect_named(sim, c("rt", "response"))
  expect_equal(nrow(sim), 500)
  expect_true(all(sim$response %in% c(0, 1)))
  expect_true(all(is.finite(sim$rt)))
  expect_true(all(sim$rt > 0.2))
  # a positive drift pushes towards the boundary coded 1
  expect_gt(mean(sim$response == 1), 0.5)
})


test_that("pcogmod_ddm integrates dcogmod_ddm", {
  # The CDF is what rcogmod_ddm() inverts to draw from, so an error here would
  # land straight in the draws. Checked against the package's own density rather
  # than against rtdists, which is the less accurate of the two: at drift -4,
  # boundary 0.6, bias 0.25 it deviates from this integral by up to 5e-4 where
  # pcogmod_ddm() is at rounding error. See ?rcogmod_ddm.
  grid <- covering_grid(
    drift = c(-3, -1, 0, 1, 3),
    boundary = c(0.6, 1.2, 2.0),
    bias = c(0.3, 0.5, 0.7),
    q = c(0.35, 0.6, 1.2, 3),
    response = c(0L, 1L)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    integrated <- stats::integrate(
      function(t) dcogmod_ddm(t, drift = g$drift, boundary = g$boundary,
                              bias = g$bias, ndt = 0.2, response = g$response,
                              poutlier = 0),
      0, g$q, rel.tol = 1e-11, subdivisions = 2000
    )$value
    expect_equal(
      pcogmod_ddm(g$q, drift = g$drift, boundary = g$boundary, bias = g$bias,
                  ndt = 0.2, response = g$response, poutlier = 0),
      integrated,
      tolerance = 1e-8,
      info = sprintf("drift %.1f boundary %.1f bias %.2f q %.2f resp %d",
                     g$drift, g$boundary, g$bias, g$q, g$response)
    )
  }
})


test_that("pcogmod_ddm behaves like a CDF", {
  q <- seq(0.05, 6, length.out = 200)
  p <- pcogmod_ddm(q, drift = -1.5, boundary = 1.1, bias = 0.45, ndt = 0.2,
                   poutlier = 0.05)
  expect_false(is.unsorted(p))
  expect_true(all(p >= 0 & p <= 1))
  # The two defective CDFs partition the marginal one.
  expect_equal(pcogmod_ddm(1.5, response = 0) + pcogmod_ddm(1.5, response = 1),
               pcogmod_ddm(1.5))
  expect_equal(pcogmod_ddm(0.8, lower.tail = FALSE), 1 - pcogmod_ddm(0.8))
  expect_equal(pcogmod_ddm(0.8, log.p = TRUE), log(pcogmod_ddm(0.8)))
  # Nothing arrives before the non-decision time without an outlier process.
  expect_equal(pcogmod_ddm(0.1, ndt = 0.2, poutlier = 0), 0)
  expect_equal(pcogmod_ddm(500), 1)
  expect_length(pcogmod_ddm(c(0.3, 0.5, 0.9), drift = c(-1, 0, 1)), 3)
  # The between-trial variability parameters would each need a quadrature layer,
  # so they are not arguments at all rather than being silently ignored.
  expect_error(pcogmod_ddm(1, sigmadrift = 0.5))
})


test_that("pcogmod_ddm is the choice probability at infinite time", {
  # The series can only approach the choice probability - every term is
  # exp(-R t) - so infinite time is taken analytically rather than from it, and
  # used to be left at the pre-allocated zero instead.
  g <- list(drift = -4, boundary = 0.6, bias = 0.25, ndt = 0.2, poutlier = 0)
  for (resp in c(0L, 1L)) {
    expect_equal(
      do.call(pcogmod_ddm, c(list(Inf, response = resp), g)),
      do.call(pcogmod_ddm, c(list(1e6, response = resp), g)),
      info = sprintf("response %d", resp)
    )
  }
  expect_equal(do.call(pcogmod_ddm, c(list(Inf), g)), 1)
  expect_equal(do.call(pcogmod_ddm, c(list(Inf, lower.tail = FALSE), g)), 0)
  # An infinite time among finite ones must not disturb them, and -Inf is still
  # before any response.
  expect_equal(
    do.call(pcogmod_ddm, c(list(c(0.5, Inf, -Inf), response = 0L), g)),
    c(do.call(pcogmod_ddm, c(list(0.5, response = 0L), g)),
      do.call(pcogmod_ddm, c(list(Inf, response = 0L), g)), 0)
  )
})


test_that("pcogmod_ddm agrees with the sampler it is inverted from", {
  skip_on_cran()
  set.seed(4)
  for (g in list(c(-2, 1.2, 0.5), c(1.5, 0.8, 0.35), c(0, 1.5, 0.6))) {
    sim <- rcogmod_ddm(60000, drift = g[1], boundary = g[2], bias = g[3],
                       ndt = 0.2, poutlier = 0)
    probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    qq <- unname(stats::quantile(sim$rt, probs))
    expect_equal(
      pcogmod_ddm(qq, drift = g[1], boundary = g[2], bias = g[3], ndt = 0.2,
                  poutlier = 0),
      probs,
      tolerance = 0.01,
      info = sprintf("drift %.1f boundary %.1f bias %.2f", g[1], g[2], g[3])
    )
  }
})


test_that("rcogmod_ddm reproduces its own density", {
  set.seed(42)
  n <- 40000
  for (g in list(
    list(v = 0.5, a = 1, w = 0.5, ndt = 0.2, sv = 0, sw = 0, st = 0),
    list(v = 0.3, a = 1.2, w = 0.45, ndt = 0.15, sv = 0.6, sw = 0.4, st = 0.06)
  )) {
    for (poutlier in c(0, 0.05)) {
      sim <- rcogmod_ddm(n, drift = g$v, boundary = g$a, bias = g$w,
                         ndt = g$ndt, sigmadrift = g$sv, sigmabias = g$sw,
                         sigmandt = g$st, poutlier = poutlier)
      dens <- function(t, k) {
        dcogmod_ddm(t, drift = g$v, boundary = g$a, bias = g$w, ndt = g$ndt,
                    response = k, sigmadrift = g$sv, sigmabias = g$sw,
                    sigmandt = g$st, poutlier = poutlier)
      }
      for (k in 0:1) {
        integrated <- stats::integrate(function(t) dens(t, k), 0, Inf,
                                       subdivisions = 4000)$value
        expect_equal(mean(sim$response == k), integrated, tolerance = 0.015,
                     info = sprintf("sv %.1f, response %d, poutlier %.2f",
                                    g$sv, k, poutlier))
      }
      for (pr in c(0.25, 0.5, 0.9)) {
        qq <- unname(stats::quantile(sim$rt, pr))
        cdf <- sum(vapply(0:1, function(k) {
          stats::integrate(function(t) dens(t, k), 0, qq,
                           subdivisions = 4000)$value
        }, numeric(1)))
        expect_equal(cdf, pr, tolerance = 0.015,
                     info = sprintf("sv %.1f, p %.2f, poutlier %.2f",
                                    g$sv, pr, poutlier))
      }
    }
  }
})


test_that("ndt is a genuine shift, not a fraction of anything observed", {
  # The case the old tau * minrt parameterization could not express: a
  # non-decision time well above what any sample minimum would allow.
  set.seed(9)
  sim <- rcogmod_ddm(20000, drift = 0.5, boundary = 1, bias = 0.5, ndt = 0.9,
                     poutlier = 0)
  expect_gt(min(sim$rt), 0.9)
})


test_that("rcogmod_ddm produces responses below ndt only via outliers", {
  set.seed(2)
  expect_true(all(rcogmod_ddm(2000, drift = 0.5, ndt = 0.4,
                              poutlier = 0)$rt > 0.4))
  sim <- rcogmod_ddm(20000, drift = 0.5, ndt = 0.4, poutlier = 0.3)
  expect_true(any(sim$rt < 0.4))
  fast <- sim$response[sim$rt < 0.02]
  expect_setequal(unique(fast), c(0, 1))
  expect_equal(mean(fast == 0), 0.5, tolerance = 0.1)
})


test_that("rcogmod_ddm errors on invalid parameters", {
  expect_error(rcogmod_ddm(10, boundary = 0), "boundary")
  expect_error(rcogmod_ddm(10, bias = 0), "bias")
  expect_error(rcogmod_ddm(10, sigmadrift = -1), "sigmadrift")
  expect_error(rcogmod_ddm(10, ndt = -1), "ndt")
  expect_error(rcogmod_ddm(10, poutlier = 2), "poutlier")
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_ddm matches dcogmod_ddm", {
  prep <- make_prep(y = c(0.5, 0.15), dec = c(1, 0), drift = 0.5, boundary = 1,
                    bias = 0.4, ndt = 0.2, poutlier = 0.03)
  for (i in 1:2) {
    expect_equal(
      log_lik_cogmod_ddm(i, prep),
      rep(dcogmod_ddm(prep$data$Y[i], 0.5, 1, 0.4, 0.2,
                      response = prep$data$dec[i], poutlier = 0.03,
                      log = TRUE), 10)
    )
  }
})


test_that("log_lik_cogmod_ddm uses the observed choice", {
  p0 <- make_prep(y = 0.5, dec = 0, drift = 1.5, boundary = 1, bias = 0.7,
                  ndt = 0.2, poutlier = 0)
  p1 <- make_prep(y = 0.5, dec = 1, drift = 1.5, boundary = 1, bias = 0.7,
                  ndt = 0.2, poutlier = 0)
  expect_false(isTRUE(all.equal(log_lik_cogmod_ddm(1, p0),
                                log_lik_cogmod_ddm(1, p1))))
})


test_that("log_lik_cogmod_ddm errors when dec is missing", {
  prep <- make_prep(y = 0.5, dec = 1, drift = 0.5, boundary = 1, bias = 0.5,
                    ndt = 0.2, poutlier = 0)
  prep$data$dec <- NULL
  expect_error(log_lik_cogmod_ddm(1, prep), "dec")
})


test_that("posterior_predict_cogmod_ddm excludes outliers by default", {
  set.seed(3)
  prep <- make_prep(y = 0.5, dec = 1, drift = 0.5, boundary = 1, bias = 0.5,
                    ndt = 0.4, poutlier = 0.9, n_draws = 500)
  out <- posterior_predict_cogmod_ddm(1, prep)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(500L, 2L))
  expect_true(all(out[, 1] > 0.4))
  expect_true(all(out[, 2] %in% c(0, 1)))

  with_out <- posterior_predict_cogmod_ddm(1, prep, predict_outliers = TRUE)
  expect_true(any(with_out[, 1] < 0.4))
})


test_that("the family flag drives predictions when no argument is given", {
  set.seed(4)
  prep <- make_prep(y = 0.5, dec = 1, drift = 0.5, boundary = 1, bias = 0.5,
                    ndt = 0.4, poutlier = 0.9, n_draws = 500,
                    family = cogmod_ddm(predict_outliers = TRUE))
  expect_true(any(posterior_predict_cogmod_ddm(1, prep)[, 1] < 0.4))
  expect_true(all(
    posterior_predict_cogmod_ddm(1, prep, predict_outliers = FALSE)[, 1] > 0.4
  ))
})


test_that("posterior_epred_cogmod_ddm returns the mean RT of the diffusion", {
  set.seed(5)
  prep <- make_prep(y = 0.5, dec = 1, drift = 0.5, boundary = 1, bias = 0.4,
                    ndt = 0.2, poutlier = 0.05, n_draws = 4)
  ep <- posterior_epred_cogmod_ddm(prep)
  sim <- rcogmod_ddm(200000, drift = 0.5, boundary = 1, bias = 0.4, ndt = 0.2,
                     poutlier = 0)
  expect_equal(unname(ep[1]), mean(sim$rt), tolerance = 0.01)

  # ... and the outlier component only when asked for, which pulls it towards
  # the half Normal's own mean
  ep_out <- posterior_epred_cogmod_ddm(prep, predict_outliers = TRUE)
  expect_equal(unname(ep_out[1]),
               0.95 * unname(ep[1]) + 0.05 * cogmod:::.mcontam(),
               tolerance = 1e-10)
})


test_that("posterior_epred_cogmod_ddm warns when the closed form is not exact", {
  prep <- make_prep(y = 0.5, dec = 1, drift = 0.5, boundary = 1, bias = 0.5,
                    ndt = 0.2, poutlier = 0, sigmadrift = 0.5)
  expect_warning(posterior_epred_cogmod_ddm(prep), "4-parameter")
})


test_that("posterior_epred_cogmod_ddm handles a drift at zero", {
  # The closed form has a removable singularity there.
  prep <- make_prep(y = 0.5, dec = 1, drift = 0, boundary = 1, bias = 0.4,
                    ndt = 0.2, poutlier = 0)
  ep <- posterior_epred_cogmod_ddm(prep)
  expect_true(all(is.finite(ep)))
  z <- 0.4 * 1
  expect_equal(unname(ep[1]), 0.2 + z * (1 - z), tolerance = 1e-8)
})


# family ------------------------------------------------------------------

test_that("cogmod_ddm() builds a valid brms custom family", {
  fam <- cogmod_ddm()
  expect_s3_class(fam, "customfamily")
  # Order matters: brms passes the dpars to cogmod_ddm_lpdf in exactly this one.
  expect_equal(fam$dpars, c("mu", "boundary", "bias", "sigmadrift",
                            "sigmabias", "sigmandt", "ndt", "poutlier"))
  expect_equal(unname(cogmod:::.family_links(fam)),
               c("identity", "softplus", "logit", "softplus", "logit", "log",
                 "log", "logit"))
  expect_equal(fam$vars, "dec[n]")
  expect_false(fam$predict_outliers)
  # the old parameterization is gone, sigmatau with it
  expect_false(any(c("tau", "minrt", "sigmatau") %in% fam$dpars))
})


test_that("cogmod_ddm is one of the outlier-mixture families", {
  expect_true("cogmod_ddm" %in% cogmod:::.OUTLIER_FAMILIES)
  expect_true("cogmod_ddm" %in% cogmod:::.cogmod_families())
  expect_length(cogmod:::.cogmod_families(),
                length(unique(cogmod:::.cogmod_families())))
})


test_that("predict_outliers survives the validation brm() performs", {
  fam <- brms:::validate_family(cogmod_ddm(predict_outliers = TRUE))
  expect_true(fam$predict_outliers)
})


# Stan code ---------------------------------------------------------------

test_that("stanvars carry the likelihood with the outlier component", {
  code <- cogmod_ddm_stanvars()[[1]]$scode
  expect_true(grepl("real cogmod_ddm_lpdf", code))
  expect_true(grepl("int dec", code))
  expect_true(grepl("log_mix\\(poutlier", code))
  expect_true(grepl("half Normal with scale 0.2", code))
  expect_true(grepl("real cogmod_ddm_decision_lpdf", code))
  # all three wiener_lpdf branches survive the change
  expect_true(grepl("wiener_lpdf\\(y \\| boundary, tau0, w, v\\)", code))
  expect_true(grepl("wiener_lpdf\\(y \\| boundary, tau0, w, v, sigmadrift\\)",
                    code))
  expect_true(grepl(
    "wiener_lpdf\\(y \\| boundary, tau0, w, v, sigmadrift, sw, sigmandt\\)",
    code
  ))
  # the old parameterization is gone
  expect_false(grepl("real tau,", code))
  expect_false(grepl("real minrt", code))
  expect_false(grepl("sigmatau", code))
})


test_that("the Stan offset is generated from the R one", {
  # They have to be the same number, so only one of them is written down.
  code <- cogmod_ddm_stanvars()[[1]]$scode
  expect_true(grepl(sprintf("real tau0 = %s;",
                            formatC(cogmod:::.DDM_TAU0, format = "g",
                                    digits = 17, width = 1)),
                    code, fixed = TRUE))
})


test_that("Stan cogmod_ddm_lpdf matches dcogmod_ddm", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_ddm")
  grid <- covering_grid(
    Y = c(0.02, 0.25, 0.5, 1.2, 4),
    mu = c(-1.5, 0, 0.8),
    boundary = c(0.6, 1.2, 2),
    bias = c(0.3, 0.5, 0.7),
    sigmadrift = c(0, 0.8),
    sigmabias = c(0, 0.5),
    sigmandt = c(0, 0.08),
    ndt = c(0, 0.15, 0.4),
    poutlier = c(0, 0.001, 0.4),
    dec = 0:1,
    # the three Stan branches - 4-, 5- and 7-parameter - all have to be swept,
    # and the 7-parameter one is where the quadrature lives
    always = function(g) {
      g$boundary == 1.2 & g$bias == 0.5 & g$ndt == 0.15 & g$poutlier == 0.001
    }
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    stan <- lpdf(g$Y, g$mu, g$boundary, g$bias, g$sigmadrift, g$sigmabias,
                 g$sigmandt, g$ndt, g$poutlier, as.integer(g$dec))
    r <- dcogmod_ddm(g$Y, g$mu, g$boundary, g$bias, g$ndt, response = g$dec,
                     sigmadrift = g$sigmadrift, sigmabias = g$sigmabias,
                     sigmandt = g$sigmandt, poutlier = g$poutlier, log = TRUE)
    if (is.infinite(r)) {
      expect_equal(stan, r)
    } else {
      # Relative, not absolute. The floor is looser than the closed-form
      # families manage because the 7-parameter branch is adaptive quadrature on
      # the Stan side and fixed-node Gauss-Legendre on the R side: the two
      # integrate the same function by different rules.
      expect_lt(abs(stan - r) / max(1, abs(r)), 1e-4)
    }
  }

  # invalid arguments are rejected on both sides
  expect_equal(lpdf(0.5, 0.5, 0, 0.5, 0, 0, 0, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0, 0, 0, 0, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 1, 0, 0, 0, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, -1, 0, 0, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, 0, 1, 0, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, 0, 0, -1, 0.2, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, 0, 0, 0, -0.1, 0.02, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, 0, 0, 0, 0.2, 1.5, 1L), -Inf)
  expect_equal(lpdf(0.5, 0.5, 1, 0.5, 0, 0, 0, 0.2, 0.02, 2L), -Inf)
  expect_equal(lpdf(0, 0.5, 1, 0.5, 0, 0, 0, 0.2, 0.02, 1L), -Inf)
})


# priors and inits --------------------------------------------------------

test_that("cogmod_priors fills ndt and poutlier for cogmod_ddm", {
  set.seed(6)
  sim <- rcogmod_ddm(150, drift = 0.5, ndt = 0.2, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150),
                  id = rep(1:10, length.out = 150))

  modelled <- brms::bf(RT | dec(Error) ~ 1, boundary ~ 1, bias ~ 1,
                       sigmadrift = 0, sigmabias = 0, sigmandt = 0, ndt ~ 1,
                       poutlier ~ 1, family = cogmod_ddm())
  p <- cogmod_priors(modelled, d)
  expect_true(any(p$dpar == "ndt" & p$class == "Intercept" &
                    p$prior == "normal(-1.2, 0.2)"))
  expect_true(any(p$dpar == "poutlier" & p$class == "Intercept" &
                    p$prior == "normal(-5, 1)"))

  omitted <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_ddm())
  p2 <- cogmod_priors(omitted, d)
  expect_false(any(grepl("uniform", p2$prior)))
  expect_true(any(p2$class == "ndt" & p2$prior == "lognormal(-1.2, 0.2)"))
  expect_true(any(p2$class == "poutlier" & p2$prior == "exponential(100)"))

  mixed <- brms::bf(RT | dec(Error) ~ Condition + (1 | id), ndt ~ Condition,
                    sigmadrift = 0, sigmabias = 0, sigmandt = 0,
                    family = cogmod_ddm())
  p3 <- cogmod_priors(mixed, d)
  expect_s3_class(p3, "brmsprior")
  code <- brms::make_stancode(mixed, data = d, prior = p3,
                              stanvars = cogmod_stanvars(mixed))
  expect_false(grepl("uniform\\(", code))
  expect_true(grepl("cogmod_ddm_lpdf", code))
})


test_that("cogmod_priors fences off the three variability parameters", {
  # Each has a floor at zero that its link reaches only at minus infinity, and
  # the likelihood stops changing well before then - the same improper-posterior
  # failure as cogmod_lba1()'s sigmabias. brms leaves all three flat.
  set.seed(11)
  sim <- rcogmod_ddm(200, drift = 0.5, ndt = 0.2, poutlier = 0.02)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 200))
  f <- brms::bf(RT | dec(Error) ~ Condition, sigmadrift ~ Condition,
                sigmabias ~ Condition, sigmandt ~ Condition, ndt ~ 1,
                family = cogmod_ddm())

  raw <- brms::get_prior(f, data = d, family = cogmod_ddm())
  flat <- raw$dpar %in% c("sigmadrift", "sigmabias", "sigmandt")
  expect_true(all(raw$prior[flat] == ""))

  p <- cogmod_priors(f, d)
  pick <- function(dp, cls) {
    p$prior[p$dpar == dp & p$class == cls & nzchar(p$prior)]
  }
  expect_equal(pick("sigmadrift", "Intercept"), "normal(0, 1)")
  expect_equal(pick("sigmabias", "Intercept"), "normal(-2, 1)")
  expect_equal(pick("sigmandt", "Intercept"), "normal(-3, 1)")
  for (dp in c("sigmadrift", "sigmabias", "sigmandt")) {
    expect_equal(pick(dp, "b"), "normal(0, 0.5)", label = paste(dp, "slope"))
  }

  # omitted, they live on their natural scales, which differ from one another
  f2 <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_ddm())
  p2 <- cogmod_priors(f2, d)
  expect_equal(p2$prior[p2$class == "sigmadrift"], "lognormal(-1, 0.75)")
  expect_equal(p2$prior[p2$class == "sigmabias"], "beta(1, 5)")
  expect_equal(p2$prior[p2$class == "sigmandt"], "lognormal(-3, 1)")

  for (form in list(f, f2)) {
    expect_silent(
      brms::make_stancode(form, data = d, family = cogmod_ddm(),
                          prior = cogmod_priors(form, d),
                          stanvars = cogmod_stanvars(form))
    )
  }
})


test_that("cogmod_inits covers the declared parameters", {
  set.seed(8)
  sim <- rcogmod_ddm(150, drift = 0.5, ndt = 0.2, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150))
  f <- brms::bf(RT | dec(Error) ~ Condition, boundary ~ 1, bias ~ 1,
                sigmadrift ~ 1, sigmabias ~ 1, sigmandt ~ 1, ndt ~ 1,
                poutlier ~ 1, family = cogmod_ddm())
  vals <- cogmod_inits(f, d)(1)
  code <- suppressWarnings(
    brms::make_stancode(f, data = d, family = cogmod_ddm())
  )
  declared <- vapply(cogmod:::.stan_param_decls(code), `[[`, character(1),
                     "name")
  expect_setequal(names(vals), declared)
  expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))))

  v0 <- cogmod_inits(f, d, jitter = 0)(1)
  expect_lt(exp(v0$Intercept_ndt), 0.3)
  expect_equal(log1p(exp(v0$Intercept_boundary)), 1, tolerance = 1e-8)
  expect_equal(plogis(v0$Intercept_bias), 0.5, tolerance = 1e-8)
  # sigmandt is in seconds now, on a log link, so it starts at a plausible
  # range rather than at a fraction of a parameter that no longer exists
  expect_equal(exp(v0$Intercept_sigmandt), 0.05, tolerance = 1e-8)
})


# a real fit --------------------------------------------------------------

test_that("cogmod_ddm recovers ndt above the fastest observed response", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n <- 2000
  sim <- rcogmod_ddm(n, drift = 0.8, boundary = 1.2, bias = 0.45, ndt = 0.25,
                     poutlier = 0.04)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  expect_lt(min(d$RT), 0.25)

  f <- brms::bf(RT | dec(Error) ~ 1, boundary ~ 1, bias ~ 1,
                sigmadrift = 0, sigmabias = 0, sigmandt = 0, ndt ~ 1,
                poutlier ~ 1, family = cogmod_ddm())
  fit <- brms::brm(f, data = d, prior = cogmod_priors(f, d),
                   init = cogmod_inits(f, d), stanvars = cogmod_stanvars(f),
                   algorithm = "pathfinder", refresh = 0,
                   backend = "cmdstanr")

  s <- brms::posterior_summary(fit)[, "Estimate"]
  expect_equal(s[["b_Intercept"]], 0.8, tolerance = 0.25)
  expect_equal(log1p(exp(s[["b_boundary_Intercept"]])), 1.2, tolerance = 0.2)
  expect_equal(plogis(s[["b_bias_Intercept"]]), 0.45, tolerance = 0.08)
  expect_equal(exp(s[["b_ndt_Intercept"]]), 0.25, tolerance = 0.05)
  expect_equal(plogis(s[["b_poutlier_Intercept"]]), 0.04, tolerance = 0.04)
  expect_gt(exp(s[["b_ndt_Intercept"]]), min(d$RT) * 10)

  ll <- brms::log_lik(fit, ndraws = 10)
  expect_equal(dim(ll), c(10L, n))
  expect_true(all(is.finite(ll)))

  po <- p_outlier(fit)
  expect_equal(nrow(po), n)
  expect_true(all(po$p_outlier[po$rt < 0.1] > 0.99))
  expect_equal(mean(po$p_outlier), 0.04, tolerance = 0.04)

  expect_true(with_outliers(fit)$family$predict_outliers)
  expect_false(without_outliers(fit)$family$predict_outliers)
})


test_that("the defective tails of pcogmod_ddm add to the response probability", {
  # `lower.tail = FALSE` used to return `1 - P(RT <= q, choice = k)`, which is
  # `P(RT > q OR choice != k)` - at q = Inf it gave the probability of the
  # *other* response instead of zero. It is now the defective survival
  # `P(RT > q, choice = k)`, integrated on the same footing as the lower tail,
  # so the two add to that response's own mass. Same convention as
  # `pcogmod_rdm()`.
  pars <- list(drift = 0.8, boundary = 1, bias = 0.45, ndt = 0.2)
  for (poutlier in c(0, 0.1)) {
    pk <- vapply(0:1, function(k) {
      do.call(pcogmod_ddm, c(list(q = Inf, response = k), pars,
                             list(poutlier = poutlier)))
    }, numeric(1))
    expect_equal(sum(pk), 1, tolerance = 1e-10)

    for (k in 0:1) {
      # the survival of a response is exhausted at infinite time
      expect_equal(
        do.call(pcogmod_ddm, c(list(q = Inf, response = k,
                                    lower.tail = FALSE), pars,
                               list(poutlier = poutlier))),
        0, tolerance = 1e-10
      )
      for (q in c(0.3, 1, 5)) {
        lo <- do.call(pcogmod_ddm, c(list(q = q, response = k), pars,
                                     list(poutlier = poutlier)))
        hi <- do.call(pcogmod_ddm, c(list(q = q, response = k,
                                          lower.tail = FALSE), pars,
                                     list(poutlier = poutlier)))
        expect_equal(lo + hi, pk[k + 1], tolerance = 1e-12,
                     info = sprintf("k = %d, q = %.1f, poutlier = %.2f",
                                    k, q, poutlier))
      }
    }
  }

  # The marginal is untouched: there the attainable mass is one.
  for (q in c(0.4, 1, 3)) {
    lo <- do.call(pcogmod_ddm, c(list(q = q), pars, list(poutlier = 0.05)))
    hi <- do.call(pcogmod_ddm, c(list(q = q, lower.tail = FALSE), pars,
                                 list(poutlier = 0.05)))
    expect_equal(lo + hi, 1, tolerance = 1e-12)
  }

  # And the lower tail still integrates the defective density.
  for (k in 0:1) {
    for (q in c(0.4, 0.8, 2)) {
      expect_equal(
        do.call(pcogmod_ddm, c(list(q = q, response = k), pars,
                               list(poutlier = 0.05))),
        stats::integrate(function(z) {
          do.call(dcogmod_ddm, c(list(x = z, response = k), pars,
                                 list(poutlier = 0.05)))
        }, lower = 0, upper = q, rel.tol = 1e-10)$value,
        tolerance = 1e-8
      )
    }
  }
})
