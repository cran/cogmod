context("LBA2 - shifted two-accumulator LBA with outliers")

# Reference implementation of the mixture, written out longhand from Brown &
# Heathcote (2008). Two things it makes explicit: the 1/K on the outlier
# component, and the conditioning on at least one positive drift.
ref_dens <- function(y, driftzero, driftone, sigmazero, sigmaone, sigmabias,
                     boundary, ndt, response, poutlier) {
  d_out <- 2 * stats::dnorm(y, 0, 0.2) / 2
  v_win <- if (response == 0) driftzero else driftone
  s_win <- if (response == 0) sigmazero else sigmaone
  v_los <- if (response == 0) driftone else driftzero
  s_los <- if (response == 0) sigmaone else sigmazero

  A <- sigmabias
  b <- sigmabias + boundary
  t <- y - ndt
  d_dec <- if (t > 0) {
    n1 <- function(v, s) (b - A - v * t) / (t * s)
    n2 <- function(v, s) (b - v * t) / (t * s)
    f <- (1 / A) * (-v_win * stats::pnorm(n1(v_win, s_win)) +
                      s_win * stats::dnorm(n1(v_win, s_win)) +
                      v_win * stats::pnorm(n2(v_win, s_win)) -
                      s_win * stats::dnorm(n2(v_win, s_win)))
    cdf <- 1 +
      ((b - A - v_los * t) / A) * stats::pnorm(n1(v_los, s_los)) -
      ((b - v_los * t) / A) * stats::pnorm(n2(v_los, s_los)) +
      ((t * s_los) / A) * (stats::dnorm(n1(v_los, s_los)) -
                             stats::dnorm(n2(v_los, s_los)))
    # the conditioning: a trial on which both drifts are negative is not a trial
    pneg <- stats::pnorm(-driftzero / sigmazero) *
      stats::pnorm(-driftone / sigmaone)
    f * (1 - cdf) / (1 - pneg)
  } else {
    0
  }
  poutlier * d_out + (1 - poutlier) * d_dec
}

make_prep <- function(y, dec, driftzero, driftone, sigmazero, sigmaone,
                      sigmabias, boundary, ndt, poutlier, n_draws = 10,
                      family = cogmod_lba2()) {
  structure(
    list(
      data = list(Y = y, dec = dec),
      family = family,
      dpars = list(
        mu = rep(driftzero, n_draws),
        driftone = rep(driftone, n_draws),
        sigmazero = rep(sigmazero, n_draws),
        sigmaone = rep(sigmaone, n_draws),
        sigmabias = rep(sigmabias, n_draws),
        boundary = rep(boundary, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}


# dcogmod_lba2 ------------------------------------------------------------

test_that("dcogmod_lba2 matches the mixture density", {
  pars <- list(driftzero = 3, driftone = 2, sigmazero = 1, sigmaone = 1,
               sigmabias = 0.5, boundary = 0.5, ndt = 0.2)
  for (poutlier in c(0, 0.02, 0.4)) {
    for (response in 0:1) {
      for (y in c(0.05, 0.15, 0.25, 0.6, 1.5, 4)) {
        expect_equal(
          do.call(dcogmod_lba2, c(list(x = y), pars,
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
  # Two things this catches. The 1/K on the outlier component - without it the
  # total comes to 1 + poutlier. And the conditioning on at least one positive
  # drift - without it the total comes to the probability of that event, which
  # at the low-drift settings below is 0.83.
  grid <- list(
    c(driftzero = 3, driftone = 2, sigmazero = 1, sigmaone = 1,
      sigmabias = 0.5, boundary = 0.5, ndt = 0.2),
    c(driftzero = 0.5, driftone = 0.2, sigmazero = 1.5, sigmaone = 1.5,
      sigmabias = 0.5, boundary = 0.5, ndt = 0.25),
    c(driftzero = -0.5, driftone = 2, sigmazero = 0.8, sigmaone = 1.2,
      sigmabias = 0.1, boundary = 1.5, ndt = 0),
    c(driftzero = 5, driftone = 5, sigmazero = 2, sigmaone = 0.5,
      sigmabias = 1.5, boundary = 0.3, ndt = 0.4)
  )
  for (g in grid) {
    for (poutlier in c(0, 0.03, 0.4)) {
      total <- sum(vapply(0:1, function(k) {
        stats::integrate(
          function(t) {
            dcogmod_lba2(t, g[["driftzero"]], g[["driftone"]],
                         g[["sigmazero"]], g[["sigmaone"]], g[["sigmabias"]],
                         g[["boundary"]], g[["ndt"]], response = k,
                         poutlier = poutlier)
          },
          lower = 0, upper = Inf, subdivisions = 4000
        )$value
      }, numeric(1)))
      expect_equal(total, 1, tolerance = 1e-4,
                   info = sprintf("driftzero = %.1f, poutlier = %.2f",
                                  g[["driftzero"]], poutlier))
    }
  }
})


test_that("the density integrates to one with both drifts negative", {
  # The case the missing normalisation was worst for: here the conditioning
  # event has probability 0.26, so an unnormalised density would come to that.
  total <- sum(vapply(0:1, function(k) {
    stats::integrate(
      function(t) {
        dcogmod_lba2(t, -2, -1.5, 1, 1, 0.5, 0.5, 0.2, response = k,
                     poutlier = 0)
      },
      lower = 0, upper = Inf, subdivisions = 4000
    )$value
  }, numeric(1)))
  expect_equal(total, 1, tolerance = 1e-4)
})


test_that("the density still integrates to one as the start-point range shrinks", {
  # The bug cogmod_lba1() was fixed for and this family shared: the density is
  # `drift * (Phi(z2) - Phi(z1)) + sigma * (phi(z1) - phi(z2))` divided by the
  # start-point range, and both differences vanish linearly in that range. The
  # survival cancels the same way. Both now switch to a Taylor branch first.
  for (A in c(0.5, 1e-2, 1e-4, 1e-6, 1e-8)) {
    total <- sum(vapply(0:1, function(k) {
      stats::integrate(
        function(t) {
          dcogmod_lba2(t, 3, 2, 1, 1, sigmabias = A, boundary = 0.5,
                       ndt = 0.2, response = k, poutlier = 0.02)
        },
        lower = 0, upper = Inf, subdivisions = 4000
      )$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-4, info = sprintf("sigmabias = %g", A))
  }
})


test_that("responses faster than ndt keep positive density (the whole point)", {
  y <- c(0.001, 0.05, 0.1, 0.19)
  for (poutlier in c(0.02, 0.4)) {
    for (response in 0:1) {
      d <- dcogmod_lba2(y, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2,
                        response = response, poutlier = poutlier)
      expect_true(all(d > 0))
      # exactly the outlier component, thinned by the two response options
      expect_equal(d, poutlier * 0.5 * 2 * stats::dnorm(y, 0, 0.2),
                   tolerance = 1e-12)
    }
  }
})


test_that("poutlier = 0 recovers the plain shifted race", {
  y <- c(0.05, 0.15, 0.2)
  expect_equal(dcogmod_lba2(y, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2, response = 0,
                            poutlier = 0),
               rep(0, length(y)))
  expect_equal(dcogmod_lba2(y, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2, response = 0,
                            poutlier = 0, log = TRUE),
               rep(-Inf, length(y)))
})


test_that("the fabricated .Machine$double.eps floor is gone", {
  # The old implementation floored the density at double.eps, so every RT below
  # the point where the density becomes representable came back with a
  # log-density of exactly log(double.eps) = -36.04: a constant the model never
  # produced, with a gradient of zero, spread over a whole stretch of the line.
  dt <- c(1e-4, 1e-2, 0.05, 0.1, 0.15, 0.2)
  lg <- dcogmod_lba2(0.2 + dt, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2, response = 0,
                     poutlier = 0, log = TRUE)
  expect_false(any(abs(lg - log(.Machine$double.eps)) < 1e-6))

  # Very close to the shift the LBA density really has underflowed - at
  # dt = 0.01 the standardized threshold is 47, so the true density is around
  # 1e-480 - and -Inf is the honest answer there, not a floor.
  expect_true(all(is.infinite(lg[dt <= 1e-2])))
  # Beyond that it is representable, and rises away from the shift.
  live <- lg[dt >= 0.05]
  expect_true(all(is.finite(live)))
  expect_true(all(diff(live) > 0))

  # And the mixture is finite throughout, which is what the outlier component
  # is for: the sampler never sees a -Inf log-likelihood.
  mix <- dcogmod_lba2(0.2 + dt, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2, response = 0,
                      poutlier = 0.02, log = TRUE)
  expect_true(all(is.finite(mix)))
})


test_that("dcogmod_lba2 is vectorized over every argument", {
  n <- 7
  d <- dcogmod_lba2(seq(0.3, 1.2, length.out = n),
                    driftzero = seq(1, 3, length.out = n),
                    driftone = seq(3, 1, length.out = n),
                    sigmazero = seq(0.5, 1.5, length.out = n),
                    sigmaone = seq(1.5, 0.5, length.out = n),
                    sigmabias = seq(0.1, 0.8, length.out = n),
                    boundary = seq(0.2, 0.9, length.out = n),
                    ndt = seq(0.05, 0.25, length.out = n),
                    response = rep(0:1, length.out = n),
                    poutlier = seq(0, 0.1, length.out = n))
  expect_length(d, n)
  expect_true(all(is.finite(d) & d >= 0))
})


test_that("dcogmod_lba2 returns 0 density for invalid parameters", {
  args <- list(x = 0.5, driftzero = 3, driftone = 2, sigmazero = 1,
               sigmaone = 1, sigmabias = 0.5, boundary = 0.5, ndt = 0.2,
               response = 0)
  for (bad in list(list(sigmazero = 0), list(sigmaone = -1),
                   list(sigmabias = -1), list(boundary = 0), list(ndt = -0.1),
                   list(poutlier = 1.5), list(response = 2))) {
    expect_warning(d <- do.call(dcogmod_lba2, modifyList(args, bad)),
                   info = names(bad))
    expect_equal(d, 0)
  }
  # a negative mean drift is not invalid: it just makes that accumulator
  # unlikely to respond
  expect_silent(d <- dcogmod_lba2(0.5, driftzero = -1, driftone = 2,
                                  response = 0))
  expect_true(is.finite(d) && d > 0)

  # `sigmabias = 0` is not invalid either: it is the nested no-start-point-
  # variability race, and the density's Taylor branch evaluates it exactly.
  expect_silent(d <- do.call(dcogmod_lba2, modifyList(args, list(sigmabias = 0))))
  expect_true(is.finite(d) && d > 0)
})


# rcogmod_lba2 ------------------------------------------------------------

test_that("rcogmod_lba2 returns rt and response", {
  set.seed(1)
  sim <- rcogmod_lba2(200, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2)
  expect_true(is.data.frame(sim))
  expect_named(sim, c("rt", "response"))
  expect_equal(nrow(sim), 200)
  expect_true(all(sim$response %in% c(0, 1)))
  expect_true(all(is.finite(sim$rt)))
  expect_true(all(sim$rt > 0.2))
  # driftzero > driftone, so accumulator 0 wins more often
  expect_gt(mean(sim$response == 0), 0.5)
})


test_that("rcogmod_lba2 reproduces its own density", {
  # The check the conditioning exists for. Before it, the empirical choice
  # proportion at the low-drift settings was 0.562 against an integral of 0.468:
  # the sampler and the density described different models.
  set.seed(42)
  n <- 40000
  for (g in list(
    list(v0 = 3, v1 = 2, s0 = 1, s1 = 1, A = 0.5, b = 0.5, ndt = 0.2),
    list(v0 = 0.5, v1 = 0.2, s0 = 1.5, s1 = 1.5, A = 0.5, b = 0.5, ndt = 0.25)
  )) {
    for (poutlier in c(0, 0.05)) {
      sim <- rcogmod_lba2(n, g$v0, g$v1, g$s0, g$s1, g$A, g$b, g$ndt,
                          poutlier = poutlier)
      for (k in 0:1) {
        integrated <- stats::integrate(
          function(t) {
            dcogmod_lba2(t, g$v0, g$v1, g$s0, g$s1, g$A, g$b, g$ndt,
                         response = k, poutlier = poutlier)
          },
          lower = 0, upper = Inf, subdivisions = 4000
        )$value
        expect_equal(mean(sim$response == k), integrated, tolerance = 0.015,
                     info = sprintf("v0 %.1f, response %d, poutlier %.2f",
                                    g$v0, k, poutlier))
      }
      for (pr in c(0.25, 0.5, 0.9)) {
        qq <- unname(stats::quantile(sim$rt, pr))
        cdf <- sum(vapply(0:1, function(k) {
          stats::integrate(
            function(t) {
              dcogmod_lba2(t, g$v0, g$v1, g$s0, g$s1, g$A, g$b, g$ndt,
                           response = k, poutlier = poutlier)
            },
            lower = 0, upper = qq, subdivisions = 4000
          )$value
        }, numeric(1)))
        expect_equal(cdf, pr, tolerance = 0.015,
                     info = sprintf("v0 %.1f, p %.2f, poutlier %.2f",
                                    g$v0, pr, poutlier))
      }
    }
  }
})


test_that("rcogmod_lba2 conditions exactly, with both drifts negative", {
  # Both means below zero, so most draws from the unconditional joint would
  # produce no response at all. The exact split - accumulator 0 positive, or
  # accumulator 0 non-positive and accumulator 1 positive - replaces the
  # rejection loop the old implementation bounded with `max_iter` and then gave
  # up on, forcing a drift positive with `abs()` and silently drawing from the
  # wrong distribution.
  set.seed(5)
  sim <- rcogmod_lba2(40000, -2, -1.5, 1, 1, 0.5, 0.5, ndt = 0.2)
  expect_true(all(is.finite(sim$rt)))
  p0 <- stats::integrate(
    function(t) dcogmod_lba2(t, -2, -1.5, 1, 1, 0.5, 0.5, 0.2, response = 0),
    lower = 0, upper = Inf, subdivisions = 4000
  )$value
  expect_equal(mean(sim$response == 0), p0, tolerance = 0.015)
})


test_that("rcogmod_lba2 produces responses below ndt only via outliers", {
  set.seed(2)
  expect_true(all(rcogmod_lba2(2000, ndt = 0.4, poutlier = 0)$rt > 0.4))
  sim <- rcogmod_lba2(20000, ndt = 0.4, poutlier = 0.3)
  expect_true(any(sim$rt < 0.4))
  fast <- sim$response[sim$rt < 0.02]
  expect_setequal(unique(fast), c(0, 1))
  expect_equal(mean(fast == 0), 0.5, tolerance = 0.1)
})


test_that("rcogmod_lba2 errors on invalid parameters", {
  expect_error(rcogmod_lba2(10, sigmazero = 0), "sigmazero")
  expect_error(rcogmod_lba2(10, sigmabias = -1), "sigmabias")
  expect_error(rcogmod_lba2(10, boundary = 0), "boundary")
  expect_error(rcogmod_lba2(10, ndt = -1), "ndt")
  expect_error(rcogmod_lba2(10, poutlier = 2), "poutlier")

  # zero is the closed end of `sigmabias`'s support, not a violation of it
  expect_silent(rcogmod_lba2(10, sigmabias = 0))
})


# the shared LBA kernels --------------------------------------------------

test_that("the shared survival is exact where 1 - CDF would cancel", {
  # S(t) = (g(z2) - g(z1)) / delta, taken directly rather than as 1 - CDF.
  # Checked against integrating the defective density it complements.
  for (prm in list(c(v = 3, s = 1, A = 0.5, k = 0.5),
                   c(v = 1, s = 1.5, A = 0.05, k = 1),
                   c(v = 0.5, s = 2, A = 1e-5, k = 0.4))) {
    v <- prm[["v"]]; s <- prm[["s"]]; A <- prm[["A"]]; k <- prm[["k"]]
    for (t in c(0.1, 0.5, 2)) {
      st <- max(s * t, 1e-10)
      surv <- cogmod:::.lba_surv_raw(v, s, st, (k - v * t) / st, A / st)
      dens <- function(u) {
        stu <- pmax(s * u, 1e-10)
        cogmod:::.lba_dens_over_A(v, s, stu, (k - v * u) / stu, A / stu)
      }
      num <- 1 - stats::integrate(dens, 0, t, subdivisions = 3000,
                                  rel.tol = 1e-11)$value
      expect_equal(surv, num, tolerance = 1e-9,
                   label = sprintf("S at v=%g A=%g t=%g", v, A, t))
    }
  }
})


test_that("the shared kernels are the ones cogmod_lba1() uses", {
  # cogmod_lba1() and cogmod_lba2() differ in how many accumulators they run,
  # not in the per-accumulator arithmetic, so they share one implementation.
  expect_true(is.function(cogmod:::.lba_dens_over_A))
  expect_true(is.function(cogmod:::.lba_surv_raw))
  expect_true(grepl("cogmod_lba_dens_over_A", cogmod:::.LBA_STAN_PRELUDE))
  expect_true(grepl("cogmod_lba_surv", cogmod:::.LBA_STAN_PRELUDE))
  # ... and cogmod_lba2()'s prelude is that one plus its own race on top
  expect_true(startsWith(cogmod:::.LBA2_STAN_PRELUDE,
                         cogmod:::.LBA_STAN_PRELUDE))
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_lba2 matches dcogmod_lba2", {
  prep <- make_prep(y = c(0.5, 0.15), dec = c(1, 0), driftzero = 3,
                    driftone = 2, sigmazero = 1, sigmaone = 1, sigmabias = 0.5,
                    boundary = 0.5, ndt = 0.2, poutlier = 0.03)
  for (i in 1:2) {
    expect_equal(
      log_lik_cogmod_lba2(i, prep),
      rep(dcogmod_lba2(prep$data$Y[i], 3, 2, 1, 1, 0.5, 0.5, 0.2,
                       response = prep$data$dec[i], poutlier = 0.03,
                       log = TRUE), 10)
    )
  }
})


test_that("log_lik_cogmod_lba2 uses the observed choice", {
  p0 <- make_prep(y = 0.5, dec = 0, driftzero = 4, driftone = 0.5,
                  sigmazero = 0.5, sigmaone = 0.5, sigmabias = 0.5,
                  boundary = 0.5, ndt = 0.2, poutlier = 0)
  p1 <- make_prep(y = 0.5, dec = 1, driftzero = 4, driftone = 0.5,
                  sigmazero = 0.5, sigmaone = 0.5, sigmabias = 0.5,
                  boundary = 0.5, ndt = 0.2, poutlier = 0)
  expect_false(isTRUE(all.equal(log_lik_cogmod_lba2(1, p0),
                                log_lik_cogmod_lba2(1, p1))))
})


test_that("log_lik_cogmod_lba2 errors when dec is missing", {
  prep <- make_prep(y = 0.5, dec = 0, driftzero = 3, driftone = 2,
                    sigmazero = 1, sigmaone = 1, sigmabias = 0.5,
                    boundary = 0.5, ndt = 0.2, poutlier = 0)
  prep$data$dec <- NULL
  expect_error(log_lik_cogmod_lba2(1, prep), "dec")
})


test_that("posterior_predict_cogmod_lba2 excludes outliers by default", {
  set.seed(3)
  prep <- make_prep(y = 0.5, dec = 0, driftzero = 3, driftone = 2,
                    sigmazero = 1, sigmaone = 1, sigmabias = 0.5,
                    boundary = 0.5, ndt = 0.4, poutlier = 0.9, n_draws = 500)
  out <- posterior_predict_cogmod_lba2(1, prep)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(500L, 2L))
  expect_true(all(out[, 1] > 0.4))
  expect_true(all(out[, 2] %in% c(0, 1)))

  with_out <- posterior_predict_cogmod_lba2(1, prep, predict_outliers = TRUE)
  expect_true(any(with_out[, 1] < 0.4))
})


test_that("the family flag drives predictions when no argument is given", {
  set.seed(4)
  prep <- make_prep(y = 0.5, dec = 0, driftzero = 3, driftone = 2,
                    sigmazero = 1, sigmaone = 1, sigmabias = 0.5,
                    boundary = 0.5, ndt = 0.4, poutlier = 0.9, n_draws = 500,
                    family = cogmod_lba2(predict_outliers = TRUE))
  expect_true(any(posterior_predict_cogmod_lba2(1, prep)[, 1] < 0.4))
  expect_true(all(
    posterior_predict_cogmod_lba2(1, prep, predict_outliers = FALSE)[, 1] > 0.4
  ))
})


test_that("posterior_epred_cogmod_lba2 refuses rather than guessing", {
  expect_error(posterior_epred_cogmod_lba2(NULL), "posterior_predict")
})


# family ------------------------------------------------------------------

test_that("cogmod_lba2() builds a valid brms custom family", {
  fam <- cogmod_lba2()
  expect_s3_class(fam, "customfamily")
  # Order matters: brms passes the dpars to cogmod_lba2_lpdf in exactly this one.
  expect_equal(fam$dpars, c("mu", "driftone", "sigmazero", "sigmaone",
                            "sigmabias", "boundary", "ndt", "poutlier"))
  expect_equal(unname(cogmod:::.family_links(fam)),
               c("identity", "identity", "softplus", "softplus", "softplus",
                 "softplus", "log", "logit"))
  expect_equal(fam$vars, "dec[n]")
  expect_false(fam$predict_outliers)
  expect_false(any(c("tau", "minrt") %in% fam$dpars))
})


test_that("cogmod_lba2 is one of the outlier-mixture families", {
  expect_true("cogmod_lba2" %in% cogmod:::.OUTLIER_FAMILIES)
  expect_true("cogmod_lba2" %in% cogmod:::.cogmod_families())
  expect_length(cogmod:::.cogmod_families(),
                length(unique(cogmod:::.cogmod_families())))
})


test_that("predict_outliers survives the validation brm() performs", {
  fam <- brms:::validate_family(cogmod_lba2(predict_outliers = TRUE))
  expect_true(fam$predict_outliers)
})


# Stan code ---------------------------------------------------------------

test_that("stanvars carry the likelihood with the outlier component", {
  code <- cogmod_lba2_stanvars()[[1]]$scode
  expect_true(grepl("real cogmod_lba2_lpdf", code))
  expect_true(grepl("int dec", code))
  expect_true(grepl("log_mix\\(poutlier", code))
  expect_true(grepl("half Normal with scale 0.2", code))
  # the shared kernels and the race built on them
  expect_true(grepl("real cogmod_lba_dens_over_A", code))
  expect_true(grepl("real cogmod_lba_surv", code))
  expect_true(grepl("real cogmod_lba2_decision_lpdf", code))
  # the conditioning
  expect_true(grepl("log1m_exp\\(lq\\)", code))
  # the old parameterization is gone
  expect_false(grepl("real tau", code))
  expect_false(grepl("real minrt", code))
})


test_that("Stan cogmod_lba2_lpdf matches dcogmod_lba2", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_lba2")
  grid <- covering_grid(
    Y = c(0.02, 0.25, 0.5, 1.2, 4),
    mu = c(-0.5, 0.5, 3),
    driftone = c(0.2, 2),
    sigmazero = c(0.3, 1, 2),
    sigmaone = c(0.5, 1.5),
    sigmabias = c(0.02, 0.5, 1.5),
    boundary = c(0.2, 1.0),
    ndt = c(0, 0.15, 0.4),
    poutlier = c(0, 0.001, 0.4),
    dec = 0:1,
    # a negative mean drift makes the conditioning bite, so that slice is swept
    # in full rather than left to the covering subset
    always = function(g) {
      g$mu == -0.5 & g$sigmazero == 1 & g$boundary == 0.2 & g$poutlier == 0.001
    }
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    stan <- lpdf(g$Y, g$mu, g$driftone, g$sigmazero, g$sigmaone, g$sigmabias,
                 g$boundary, g$ndt, g$poutlier, as.integer(g$dec))
    r <- dcogmod_lba2(g$Y, g$mu, g$driftone, g$sigmazero, g$sigmaone,
                      g$sigmabias, g$boundary, g$ndt, response = g$dec,
                      poutlier = g$poutlier, log = TRUE)
    if (is.infinite(r)) {
      expect_equal(stan, r)
    } else {
      # Relative, not absolute: these log-densities reach large magnitudes in
      # the far tails, where one ulp is already ~1e-7 absolute.
      expect_lt(abs(stan - r) / max(1, abs(r)), 1e-9)
    }
  }

  # invalid arguments are rejected on both sides
  expect_equal(lpdf(0.5, 3, 2, 0, 1, 0.5, 0.5, 0.2, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 0, 0.5, 0.5, 0.2, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 1, -1, 0.5, 0.2, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 1, 0.5, 0, 0.2, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 1, 0.5, 0.5, -0.1, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 1, 0.5, 0.5, 0.2, 1.5, 0L), -Inf)
  expect_equal(lpdf(0.5, 3, 2, 1, 1, 0.5, 0.5, 0.2, 0.02, 2L), -Inf)
  expect_equal(lpdf(0, 3, 2, 1, 1, 0.5, 0.5, 0.2, 0.02, 0L), -Inf)

  # `sigmabias = 0` is the nested no-start-point-variability race, not an
  # invalid argument, and the two implementations have to agree there too -
  # that shared likelihood is what makes loo_compare() against the general
  # model like-for-like rather than a comparison of two different densities.
  for (resp in 0:1) {
    s0 <- lpdf(0.5, 3, 2, 1, 1, 0, 0.5, 0.2, 0.02, as.integer(resp))
    r0 <- dcogmod_lba2(0.5, 3, 2, 1, 1, sigmabias = 0, boundary = 0.5,
                       ndt = 0.2, response = resp, poutlier = 0.02, log = TRUE)
    expect_true(is.finite(s0))
    expect_lt(abs(s0 - r0) / max(1, abs(r0)), 1e-9)
  }
})


# priors and inits --------------------------------------------------------

test_that("cogmod_priors fills ndt and poutlier for cogmod_lba2", {
  set.seed(6)
  sim <- rcogmod_lba2(150, ndt = 0.2, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150),
                  id = rep(1:10, length.out = 150))

  modelled <- brms::bf(RT | dec(Error) ~ 1, driftone ~ 1, sigmazero ~ 1,
                       sigmaone ~ 1, sigmabias ~ 1, boundary ~ 1, ndt ~ 1,
                       poutlier ~ 1, family = cogmod_lba2())
  p <- cogmod_priors(modelled, d)
  expect_true(any(p$dpar == "ndt" & p$class == "Intercept" &
                    p$prior == "normal(-1.2, 0.2)"))
  expect_true(any(p$dpar == "poutlier" & p$class == "Intercept" &
                    p$prior == "normal(-5, 1)"))

  omitted <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_lba2())
  expect_true(any(grepl("uniform",
                        brms::get_prior(omitted, data = d,
                                        family = cogmod_lba2())$prior)))
  p2 <- cogmod_priors(omitted, d)
  expect_false(any(grepl("uniform", p2$prior)))
  expect_true(any(p2$class == "ndt" & p2$prior == "lognormal(-1.2, 0.2)"))
  expect_true(any(p2$class == "poutlier" & p2$prior == "exponential(100)"))

  mixed <- brms::bf(RT | dec(Error) ~ Condition + (1 | id), ndt ~ Condition,
                    sigmazero = 1, family = cogmod_lba2())
  p3 <- cogmod_priors(mixed, d)
  expect_s3_class(p3, "brmsprior")
  # sigmazero is fixed in the formula, so it is not a parameter at all
  expect_false(any(p3$dpar == "sigmazero" | p3$class == "sigmazero"))
  code <- brms::make_stancode(mixed, data = d, prior = p3,
                              stanvars = cogmod_stanvars(mixed))
  expect_false(grepl("uniform\\(", code))
  expect_true(grepl("cogmod_lba2_lpdf", code))
})


test_that("cogmod_priors fences off the scale ray and the threshold ridge", {
  # Multiply the drifts, their SDs, the start-point range and the threshold all
  # by c > 0 and every finishing time (b - z) / v is unchanged, so the
  # likelihood is *exactly* constant along that ray - and it runs to infinity.
  set.seed(11)
  sim <- rcogmod_lba2(600, 3, 2, 1, 1, 0.5, 0.5, ndt = 0.2, poutlier = 0.02)
  ll <- function(cc) {
    sum(dcogmod_lba2(sim$rt, 3 * cc, 2 * cc, 1 * cc, 1 * cc, 0.5 * cc,
                     0.5 * cc, ndt = 0.2, response = sim$response,
                     poutlier = 0.02, log = TRUE))
  }
  expect_equal(ll(1), ll(3), tolerance = 1e-8)
  expect_equal(ll(1), ll(0.1), tolerance = 1e-8)

  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = nrow(sim)))
  f <- brms::bf(RT | dec(Error) ~ Condition, driftone ~ Condition,
                sigmazero ~ Condition, sigmaone ~ Condition,
                sigmabias ~ Condition, boundary ~ Condition, ndt ~ 1,
                family = cogmod_lba2())

  # brms leaves all four flat on its own
  raw <- brms::get_prior(f, data = d, family = cogmod_lba2())
  flat <- raw$dpar %in% c("sigmazero", "sigmaone", "sigmabias", "boundary")
  expect_true(all(raw$prior[flat] == ""))

  p <- cogmod_priors(f, d)
  pick <- function(dp, cls) {
    p$prior[p$dpar == dp & p$class == cls & nzchar(p$prior)]
  }
  for (dp in c("sigmazero", "sigmaone", "sigmabias", "boundary")) {
    expect_equal(pick(dp, "Intercept"), "normal(0, 1)", label = dp)
    expect_equal(pick(dp, "b"), "normal(0, 0.5)", label = paste(dp, "slope"))
  }

  f2 <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_lba2())
  p2 <- cogmod_priors(f2, d)
  for (dp in c("sigmazero", "sigmaone", "sigmabias", "boundary")) {
    expect_equal(p2$prior[p2$class == dp], "lognormal(-0.7, 0.75)", label = dp)
  }

  for (form in list(f, f2)) {
    # Both leave the evidence scale free, which is the whole point of this
    # test, so cogmod_stanvars() warns about it - see the test below. What is
    # being checked here is that the generated Stan itself is clean.
    sv <- suppressWarnings(cogmod_stanvars(form))
    expect_silent(
      brms::make_stancode(form, data = d, family = cogmod_lba2(),
                          prior = cogmod_priors(form, d), stanvars = sv)
    )
  }
})


test_that("cogmod_stanvars warns when the evidence scale is left free", {
  free <- brms::bf(RT | dec(Error) ~ 1, driftone ~ 1, sigmazero ~ 1,
                   sigmaone ~ 1, sigmabias ~ 1, boundary ~ 1, ndt ~ 1,
                   poutlier ~ 1, family = cogmod_lba2())
  expect_warning(cogmod_stanvars(free), "evidence scale is arbitrary")
  expect_warning(cogmod_stanvars(free), "sigmazero = 1")

  # Omitting a dpar from bf() does not fix it - brms estimates it anyway - so
  # the ray is exactly as free as before. This is the trap the warning is for.
  omitted <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_lba2())
  expect_warning(cogmod_stanvars(omitted), "evidence scale is arbitrary")

  # cogmod_lba1() has the same ray, and names its own convention.
  expect_warning(
    cogmod_stanvars(brms::bf(RT ~ 1, sigma ~ 1, sigmabias ~ 1, boundary ~ 1,
                             family = cogmod_lba1())),
    "sigma = 1"
  )

  # Any single member of the ray pins it, not just the conventional one.
  for (pinned in list(
    brms::bf(RT | dec(Error) ~ 1, sigmazero = 1, family = cogmod_lba2()),
    brms::bf(RT | dec(Error) ~ 1, sigmaone = 1, family = cogmod_lba2()),
    brms::bf(RT | dec(Error) ~ 1, boundary = 1, family = cogmod_lba2()),
    brms::bf(RT | dec(Error) ~ 1, sigmabias = 0.5, family = cogmod_lba2())
  )) {
    expect_no_warning(cogmod_stanvars(pinned))
  }

  # Except at ZERO. Scaling the ray multiplies every member by c, and c * 0 = 0,
  # so `sigmabias = 0` - the recinormal, a legitimate model since the bound was
  # closed - takes that parameter off the ray instead of pinning the ray. The
  # remaining members are as free as before, and the warning has to survive.
  for (zero in list(
    brms::bf(RT | dec(Error) ~ 1, sigmabias = 0, family = cogmod_lba2()),
    brms::bf(RT ~ 1, sigmabias = 0, family = cogmod_lba1())
  )) {
    expect_warning(cogmod_stanvars(zero), "evidence scale is arbitrary")
    expect_warning(cogmod_stanvars(zero), "does not count")
  }

  # A second, non-zero pin alongside it does the job.
  expect_no_warning(cogmod_stanvars(
    brms::bf(RT ~ 1, sigmabias = 0, boundary = 1, family = cogmod_lba1())))
  expect_no_warning(cogmod_stanvars(
    brms::bf(RT | dec(Error) ~ 1, sigmabias = 0, sigmazero = 1,
             family = cogmod_lba2())))

  # The families whose scale is pinned by construction stay quiet: the RDM and
  # the DDM both have a unit diffusion coefficient baked into the likelihood.
  expect_no_warning(cogmod_stanvars(
    brms::bf(RT | dec(Error) ~ 1, family = cogmod_rdm())))
  expect_no_warning(cogmod_stanvars(
    brms::bf(RT | dec(Error) ~ 1, family = cogmod_ddm())))
  expect_no_warning(cogmod_stanvars(
    brms::bf(RT ~ 1, family = cogmod_lognormal())))

  # Nothing to read a formula off, so nothing to say.
  expect_no_warning(cogmod_stanvars(cogmod_lba2()))
})


test_that("cogmod_inits covers the declared parameters", {
  set.seed(8)
  sim <- rcogmod_lba2(150, ndt = 0.2, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150))
  f <- brms::bf(RT | dec(Error) ~ Condition, driftone ~ 1, sigmazero ~ 1,
                sigmaone ~ 1, sigmabias ~ 1, boundary ~ 1, ndt ~ 1,
                poutlier ~ 1, family = cogmod_lba2())
  vals <- cogmod_inits(f, d)(1)
  code <- suppressWarnings(
    brms::make_stancode(f, data = d, family = cogmod_lba2())
  )
  declared <- vapply(cogmod:::.stan_param_decls(code), `[[`, character(1),
                     "name")
  expect_setequal(names(vals), declared)
  expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))))

  v0 <- cogmod_inits(f, d, jitter = 0)(1)
  expect_lt(exp(v0$Intercept_ndt), 0.3)
  expect_equal(log1p(exp(v0$Intercept_sigmazero)), 1, tolerance = 1e-8)
  expect_equal(log1p(exp(v0$Intercept_sigmabias)), 0.5, tolerance = 1e-8)
  expect_equal(v0$Intercept_driftone, 3, tolerance = 1e-8)
})


# a real fit --------------------------------------------------------------

test_that("cogmod_lba2 recovers ndt above the fastest observed response", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n <- 2000
  sim <- rcogmod_lba2(n, driftzero = 3, driftone = 2, sigmazero = 1,
                      sigmaone = 1, sigmabias = 0.5, boundary = 0.5,
                      ndt = 0.25, poutlier = 0.04)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  expect_lt(min(d$RT), 0.25)

  # sigmazero fixed to 1: the evidence scale is otherwise unidentified
  f <- brms::bf(RT | dec(Error) ~ 1, driftone ~ 1, sigmazero = 1,
                sigmaone ~ 1, sigmabias ~ 1, boundary ~ 1, ndt ~ 1,
                poutlier ~ 1, family = cogmod_lba2())
  fit <- brms::brm(f, data = d, prior = cogmod_priors(f, d),
                   init = cogmod_inits(f, d), stanvars = cogmod_stanvars(f),
                   algorithm = "pathfinder", refresh = 0,
                   backend = "cmdstanr")

  sp <- function(x) log1p(exp(x))
  s <- brms::posterior_summary(fit)[, "Estimate"]
  expect_equal(s[["b_Intercept"]], 3, tolerance = 0.4)
  expect_equal(s[["b_driftone_Intercept"]], 2, tolerance = 0.4)
  expect_equal(sp(s[["b_sigmaone_Intercept"]]), 1, tolerance = 0.3)
  # sigmabias and boundary trade off along a ridge; their sum does not
  expect_equal(sp(s[["b_boundary_Intercept"]]) +
                 sp(s[["b_sigmabias_Intercept"]]), 1, tolerance = 0.2)
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
