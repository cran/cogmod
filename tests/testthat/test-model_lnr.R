context("LNR - shifted Log-Normal Race with outliers")

# Reference implementation of the mixture, written out longhand. The 1/K on the
# outlier component is what keeps the density summing to one over the response
# options; without it the total comes to 1 + poutlier.
ref_dens <- function(y, nuzero, nuone, sigmazero, sigmaone, ndt, response,
                     poutlier) {
  d_out <- 2 * stats::dnorm(y, 0, 0.2) / 2
  ml_win <- if (response == 0) -nuzero else -nuone
  sd_win <- if (response == 0) sigmazero else sigmaone
  ml_los <- if (response == 0) -nuone else -nuzero
  sd_los <- if (response == 0) sigmaone else sigmazero
  t <- y - ndt
  d_dec <- if (t > 0) {
    stats::dlnorm(t, ml_win, sd_win) *
      stats::plnorm(t, ml_los, sd_los, lower.tail = FALSE)
  } else {
    0
  }
  poutlier * d_out + (1 - poutlier) * d_dec
}

make_prep <- function(y, dec, nuzero, nuone, sigmazero, sigmaone, ndt,
                      poutlier, n_draws = 10, family = cogmod_lnr()) {
  structure(
    list(
      data = list(Y = y, dec = dec),
      family = family,
      dpars = list(
        mu = rep(nuzero, n_draws),
        nuone = rep(nuone, n_draws),
        sigmazero = rep(sigmazero, n_draws),
        sigmaone = rep(sigmaone, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}


# dcogmod_lnr -------------------------------------------------------------

test_that("dcogmod_lnr matches the mixture density", {
  pars <- list(nuzero = 0.5, nuone = 0.2, sigmazero = 0.8, sigmaone = 1.0,
               ndt = 0.25)
  for (poutlier in c(0, 0.02, 0.4)) {
    for (response in 0:1) {
      for (y in c(0.05, 0.2, 0.3, 0.6, 1.5, 4)) {
        expect_equal(
          do.call(dcogmod_lnr, c(list(x = y), pars,
                                 list(response = response,
                                      poutlier = poutlier))),
          do.call(ref_dens, c(list(y = y), pars,
                              list(response = response,
                                   poutlier = poutlier))),
          tolerance = 1e-12,
          info = sprintf("y = %.2f, response = %d, poutlier = %.2f",
                         y, response, poutlier)
        )
      }
    }
  }
})


test_that("the density sums to one over responses and integrates to one", {
  # The check the 1/K exists for: without it the total comes to 1 + poutlier,
  # and every individual density still looks perfectly reasonable.
  grid <- list(
    c(nuzero = 0.5, nuone = 0.2, sigmazero = 0.8, sigmaone = 1.0, ndt = 0.25),
    c(nuzero = -0.5, nuone = 0.8, sigmazero = 0.5, sigmaone = 0.4, ndt = 0.10),
    c(nuzero = 1.2, nuone = 1.2, sigmazero = 1.2, sigmaone = 0.6, ndt = 0.40)
  )
  for (g in grid) {
    for (poutlier in c(0, 0.03, 0.4)) {
      total <- sum(vapply(0:1, function(k) {
        stats::integrate(
          function(t) {
            dcogmod_lnr(t, g[["nuzero"]], g[["nuone"]], g[["sigmazero"]],
                        g[["sigmaone"]], g[["ndt"]], response = k,
                        poutlier = poutlier)
          },
          lower = 0, upper = Inf, subdivisions = 2000
        )$value
      }, numeric(1)))
      expect_equal(total, 1, tolerance = 1e-4,
                   info = sprintf("ndt = %.2f, poutlier = %.2f",
                                  g[["ndt"]], poutlier))
    }
  }
})


test_that("responses faster than ndt keep positive density (the whole point)", {
  y <- c(0.001, 0.05, 0.1, 0.19)
  for (poutlier in c(0.02, 0.4)) {
    for (response in 0:1) {
      d <- dcogmod_lnr(y, 0.5, 0.2, 0.8, 1.0, ndt = 0.2, response = response,
                       poutlier = poutlier)
      expect_true(all(d > 0))
      # exactly the outlier component, thinned by the two response options
      expect_equal(d, poutlier * 0.5 * 2 * stats::dnorm(y, 0, 0.2),
                   tolerance = 1e-12)
    }
  }
})


test_that("poutlier = 0 recovers the plain shifted race", {
  y <- c(0.05, 0.15, 0.2)
  expect_equal(dcogmod_lnr(y, 0.5, 0.2, 0.8, 1.0, ndt = 0.2, response = 0,
                           poutlier = 0),
               rep(0, length(y)))
  expect_equal(dcogmod_lnr(y, 0.5, 0.2, 0.8, 1.0, ndt = 0.2, response = 0,
                           poutlier = 0, log = TRUE),
               rep(-Inf, length(y)))

  x <- seq(0.25, 4, length.out = 50)
  for (response in 0:1) {
    expect_equal(
      dcogmod_lnr(x, 0.5, 0.2, 0.8, 1.0, ndt = 0.2, response = response),
      vapply(x, ref_dens, numeric(1), nuzero = 0.5, nuone = 0.2,
             sigmazero = 0.8, sigmaone = 1.0, ndt = 0.2, response = response,
             poutlier = 0),
      tolerance = 1e-12
    )
  }
})


test_that("dcogmod_lnr is vectorized over every argument", {
  n <- 7
  d <- dcogmod_lnr(seq(0.3, 1.2, length.out = n),
                   nuzero = seq(0, 1, length.out = n),
                   nuone = seq(1, 0, length.out = n),
                   sigmazero = seq(0.5, 1, length.out = n),
                   sigmaone = seq(1, 0.5, length.out = n),
                   ndt = seq(0.05, 0.25, length.out = n),
                   response = rep(0:1, length.out = n),
                   poutlier = seq(0, 0.1, length.out = n))
  expect_length(d, n)
  expect_true(all(is.finite(d) & d >= 0))
})


test_that("dcogmod_lnr returns 0 density for invalid parameters", {
  args <- list(x = 0.5, nuzero = 0.5, nuone = 0.2, sigmazero = 0.8,
               sigmaone = 1.0, ndt = 0.2, response = 0)
  expect_warning(d <- do.call(dcogmod_lnr, modifyList(args, list(sigmazero = -1))))
  expect_equal(d, 0)
  expect_warning(d <- do.call(dcogmod_lnr, modifyList(args, list(sigmaone = 0))))
  expect_equal(d, 0)
  expect_warning(d <- do.call(dcogmod_lnr, modifyList(args, list(ndt = -0.1))))
  expect_equal(d, 0)
  expect_warning(d <- do.call(dcogmod_lnr, modifyList(args, list(poutlier = 1.5))))
  expect_equal(d, 0)
  # `response` is checked against the K options the family declares
  expect_warning(d <- do.call(dcogmod_lnr, modifyList(args, list(response = 2))))
  expect_equal(d, 0)
})


# rcogmod_lnr -------------------------------------------------------------

test_that("rcogmod_lnr returns rt and response", {
  set.seed(1)
  sim <- rcogmod_lnr(100, 0.5, 0.2, 0.8, 1.0, ndt = 0.2)
  expect_true(is.data.frame(sim))
  expect_named(sim, c("rt", "response"))
  expect_equal(nrow(sim), 100)
  expect_true(all(sim$response %in% c(0, 1)))
  # nuzero > nuone, so accumulator 0 wins more often
  expect_gt(mean(sim$response == 0), 0.5)
})


test_that("rcogmod_lnr reproduces its own density", {
  set.seed(42)
  n <- 20000
  pars <- list(nuzero = 0.5, nuone = 0.2, sigmazero = 0.8, sigmaone = 1.0,
               ndt = 0.25)
  for (poutlier in c(0, 0.05)) {
    sim <- do.call(rcogmod_lnr, c(list(n = n), pars,
                                  list(poutlier = poutlier)))
    for (k in 0:1) {
      integrated <- stats::integrate(
        function(t) {
          do.call(dcogmod_lnr, c(list(x = t), pars,
                                 list(response = k, poutlier = poutlier)))
        },
        lower = 0, upper = Inf, subdivisions = 2000
      )$value
      expect_equal(mean(sim$response == k), integrated, tolerance = 0.01,
                   info = sprintf("response %d, poutlier %.2f", k, poutlier))
    }
    # marginal CDF at the median
    med <- stats::median(sim$rt)
    cdf <- sum(vapply(0:1, function(k) {
      stats::integrate(
        function(t) {
          do.call(dcogmod_lnr, c(list(x = t), pars,
                                 list(response = k, poutlier = poutlier)))
        },
        lower = 0, upper = med, subdivisions = 2000
      )$value
    }, numeric(1)))
    expect_equal(cdf, 0.5, tolerance = 0.01)
  }
})


test_that("rcogmod_lnr produces responses below ndt only via outliers", {
  set.seed(2)
  expect_true(all(rcogmod_lnr(2000, ndt = 0.4, poutlier = 0)$rt > 0.4))
  sim <- rcogmod_lnr(20000, ndt = 0.4, poutlier = 0.3)
  expect_true(any(sim$rt < 0.4))
  # the contaminant produces a choice as well as a time, uniform over the two
  fast <- sim$response[sim$rt < 0.02]
  expect_setequal(unique(fast), c(0, 1))
  expect_equal(mean(fast == 0), 0.5, tolerance = 0.1)
})


test_that("rcogmod_lnr errors on invalid parameters", {
  expect_error(rcogmod_lnr(10, sigmazero = -1), "sigmazero")
  expect_error(rcogmod_lnr(10, ndt = -1), "ndt")
  expect_error(rcogmod_lnr(10, poutlier = 2), "poutlier")
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_lnr matches dcogmod_lnr", {
  prep <- make_prep(y = c(0.5, 0.15), dec = c(1, 0), nuzero = 0.5, nuone = 0.2,
                    sigmazero = 0.8, sigmaone = 1.0, ndt = 0.25,
                    poutlier = 0.03)
  for (i in 1:2) {
    expect_equal(
      log_lik_cogmod_lnr(i, prep),
      rep(dcogmod_lnr(prep$data$Y[i], 0.5, 0.2, 0.8, 1.0, 0.25,
                      response = prep$data$dec[i], poutlier = 0.03,
                      log = TRUE), 10)
    )
  }
})


test_that("log_lik_cogmod_lnr uses the observed choice", {
  # The two responses have different densities, so passing the wrong one is
  # detectable - this is what catches `dec` being dropped on the floor.
  p0 <- make_prep(y = 0.5, dec = 0, nuzero = 1.5, nuone = -1.5,
                  sigmazero = 0.4, sigmaone = 0.4, ndt = 0.2, poutlier = 0)
  p1 <- make_prep(y = 0.5, dec = 1, nuzero = 1.5, nuone = -1.5,
                  sigmazero = 0.4, sigmaone = 0.4, ndt = 0.2, poutlier = 0)
  expect_false(isTRUE(all.equal(log_lik_cogmod_lnr(1, p0),
                                log_lik_cogmod_lnr(1, p1))))
})


test_that("log_lik_cogmod_lnr errors when dec is missing", {
  prep <- make_prep(y = 0.5, dec = 0, nuzero = 0.5, nuone = 0.2,
                    sigmazero = 0.8, sigmaone = 1.0, ndt = 0.2, poutlier = 0)
  prep$data$dec <- NULL
  expect_error(log_lik_cogmod_lnr(1, prep), "dec")
})


test_that("posterior_predict_cogmod_lnr excludes outliers by default", {
  set.seed(3)
  prep <- make_prep(y = 0.5, dec = 0, nuzero = 0.5, nuone = 0.2,
                    sigmazero = 0.8, sigmaone = 1.0, ndt = 0.4,
                    poutlier = 0.9, n_draws = 500)
  out <- posterior_predict_cogmod_lnr(1, prep)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(500L, 2L))
  # poutlier is 0.9, so with the outlier component on nearly every draw would
  # be below ndt; with it off, none can be
  expect_true(all(out[, 1] > 0.4))
  expect_true(all(out[, 2] %in% c(0, 1)))

  with_out <- posterior_predict_cogmod_lnr(1, prep, predict_outliers = TRUE)
  expect_true(any(with_out[, 1] < 0.4))
})


test_that("the family flag drives predictions when no argument is given", {
  set.seed(4)
  prep <- make_prep(y = 0.5, dec = 0, nuzero = 0.5, nuone = 0.2,
                    sigmazero = 0.8, sigmaone = 1.0, ndt = 0.4,
                    poutlier = 0.9, n_draws = 500,
                    family = cogmod_lnr(predict_outliers = TRUE))
  expect_true(any(posterior_predict_cogmod_lnr(1, prep)[, 1] < 0.4))
  # an explicit argument still wins
  expect_true(all(
    posterior_predict_cogmod_lnr(1, prep, predict_outliers = FALSE)[, 1] > 0.4
  ))
})


test_that("posterior_epred_cogmod_lnr refuses rather than guessing", {
  expect_error(posterior_epred_cogmod_lnr(NULL), "posterior_predict")
})


# family ------------------------------------------------------------------

test_that("cogmod_lnr() builds a valid brms custom family", {
  fam <- cogmod_lnr()
  expect_s3_class(fam, "customfamily")
  expect_equal(fam$dpars,
               c("mu", "nuone", "sigmazero", "sigmaone", "ndt", "poutlier"))
  expect_equal(unname(cogmod:::.family_links(fam)),
               c("identity", "identity", "softplus", "softplus", "log",
                 "logit"))
  expect_equal(fam$vars, "dec[n]")
  expect_false(fam$predict_outliers)
  # no tau / minrt dpars survive the migration
  expect_false(any(c("tau", "minrt") %in% fam$dpars))
})


test_that("cogmod_lnr is one of the outlier-mixture families", {
  expect_true("cogmod_lnr" %in% cogmod:::.OUTLIER_FAMILIES)
  expect_true("cogmod_lnr" %in% cogmod:::.cogmod_families())
  expect_length(cogmod:::.cogmod_families(),
                length(unique(cogmod:::.cogmod_families())))
})


test_that("predict_outliers survives the validation brm() performs", {
  fam <- brms:::validate_family(cogmod_lnr(predict_outliers = TRUE))
  expect_true(fam$predict_outliers)
})


# Stan code ---------------------------------------------------------------

test_that("stanvars carry the likelihood with the outlier component", {
  code <- cogmod_lnr_stanvars()[[1]]$scode
  expect_true(grepl("real cogmod_lnr_lpdf", code))
  expect_true(grepl("int dec", code))
  expect_true(grepl("log_mix\\(poutlier", code))
  expect_true(grepl("half Normal with scale 0.2", code))
  # the old parameterization is gone
  expect_false(grepl("real tau", code))
  expect_false(grepl("real minrt", code))
})


test_that("the 1/K is folded into the generated outlier constant", {
  # The constant is log(2 * dnorm(0, 0, 0.2) / K): the outlier
  # log-density at Y = 0, thinned by the K response options. Getting the 1/K
  # wrong is the silent failure this whole file exists to catch. (It is written
  # as a limit rather than via .dcontam(), which is -Inf at 0 exactly - the
  # component's support is open there, and Stan rejects Y <= 0 anyway.)
  code <- cogmod_lnr_stanvars()[[1]]$scode
  const <- as.numeric(sub(
    ".*real lp_out = ([-0-9.e+]+) -.*", "\\1",
    gsub("\n", " ", code)
  ))
  expect_equal(const, log(2 * stats::dnorm(0, 0, 0.2) / 2),
               tolerance = 1e-12)
  # and it agrees with .dcontam() just inside the support
  expect_equal(
    const - 12.5 * (1e-8)^2,
    log(cogmod:::.dcontam(1e-8) / 2),
    tolerance = 1e-12
  )
})


test_that("Stan cogmod_lnr_lpdf matches dcogmod_lnr", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_lnr")
  grid <- expand.grid(
    Y = c(0.02, 0.2, 0.5, 1.5, 8),
    nuzero = c(-0.5, 1.2),
    nuone = c(-0.2, 0.8),
    sigmazero = c(0.4, 1.1),
    sigmaone = c(0.6, 1.3),
    ndt = c(0, 0.15, 0.4),
    poutlier = c(0, 0.001, 0.4),
    dec = 0:1
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    stan <- lpdf(g$Y, g$nuzero, g$nuone, g$sigmazero, g$sigmaone, g$ndt,
                 g$poutlier, as.integer(g$dec))
    r <- dcogmod_lnr(g$Y, g$nuzero, g$nuone, g$sigmazero, g$sigmaone, g$ndt,
                     response = g$dec, poutlier = g$poutlier, log = TRUE)
    if (is.infinite(r)) {
      expect_equal(stan, r)
    } else {
      # Relative, not absolute: these log-densities reach large magnitudes in
      # the far tails, where one ulp is already ~1e-7 absolute.
      expect_lt(abs(stan - r) / max(1, abs(r)), 1e-10)
    }
  }

  # invalid arguments are rejected on both sides
  expect_equal(lpdf(0.5, 0, 0, 0, 1, 0.2, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 0, 0, 1, 1, 0.2, 1.5, 0L), -Inf)
  expect_equal(lpdf(0.5, 0, 0, 1, 1, -0.1, 0.02, 0L), -Inf)
  expect_equal(lpdf(0.5, 0, 0, 1, 1, 0.2, 0.02, 2L), -Inf)
  expect_equal(lpdf(0, 0, 0, 1, 1, 0.2, 0.02, 0L), -Inf)
})


# priors and inits --------------------------------------------------------

test_that("brms leaves the ndt and poutlier intercepts flat", {
  # The reason cogmod_priors() is not a convenience for this family.
  set.seed(5)
  sim <- rcogmod_lnr(100, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  f <- brms::bf(RT | dec(Error) ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_lnr())
  p <- brms::get_prior(f, data = d, family = cogmod_lnr())
  expect_true(all(!nzchar(p$prior[p$dpar %in% c("ndt", "poutlier") &
                                    p$class == "Intercept"])))
})


test_that("cogmod_priors fills ndt and poutlier for cogmod_lnr", {
  set.seed(6)
  sim <- rcogmod_lnr(150, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150),
                  id = rep(1:10, length.out = 150))

  modelled <- brms::bf(RT | dec(Error) ~ 1, nuone ~ 1, sigmazero ~ 1,
                       sigmaone ~ 1, ndt ~ 1, poutlier ~ 1,
                       family = cogmod_lnr())
  p <- cogmod_priors(modelled, d)
  expect_true(any(p$dpar == "ndt" & p$class == "Intercept" &
                    p$prior == "normal(-1.2, 0.2)"))
  expect_true(any(p$dpar == "poutlier" & p$class == "Intercept" &
                    p$prior == "normal(-5, 1)"))

  # Omitted dpars live on the natural scale and need different priors. brms
  # recognises the name `ndt` and supplies uniform(0, min_Y) - precisely the
  # min-RT bound this parameterization exists to remove - so that row arrives
  # non-empty and has to be overridden rather than filled.
  omitted <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_lnr())
  expect_true(any(grepl("uniform",
                        brms::get_prior(omitted, data = d,
                                        family = cogmod_lnr())$prior)))
  p2 <- cogmod_priors(omitted, d)
  expect_false(any(grepl("uniform", p2$prior)))
  expect_true(any(p2$class == "ndt" & p2$prior == "lognormal(-1.2, 0.2)"))
  expect_true(any(p2$class == "poutlier" & p2$prior == "exponential(100)"))

  # and a mixed formula with group-level terms still builds a Stan program
  mixed <- brms::bf(RT | dec(Error) ~ Condition + (1 | id), ndt ~ Condition,
                    family = cogmod_lnr())
  p3 <- cogmod_priors(mixed, d)
  expect_s3_class(p3, "brmsprior")
  code <- brms::make_stancode(mixed, data = d, prior = p3,
                              stanvars = cogmod_stanvars(mixed))
  expect_false(grepl("uniform\\(", code))
  expect_true(grepl("cogmod_lnr_lpdf", code))
})


test_that("cogmod_priors fences off the losing accumulator's flat direction", {
  # Push an accumulator's rate down far enough and it stops finishing first
  # ever. The density then depends on it only through the loser's survival
  # term, which has already saturated at 1, so the log-likelihood goes exactly
  # flat - and that accumulator's sigma goes with it, having nothing left to act
  # on. Both are directions brms leaves improper, the same failure the function
  # already prevents for `ndt` and `poutlier`.
  set.seed(11)
  sim <- rcogmod_lnr(400, nuzero = 1.2, nuone = 0.2, sigmazero = 0.5,
                     sigmaone = 0.6, ndt = 0.25, poutlier = 0.02)
  ll <- function(d, nuone, sigmaone = 0.6) {
    sum(dcogmod_lnr(d$rt, nuzero = 1.2, nuone = nuone, sigmazero = 0.5,
                    sigmaone = sigmaone, ndt = 0.25, response = d$response,
                    poutlier = 0.02, log = TRUE))
  }

  # the plateau, and the fact that it is reached at a finite nuone
  only0 <- sim[sim$response == 0, ]
  expect_equal(ll(only0, -10), ll(only0, -50), tolerance = 1e-10)
  # sigmaone is unidentified out there, so a prior on nuone alone is not enough
  expect_equal(ll(only0, -20, sigmaone = 0.2), ll(only0, -20, sigmaone = 2),
               tolerance = 1e-10)
  # and the outlier component means both responses being observed does not save
  # it: the trials the retreating accumulator can no longer explain are floored
  expect_equal(ll(sim, -10), ll(sim, -20), tolerance = 1e-10)

  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = nrow(sim)))
  f <- brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition,
                sigmazero ~ Condition, sigmaone ~ Condition, ndt ~ 1,
                family = cogmod_lnr())

  # brms leaves all three flat on its own
  raw <- brms::get_prior(f, data = d, family = cogmod_lnr())
  flat <- raw$dpar %in% c("nuone", "sigmazero", "sigmaone")
  expect_true(all(raw$prior[flat] == ""))

  p <- cogmod_priors(f, d)
  pick <- function(dp, cls) p$prior[p$dpar == dp & p$class == cls & nzchar(p$prior)]
  expect_equal(pick("nuone", "Intercept"), "normal(0.7, 1.5)")
  expect_equal(pick("sigmazero", "Intercept"), "normal(0, 1)")
  expect_equal(pick("sigmaone", "Intercept"), "normal(0, 1)")
  for (dp in c("nuone", "sigmazero", "sigmaone")) {
    expect_equal(pick(dp, "b"), "normal(0, 0.5)", label = paste(dp, "slope"))
  }

  # `mu` is nuzero and has the mirror-image plateau, but it is the response's
  # own intercept, so brms already supplies a proper default there and
  # cogmod_priors leaves it alone
  expect_match(p$prior[p$class == "Intercept" & !nzchar(p$dpar)], "student_t")

  # omitted dpars live on their natural scale: identity for nuone, positive for
  # the two sigmas
  f2 <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_lnr())
  p2 <- cogmod_priors(f2, d)
  expect_equal(p2$prior[p2$class == "nuone"], "normal(0.7, 1.5)")
  expect_equal(p2$prior[p2$class == "sigmazero"], "lognormal(-0.7, 0.75)")
  expect_equal(p2$prior[p2$class == "sigmaone"], "lognormal(-0.7, 0.75)")

  for (form in list(f, f2)) {
    expect_silent(
      brms::make_stancode(form, data = d, family = cogmod_lnr(),
                          prior = cogmod_priors(form, d),
                          stanvars = cogmod_stanvars(form))
    )
  }
})


test_that("cogmod_inits covers the declared parameters", {
  set.seed(8)
  sim <- rcogmod_lnr(150, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150))
  f <- brms::bf(RT | dec(Error) ~ Condition, nuone ~ 1, sigmazero ~ 1,
                sigmaone ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_lnr())
  vals <- cogmod_inits(f, d)(1)
  code <- suppressWarnings(
    brms::make_stancode(f, data = d, family = cogmod_lnr())
  )
  declared <- vapply(cogmod:::.stan_param_decls(code), `[[`, character(1),
                     "name")
  expect_setequal(names(vals), declared)
  expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))))
  # ndt starts small rather than at exp(0) = 1s, which is the point
  expect_lt(exp(cogmod_inits(f, d, jitter = 0)(1)$Intercept_ndt), 0.3)
})


# a real fit --------------------------------------------------------------

test_that("cogmod_lnr recovers ndt above the fastest observed response", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n <- 2000
  sim <- rcogmod_lnr(n, nuzero = 0.6, nuone = 0.1, sigmazero = 0.7,
                     sigmaone = 0.9, ndt = 0.25, poutlier = 0.04)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  # the case the old tau * minrt parameterization could not express
  expect_lt(min(d$RT), 0.25)

  f <- brms::bf(RT | dec(Error) ~ 1, nuone ~ 1, sigmazero ~ 1, sigmaone ~ 1,
                ndt ~ 1, poutlier ~ 1, family = cogmod_lnr())
  fit <- brms::brm(f, data = d, prior = cogmod_priors(f, d),
                   init = cogmod_inits(f, d), stanvars = cogmod_stanvars(f),
                   algorithm = "pathfinder", refresh = 0,
                   backend = "cmdstanr")

  s <- brms::posterior_summary(fit)[, "Estimate"]
  expect_equal(s[["b_Intercept"]], 0.6, tolerance = 0.15)
  expect_equal(s[["b_nuone_Intercept"]], 0.1, tolerance = 0.15)
  expect_equal(log1p(exp(s[["b_sigmazero_Intercept"]])), 0.7, tolerance = 0.15)
  expect_equal(log1p(exp(s[["b_sigmaone_Intercept"]])), 0.9, tolerance = 0.15)
  expect_equal(exp(s[["b_ndt_Intercept"]]), 0.25, tolerance = 0.05)
  expect_equal(plogis(s[["b_poutlier_Intercept"]]), 0.04, tolerance = 0.04)
  # ndt is estimated well above the fastest response, not pinned to it
  expect_gt(exp(s[["b_ndt_Intercept"]]), min(d$RT) * 10)

  ll <- brms::log_lik(fit, ndraws = 10)
  expect_equal(dim(ll), c(10L, n))
  expect_true(all(is.finite(ll)))

  po <- p_outlier(fit)
  expect_equal(nrow(po), n)
  expect_true(all(po$p_outlier >= 0 & po$p_outlier <= 1))
  # responses faster than ndt can only have come from the outlier component
  expect_true(all(po$p_outlier[po$rt < 0.1] > 0.99))
  expect_equal(mean(po$p_outlier), 0.04, tolerance = 0.04)

  expect_true(with_outliers(fit)$family$predict_outliers)
  expect_false(without_outliers(fit)$family$predict_outliers)
})
