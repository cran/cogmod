context("RDM - shifted Racing Diffusion Model with outliers")

# Reference implementation of the mixture, written out longhand. The 1/K on the
# outlier component is what keeps the density summing to one over the response
# options; without it the total comes to 1 + poutlier.
ref_dens <- function(y, vzero, vone, boundary, bias, ndt, response, poutlier) {
  d_out <- 2 * stats::dnorm(y, 0, 0.2) / 2
  v_win <- if (response == 0) vzero else vone
  v_los <- if (response == 0) vone else vzero
  d_dec <- if (y - ndt > 0) {
    .dwald(y, v_win, boundary, bias, ndt) * .swald(y, v_los, boundary, bias, ndt)
  } else {
    0
  }
  poutlier * d_out + (1 - poutlier) * d_dec
}

make_prep <- function(y, dec, vzero, vone, boundary, bias, ndt, poutlier,
                      n_draws = 10, family = cogmod_rdm()) {
  structure(
    list(
      data = list(Y = y, dec = dec),
      family = family,
      dpars = list(
        mu = rep(vzero, n_draws),
        driftone = rep(vone, n_draws),
        sigmabias = rep(bias, n_draws),
        boundary = rep(boundary, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}


# dcogmod_rdm -------------------------------------------------------------

test_that("dcogmod_rdm matches the mixture density", {
  pars <- list(vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2, ndt = 0.25)
  for (poutlier in c(0, 0.02, 0.4)) {
    for (response in 0:1) {
      for (y in c(0.05, 0.2, 0.3, 0.6, 1.5, 4)) {
        expect_equal(
          do.call(dcogmod_rdm, c(list(x = y), pars,
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
    c(vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2, ndt = 0.25),
    c(vzero = 0.8, vone = 0.8, boundary = 1.0, bias = 0.05, ndt = 0.10),
    c(vzero = 4.0, vone = 0.0, boundary = 0.3, bias = 0.5, ndt = 0.40),
    c(vzero = 1.2, vone = 3.0, boundary = 2.0, bias = 1.0, ndt = 0.00)
  )
  for (g in grid) {
    for (poutlier in c(0, 0.03, 0.4)) {
      total <- sum(vapply(0:1, function(k) {
        stats::integrate(
          function(t) {
            dcogmod_rdm(t, g[["vzero"]], g[["vone"]], g[["boundary"]],
                        g[["bias"]], g[["ndt"]], response = k,
                        poutlier = poutlier)
          },
          lower = 0, upper = Inf, subdivisions = 3000
        )$value
      }, numeric(1)))
      expect_equal(total, 1, tolerance = 1e-4,
                   info = sprintf("ndt = %.2f, bias = %.2f, poutlier = %.2f",
                                  g[["ndt"]], g[["bias"]], poutlier))
    }
  }
})


test_that("the density still integrates to one as the start-point range shrinks", {
  # The failure mode cogmod_lba1() had before its own migration: the density is
  # a difference quotient in the threshold, divided by the start-point range, so
  # it cancels as that range goes to zero. `.wald_lpdf_core()` switches to the
  # midpoint (plain Wald) limit before that happens, and this is the check that
  # the switch is in the right place.
  for (A in c(0.5, 1e-2, 1e-4, 1e-6)) {
    total <- sum(vapply(0:1, function(k) {
      stats::integrate(
        function(t) {
          dcogmod_rdm(t, vzero = 2.5, vone = 1.6, boundary = 0.5, bias = A,
                      ndt = 0.2, response = k, poutlier = 0.02)
        },
        lower = 0, upper = Inf, subdivisions = 3000
      )$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-4,
                 info = sprintf("bias = %g", A))
  }
})


test_that("responses faster than ndt keep positive density (the whole point)", {
  y <- c(0.001, 0.05, 0.1, 0.19)
  for (poutlier in c(0.02, 0.4)) {
    for (response in 0:1) {
      d <- dcogmod_rdm(y, 2.5, 1.6, 0.5, 0.2, ndt = 0.2, response = response,
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
  expect_equal(dcogmod_rdm(y, 2.5, 1.6, 0.5, 0.2, ndt = 0.2, response = 0,
                           poutlier = 0),
               rep(0, length(y)))
  expect_equal(dcogmod_rdm(y, 2.5, 1.6, 0.5, 0.2, ndt = 0.2, response = 0,
                           poutlier = 0, log = TRUE),
               rep(-Inf, length(y)))

  x <- seq(0.25, 4, length.out = 50)
  for (response in 0:1) {
    expect_equal(
      dcogmod_rdm(x, 2.5, 1.6, 0.5, 0.2, ndt = 0.2, response = response),
      vapply(x, ref_dens, numeric(1), vzero = 2.5, vone = 1.6, boundary = 0.5,
             bias = 0.2, ndt = 0.2, response = response, poutlier = 0),
      tolerance = 1e-12
    )
  }
})


test_that("the log density stays finite just above ndt, where the density is 0", {
  # The property the sampler depends on: immediately above the shift the density
  # really is below the smallest double, but its log is still usable.
  lg <- dcogmod_rdm(0.15 + c(1e-8, 1e-6, 1e-4), vzero = 0.8, vone = 0.6,
                    boundary = 0.5, bias = 0.2, ndt = 0.15, response = 0,
                    poutlier = 0, log = TRUE)
  expect_true(all(is.finite(lg)))
  expect_true(all(diff(lg) > 0))
})


test_that("the marginal density is the sum of the two defective ones", {
  pars <- list(vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2, ndt = 0.2)
  xs <- c(0.05, 0.25, 0.4, 0.8, 1.5, 3)
  for (poutlier in c(0, 0.03)) {
    marg <- do.call(dcogmod_rdm, c(list(x = xs), pars,
                                   list(poutlier = poutlier)))
    d0 <- do.call(dcogmod_rdm, c(list(x = xs), pars,
                                 list(response = 0, poutlier = poutlier)))
    d1 <- do.call(dcogmod_rdm, c(list(x = xs), pars,
                                 list(response = 1, poutlier = poutlier)))
    expect_equal(d0 + d1, marg, tolerance = 1e-12)
    expect_equal(
      do.call(dcogmod_rdm, c(list(x = xs, log = TRUE), pars,
                             list(poutlier = poutlier))),
      log(marg), tolerance = 1e-10
    )
  }
})


test_that("dcogmod_rdm is vectorized over every argument", {
  n <- 7
  d <- dcogmod_rdm(seq(0.3, 1.2, length.out = n),
                   vzero = seq(1, 3, length.out = n),
                   vone = seq(3, 1, length.out = n),
                   boundary = seq(0.3, 0.8, length.out = n),
                   bias = seq(0.05, 0.4, length.out = n),
                   ndt = seq(0.05, 0.25, length.out = n),
                   response = rep(0:1, length.out = n),
                   poutlier = seq(0, 0.1, length.out = n))
  expect_length(d, n)
  expect_true(all(is.finite(d) & d >= 0))
})


test_that("dcogmod_rdm returns 0 density for invalid parameters", {
  args <- list(x = 0.5, vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2,
               ndt = 0.2, response = 0)
  for (bad in list(list(vzero = -0.1), list(vone = -1), list(boundary = 0),
                   list(boundary = -1), list(bias = -1),
                   list(ndt = -0.1), list(poutlier = 1.5),
                   list(response = 2))) {
    expect_warning(d <- do.call(dcogmod_rdm, modifyList(args, bad)),
                   info = names(bad))
    expect_equal(d, 0)
  }
})


test_that("a drift of exactly zero is a legal accumulator, not an invalid one", {
  # Driftless Brownian motion still reaches the threshold with probability one,
  # so a zero-drift accumulator is slow rather than absent. It is the one closed
  # lower bound in either registry.
  expect_silent(d <- dcogmod_rdm(0.5, vzero = 0, vone = 1.6, boundary = 0.5,
                                 bias = 0.2, ndt = 0.2, response = 0))
  expect_true(is.finite(d) && d > 0)
  expect_true(all(is.finite(
    dcogmod_rdm(c(0.3, 1, 5), vzero = 2, vone = 0, boundary = 0.5, bias = 0.2,
                ndt = 0.15, response = 1, log = TRUE)
  )))
})


# pcogmod_rdm -------------------------------------------------------------

test_that("pcogmod_rdm agrees with integrating the marginal density", {
  pars <- list(vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2, ndt = 0.25)
  for (poutlier in c(0, 0.05)) {
    for (q in c(0.15, 0.3, 0.6, 1.5)) {
      lo <- do.call(pcogmod_rdm, c(list(q = q), pars,
                                   list(poutlier = poutlier)))
      num <- stats::integrate(
        function(z) {
          do.call(dcogmod_rdm, c(list(x = z), pars,
                                 list(poutlier = poutlier)))
        },
        lower = 0, upper = q, subdivisions = 3000, rel.tol = 1e-10
      )$value
      expect_equal(lo, num, tolerance = 1e-6,
                   info = sprintf("q = %.2f, poutlier = %.2f", q, poutlier))
    }
  }

  # CDF and survival are complementary, and log.p is consistent
  q <- c(0.3, 0.8, 2.0)
  lo <- do.call(pcogmod_rdm, c(list(q = q), pars, list(poutlier = 0.05)))
  hi <- do.call(pcogmod_rdm, c(list(q = q, lower.tail = FALSE), pars,
                               list(poutlier = 0.05)))
  expect_equal(lo + hi, rep(1, length(q)), tolerance = 1e-12)
  expect_equal(
    do.call(pcogmod_rdm, c(list(q = q, log.p = TRUE), pars,
                           list(poutlier = 0.05))),
    log(lo), tolerance = 1e-10
  )
})


test_that("pcogmod_rdm keeps the far-tail survival in log space", {
  # The mixture must not cost the property `.swald()` was rewritten for: at
  # poutlier = 0 the survival is still finite on the log scale where 1 - CDF is
  # exactly zero on the natural one.
  ls <- pcogmod_rdm(8, vzero = 6, vone = 6, boundary = 0.5, bias = 0.2,
                    ndt = 0.1, poutlier = 0, lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(ls))
  expect_lt(ls, -50)
})


test_that("the defective CDF integrates the defective density", {
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  for (poutlier in c(0, 0.05)) {
    for (k in 0:1) {
      for (q in c(0.3, 0.5, 0.9, 2)) {
        got <- do.call(pcogmod_rdm, c(list(q = q, response = k), pars,
                                      list(poutlier = poutlier)))
        num <- stats::integrate(
          function(z) {
            do.call(dcogmod_rdm, c(list(x = z, response = k), pars,
                                   list(poutlier = poutlier)))
          },
          lower = 0, upper = q, subdivisions = 3000, rel.tol = 1e-11
        )$value
        expect_equal(got, num, tolerance = 1e-6,
                     info = sprintf("k = %d, q = %.2f, poutlier = %.2f",
                                    k, q, poutlier))
      }
    }
  }
})


test_that("the defective CDFs sum to the marginal one", {
  # The marginal CDF is a closed form and the defective ones are quadrature, so
  # this checks the two paths against each other rather than each against
  # itself.
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  for (poutlier in c(0, 0.05)) {
    for (q in c(0.3, 0.6, 1.5)) {
      per_k <- vapply(0:1, function(k) {
        do.call(pcogmod_rdm, c(list(q = q, response = k), pars,
                               list(poutlier = poutlier)))
      }, numeric(1))
      marg <- do.call(pcogmod_rdm, c(list(q = q), pars,
                                     list(poutlier = poutlier)))
      expect_equal(sum(per_k), marg, tolerance = 1e-7,
                   info = sprintf("q = %.2f, poutlier = %.2f", q, poutlier))
    }
  }
})


test_that("the defective CDF tends to the response probability", {
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  for (poutlier in c(0, 0.1)) {
    pk <- vapply(0:1, function(k) {
      do.call(pcogmod_rdm, c(list(q = Inf, response = k), pars,
                             list(poutlier = poutlier)))
    }, numeric(1))
    expect_equal(sum(pk), 1, tolerance = 1e-7)

    # The upper tail is integrated from q upwards rather than subtracted, so the
    # two sides have to be checked to add back up.
    for (k in 0:1) {
      for (q in c(0.4, 1, 3)) {
        lo <- do.call(pcogmod_rdm, c(list(q = q, response = k), pars,
                                     list(poutlier = poutlier)))
        hi <- do.call(pcogmod_rdm, c(list(q = q, response = k,
                                          lower.tail = FALSE), pars,
                                     list(poutlier = poutlier)))
        expect_equal(lo + hi, pk[k + 1], tolerance = 1e-7)
      }
    }
  }
  # log.p is the log of it.
  expect_equal(
    do.call(pcogmod_rdm, c(list(q = 0.6, response = 1, log.p = TRUE), pars)),
    log(do.call(pcogmod_rdm, c(list(q = 0.6, response = 1), pars))),
    tolerance = 1e-12
  )
})


test_that("pcogmod_rdm is 0 below the support and NA for a missing quantile", {
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  q <- c(-1, 0, NA)
  # Marginally, `.log_mix()` maps NA onto a log-survival of -Inf, which would
  # otherwise come back out here as a CDF of exactly 1.
  expect_equal(do.call(pcogmod_rdm, c(list(q = q), pars,
                                      list(poutlier = 0.02))),
               c(0, 0, NA))
  expect_equal(do.call(pcogmod_rdm, c(list(q = q, response = 0), pars,
                                      list(poutlier = 0.02))),
               c(0, 0, NA))
  # And the survival at the other end of the same range.
  expect_equal(do.call(pcogmod_rdm, c(list(q = q, lower.tail = FALSE), pars,
                                      list(poutlier = 0.02))),
               c(1, 1, NA))
})


test_that("the defective CDF matches the simulated joint distribution", {
  skip_on_cran()
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  set.seed(7)
  d <- do.call(rcogmod_rdm, c(list(n = 2e5), pars, list(poutlier = 0.02)))
  for (k in 0:1) {
    for (q in c(0.4, 0.6, 1, 2)) {
      expect_equal(
        mean(d$rt <= q & d$response == k),
        do.call(pcogmod_rdm, c(list(q = q, response = k), pars,
                               list(poutlier = 0.02))),
        tolerance = 0.005,
        info = sprintf("k = %d, q = %.2f", k, q)
      )
    }
  }
})


# qcogmod_rdm -------------------------------------------------------------

test_that("qcogmod_rdm inverts pcogmod_rdm", {
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  pr <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  # Marginal: p is an ordinary probability.
  q <- do.call(qcogmod_rdm, c(list(p = pr), pars, list(poutlier = 0.02)))
  expect_equal(do.call(pcogmod_rdm, c(list(q = q), pars,
                                      list(poutlier = 0.02))),
               pr, tolerance = 1e-6)
  expect_false(is.unsorted(q))

  # Per response, scaled: p is a fraction of that response's own mass, which is
  # what a quantile-probability plot asks for.
  for (k in 0:1) {
    q <- do.call(qcogmod_rdm, c(list(p = pr, response = k, scale_p = TRUE),
                                pars, list(poutlier = 0.02)))
    mass <- do.call(pcogmod_rdm, c(list(q = Inf, response = k), pars,
                                   list(poutlier = 0.02)))
    got <- do.call(pcogmod_rdm, c(list(q = q, response = k), pars,
                                  list(poutlier = 0.02)))
    expect_equal(got / mass, pr, tolerance = 1e-6)
  }

  # Unscaled, p is read off the defective CDF directly.
  k <- 1
  mass <- do.call(pcogmod_rdm, c(list(q = Inf, response = k), pars))
  target <- c(0.05, 0.15, 0.25)
  q <- do.call(qcogmod_rdm, c(list(p = target, response = k), pars))
  expect_equal(do.call(pcogmod_rdm, c(list(q = q, response = k), pars)),
               target, tolerance = 1e-6)
  expect_lt(max(target), mass)
})


test_that("qcogmod_rdm handles the edges of the attainable range", {
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)

  # Zero is the infimum, and the supremum is only reached in the limit.
  expect_equal(do.call(qcogmod_rdm, c(list(p = c(0, 1), response = 0,
                                           scale_p = TRUE), pars)),
               c(0, Inf))

  # Above the response's own mass there is no quantile at all.
  expect_warning(
    got <- do.call(qcogmod_rdm, c(list(p = 0.9, response = 1), pars)),
    "exceeds the attainable probability"
  )
  expect_true(is.na(got))

  # The bracket only affects speed: an ndt above it, and quantiles well past it,
  # both still come back.
  slow <- list(vzero = 0.6, vone = 0.5, boundary = 1.5, bias = 0.3, ndt = 12)
  q <- do.call(qcogmod_rdm, c(list(p = 0.9, interval = c(0, 1)), slow))
  expect_gt(q, 12)
  expect_equal(do.call(pcogmod_rdm, c(list(q = q), slow)), 0.9,
               tolerance = 1e-6)
  expect_error(do.call(qcogmod_rdm, c(list(p = 0.5, interval = 3), pars)),
               "two increasing")

  # lower.tail and log.p are the same request, differently written.
  expect_equal(
    do.call(qcogmod_rdm, c(list(p = 0.75), pars)),
    do.call(qcogmod_rdm, c(list(p = 0.25, lower.tail = FALSE), pars))
  )
  expect_equal(
    do.call(qcogmod_rdm, c(list(p = 0.4), pars)),
    do.call(qcogmod_rdm, c(list(p = log(0.4), log.p = TRUE), pars))
  )
})


test_that("qcogmod_rdm recovers the quantiles of simulated data", {
  skip_on_cran()
  pars <- list(vzero = 2.5, vone = 1.5, boundary = 0.5, bias = 0.3, ndt = 0.25)
  pr <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  set.seed(11)
  d <- do.call(rcogmod_rdm, c(list(n = 2e5), pars, list(poutlier = 0.02)))
  for (k in 0:1) {
    got <- do.call(qcogmod_rdm, c(list(p = pr, response = k, scale_p = TRUE),
                                  pars, list(poutlier = 0.02)))
    emp <- as.numeric(stats::quantile(d$rt[d$response == k], pr))
    expect_equal(got, emp, tolerance = 0.01,
                 info = sprintf("response %d", k))
  }
})


# Agreement with rtdists --------------------------------------------------

test_that("the density and defective CDF agree with rtdists", {
  # rtdists gained the RDM in 2026 (rtdists/rtdists#23), from the same Tillman
  # et al. (2020) equations but a separate lineage - Van Zandt's digt/pigt, on
  # the natural scale. An independent implementation is worth more than any
  # self-consistency check, so compare where both are accurate.
  #
  # Their `b` is the absolute threshold, ours the offset: b = boundary + bias.
  # Their responses are 1-based, ours 0-based. Nothing else differs.
  #
  # Only ordinary parameters are compared. Their natural-scale arithmetic
  # returns exactly 0 for a density that underflows - at drift 4 it does so from
  # about t = 5 - which is the whole reason ours is in log space, so the tails
  # are not a like-for-like comparison and are left out.
  skip_on_cran()
  skip_if_not_installed("rtdists")
  skip_if_not(exists("dRDM", envir = asNamespace("rtdists")),
              "installed rtdists predates the RDM")

  dRDM <- get("dRDM", envir = asNamespace("rtdists"))
  pRDM <- get("pRDM", envir = asNamespace("rtdists"))

  A <- 0.3
  boundary <- 0.5
  ndt <- 0.25
  v <- c(2.5, 1.5)
  rt <- seq(0.35, 2.5, length.out = 25)

  for (k in 0:1) {
    theirs <- dRDM(rt, rep(k + 1L, length(rt)), A = A, b = boundary + A,
                   t0 = ndt, v = v, silent = TRUE)
    ours <- dcogmod_rdm(rt, vzero = v[1], vone = v[2], boundary = boundary,
                        bias = A, ndt = ndt, response = k, poutlier = 0)
    expect_equal(ours, theirs, tolerance = 1e-8,
                 info = sprintf("density, response %d", k))

    theirs <- pRDM(rt, rep(k + 1L, length(rt)), A = A, b = boundary + A,
                   t0 = ndt, v = v, silent = TRUE)
    ours <- pcogmod_rdm(rt, vzero = v[1], vone = v[2], boundary = boundary,
                        bias = A, ndt = ndt, response = k, poutlier = 0)
    expect_equal(ours, theirs, tolerance = 1e-6,
                 info = sprintf("defective CDF, response %d", k))
  }

  # A zero start-point range and a zero drift are legal for both.
  expect_equal(
    dcogmod_rdm(rt, vzero = v[1], vone = v[2], boundary = boundary, bias = 0,
                ndt = ndt, response = 0, poutlier = 0),
    dRDM(rt, rep(1L, length(rt)), A = 0, b = boundary, t0 = ndt, v = v,
         silent = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dcogmod_rdm(rt, vzero = 0, vone = v[2], boundary = boundary, bias = A,
                ndt = ndt, response = 1, poutlier = 0),
    dRDM(rt, rep(2L, length(rt)), A = A, b = boundary + A, t0 = ndt,
         v = c(0, v[2]), silent = TRUE),
    tolerance = 1e-8
  )
})


# rcogmod_rdm -------------------------------------------------------------

test_that("rcogmod_rdm returns rt and response", {
  set.seed(1)
  sim <- rcogmod_rdm(100, 2.5, 1.6, 0.5, 0.2, ndt = 0.2)
  expect_true(is.data.frame(sim))
  expect_named(sim, c("rt", "response"))
  expect_equal(nrow(sim), 100)
  expect_true(all(sim$response %in% c(0, 1)))
  expect_true(all(is.finite(sim$rt)))
  # vzero > vone, so accumulator 0 wins more often
  expect_gt(mean(sim$response == 0), 0.5)
})


test_that("rcogmod_rdm handles a zero drift", {
  # A zero-drift accumulator is slow, but it still finishes, so it can win --
  # and the RTs stay finite.
  set.seed(456)
  sim <- rcogmod_rdm(10000, vzero = 0.8, vone = 0, boundary = 0.5, bias = 0.2,
                     ndt = 0.15)
  expect_gt(mean(sim$response == 0), 0.5)
  expect_gt(mean(sim$response == 1), 0)
  expect_true(all(sim$rt >= 0.15))
  expect_true(all(is.finite(sim$rt)))
})


test_that("rcogmod_rdm reproduces its own density", {
  set.seed(42)
  n <- 30000
  pars <- list(vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2, ndt = 0.25)
  for (poutlier in c(0, 0.05)) {
    sim <- do.call(rcogmod_rdm, c(list(n = n), pars,
                                  list(poutlier = poutlier)))
    for (k in 0:1) {
      integrated <- stats::integrate(
        function(t) {
          do.call(dcogmod_rdm, c(list(x = t), pars,
                                 list(response = k, poutlier = poutlier)))
        },
        lower = 0, upper = Inf, subdivisions = 3000
      )$value
      expect_equal(mean(sim$response == k), integrated, tolerance = 0.015,
                   info = sprintf("response %d, poutlier %.2f", k, poutlier))
    }
    # empirical quantiles against the numeric CDF
    for (pr in c(0.1, 0.5, 0.9)) {
      qq <- unname(stats::quantile(sim$rt, pr))
      cdf <- do.call(pcogmod_rdm, c(list(q = qq), pars,
                                    list(poutlier = poutlier)))
      expect_equal(cdf, pr, tolerance = 0.015,
                   info = sprintf("p %.1f, poutlier %.2f", pr, poutlier))
    }
  }
})


test_that("rcogmod_rdm produces responses below ndt only via outliers", {
  set.seed(2)
  expect_true(all(rcogmod_rdm(2000, ndt = 0.4, poutlier = 0)$rt > 0.4))
  sim <- rcogmod_rdm(20000, ndt = 0.4, poutlier = 0.3)
  expect_true(any(sim$rt < 0.4))
  # the contaminant produces a choice as well as a time, uniform over the two
  fast <- sim$response[sim$rt < 0.02]
  expect_setequal(unique(fast), c(0, 1))
  expect_equal(mean(fast == 0), 0.5, tolerance = 0.1)
})


test_that("rcogmod_rdm errors on invalid parameters", {
  expect_error(rcogmod_rdm(10, vzero = -1), "mu")
  expect_error(rcogmod_rdm(10, boundary = 0), "boundary")
  expect_error(rcogmod_rdm(10, bias = -1), "sigmabias")
  expect_error(rcogmod_rdm(10, ndt = -1), "ndt")
  expect_error(rcogmod_rdm(10, poutlier = 2), "poutlier")
})


# The Wald building blocks ------------------------------------------------

test_that("the Wald survival stays finite where 1 - CDF would underflow", {
  # This is the regression that motivated the log-space rewrite. Computing the
  # survival as 1 - CDF returns exactly 0 here, which would give the sampler a
  # -Inf log-likelihood with a zero gradient.
  n1 <- function(v) rep(v, 1)
  for (drift in c(3, 6, 10)) {
    for (t in c(4, 8, 20)) {
      ls <- .swald(t, n1(drift), n1(0.5), n1(0.2), n1(0), log.p = TRUE)
      expect_true(is.finite(ls),
                  label = sprintf("log survival finite at drift=%g, t=%g", drift, t))
      expect_lt(ls, 0)
    }
  }
  # ... and it is still the right number, matching numerical integration.
  dens <- function(z) .dwald(z, rep(3, length(z)), rep(0.5, length(z)),
                             rep(0.2, length(z)), rep(0, length(z)))
  num <- integrate(dens, 3, Inf, rel.tol = 1e-11, subdivisions = 2000)$value
  expect_equal(.swald(3, 3, 0.5, 0.2, 0), num, tolerance = 1e-8)
})


test_that("the Wald density and CDF survive extreme drift x threshold", {
  # exp(2 * drift * b) overflows on its own here; folding it into a single
  # exponent with the log-CDF keeps the product representable.
  for (prm in list(c(80, 4, 1), c(200, 1.6, 0.4), c(500, 4, 1))) {
    ld <- .dwald(0.5, prm[1], prm[2], prm[3], 0.1, log = TRUE)
    ls <- .swald(0.5, prm[1], prm[2], prm[3], 0.1, log.p = TRUE)
    expect_true(is.finite(ld) && !is.nan(ld))
    expect_true(is.finite(ls) && !is.nan(ls))
  }
})


test_that("the driftless Wald CDF is correct (regression: was off by a factor of 2)", {
  # `.pwald()`'s drift == 0 branch used to return half the right answer. It was
  # unreachable from dcogmod_rdm() at the time, but the Stan version reaches it
  # whenever the drift approaches zero.
  k <- 0.5
  A <- 0.2
  for (t in c(0.3, 1.0, 3.0)) {
    num <- integrate(function(z) .dwald(z, rep(0, length(z)), rep(k, length(z)),
                                        rep(A, length(z)), rep(0, length(z))),
                     0, t, rel.tol = 1e-11, subdivisions = 2000)$value
    expect_equal(.pwald(t, 0, k, A, 0), num, tolerance = 1e-8,
                 label = sprintf("driftless CDF at t=%g", t))
  }
  # And it is continuous as the drift crosses zero
  s <- vapply(c(1e-5, 1e-7, 1e-9, 0, -1e-9), function(v) {
    .swald(1.1, v, 1, 0.3, 0.1, log.p = TRUE)
  }, numeric(1))
  expect_lt(max(abs(diff(s))), 1e-5)
})


test_that("the Wald density is the published one (Tillman et al., 2020)", {
  # Everything else in this file checks the implementation against itself - it
  # integrates to one, the survival complements the CDF, Stan matches R. All of
  # that would pass a log-space rewrite that was self-consistently wrong. These
  # three pin it to the equations as printed, in the paper's own notation, where
  # `b` is the absolute threshold and `A` the start-point range, so
  # `boundary = b - A`.
  rt <- seq(0.05, 3, length.out = 50)
  b <- 1
  A <- 0.4
  v <- 2
  ours <- function(t, A, b, v) {
    n <- length(t)
    .dwald(t, rep(v, n), rep(b - A, n), rep(A, n), rep(0, n))
  }

  # Equation 5: with start-point variability.
  alpha <- (b - A - rt * v) / sqrt(rt)
  beta <- (b - rt * v) / sqrt(rt)
  eq5 <- (1 / A) * (-v * pnorm(alpha) + dnorm(alpha) / sqrt(rt) +
                      v * pnorm(beta) - dnorm(beta) / sqrt(rt))
  expect_equal(ours(rt, A, b, v), eq5, tolerance = 1e-12)

  # Equation 6: start-point variability with a drift rate of zero.
  eq6 <- (1 / A) * (dnorm((b - A) / sqrt(rt)) / sqrt(rt) -
                      dnorm(b / sqrt(rt)) / sqrt(rt))
  expect_equal(ours(rt, A, b, 0), eq6, tolerance = 1e-12)

  # Equation 2: the plain shifted Wald, which is what the model becomes when the
  # start point stops varying. Exact at bias = 0, and approached from above.
  eq2 <- b * (2 * pi * rt^3)^(-1 / 2) * exp(-(1 / (2 * rt)) * (v * rt - b)^2)
  expect_equal(ours(rt, 0, b, v), eq2, tolerance = 1e-12)
  expect_lt(max(abs(ours(rt, 1e-5, b, v) - eq2) / eq2), 1e-4)
})


test_that("a start-point range of exactly zero is the plain Wald", {
  # `bias = 0` is a model - both accumulators start at 0 every trial - not an
  # invalid parameter, exactly as for cogmod_lba1() and cogmod_lba2().
  k <- 1
  v <- 2
  t <- c(0.05, 0.3, 0.7, 1.5, 3)
  n <- length(t)
  z0 <- rep(0, n)

  # Density, survival and CDF all reduce to the textbook Wald.
  expect_equal(
    .dwald(t, rep(v, n), rep(k, n), z0, z0),
    k * (2 * pi * t^3)^(-1 / 2) * exp(-(1 / (2 * t)) * (v * t - k)^2),
    tolerance = 1e-12
  )
  wald_surv <- pnorm((k - v * t) / sqrt(t)) -
    exp(2 * v * k) * pnorm(-(k + v * t) / sqrt(t))
  expect_equal(.swald(t, rep(v, n), rep(k, n), z0, z0), wald_surv,
               tolerance = 1e-10)
  expect_equal(.pwald(t, rep(v, n), rep(k, n), z0, z0), 1 - wald_surv,
               tolerance = 1e-10)

  # The driftless case goes through `.wald_lcdf_direct()`, which had no branch
  # that could survive a zero range - both of its others divide by it.
  expect_equal(.pwald(t, z0, rep(k, n), z0, z0), 2 * pnorm(-k / sqrt(t)),
               tolerance = 1e-10)

  # And it is continuous coming down to it, rather than a special case bolted on.
  near <- vapply(c(1e-5, 1e-7, 1e-9, 0), function(a) {
    dcogmod_rdm(0.7, bias = a, response = 0, poutlier = 0, log = TRUE)
  }, numeric(1))
  expect_lt(max(abs(diff(near))), 1e-4)

  # It is a legal argument end to end, not just in the internals.
  expect_equal(
    integrate(function(z) dcogmod_rdm(z, bias = 0, poutlier = 0), 0, Inf)$value,
    1, tolerance = 1e-6
  )
  d <- rcogmod_rdm(50, bias = 0)
  expect_true(all(is.finite(d$rt)) && all(d$rt > 0))

  # A negative range is still rejected.
  expect_warning(expect_equal(dcogmod_rdm(0.7, bias = -1, response = 0), 0),
                 "sigmabias")
})


test_that("the Wald density integrates to one across the parameter grid", {
  for (prm in list(c(3, 0.5, 0.2), c(1, 1, 0.5), c(0.5, 2, 0.1),
                   c(6, 0.3, 0.05), c(0, 1, 0.3))) {
    nu <- prm[1]; k <- prm[2]; A <- prm[3]
    tot <- integrate(function(z) .dwald(z, rep(nu, length(z)), rep(k, length(z)),
                                        rep(A, length(z)), rep(0, length(z))),
                     0, Inf, rel.tol = 1e-10, subdivisions = 3000)$value
    expect_equal(tot, 1, tolerance = 1e-6,
                 label = sprintf("total mass at drift=%g, boundary=%g, bias=%g", nu, k, A))
  }
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_rdm matches dcogmod_rdm", {
  prep <- make_prep(y = c(0.5, 0.15), dec = c(1, 0), vzero = 2.5, vone = 1.6,
                    boundary = 0.5, bias = 0.2, ndt = 0.25, poutlier = 0.03)
  for (i in 1:2) {
    expect_equal(
      log_lik_cogmod_rdm(i, prep),
      rep(dcogmod_rdm(prep$data$Y[i], 2.5, 1.6, 0.5, 0.2, 0.25,
                      response = prep$data$dec[i], poutlier = 0.03,
                      log = TRUE), 10)
    )
  }
})


test_that("log_lik_cogmod_rdm uses the observed choice", {
  # The two responses have different densities, so passing the wrong one is
  # detectable - this is what catches `dec` being dropped on the floor.
  p0 <- make_prep(y = 0.5, dec = 0, vzero = 4, vone = 0.5, boundary = 0.5,
                  bias = 0.2, ndt = 0.2, poutlier = 0)
  p1 <- make_prep(y = 0.5, dec = 1, vzero = 4, vone = 0.5, boundary = 0.5,
                  bias = 0.2, ndt = 0.2, poutlier = 0)
  expect_false(isTRUE(all.equal(log_lik_cogmod_rdm(1, p0),
                                log_lik_cogmod_rdm(1, p1))))
})


test_that("log_lik_cogmod_rdm errors when dec is missing", {
  prep <- make_prep(y = 0.5, dec = 0, vzero = 2.5, vone = 1.6, boundary = 0.5,
                    bias = 0.2, ndt = 0.2, poutlier = 0)
  prep$data$dec <- NULL
  expect_error(log_lik_cogmod_rdm(1, prep), "dec")
})


test_that("posterior_predict_cogmod_rdm excludes outliers by default", {
  set.seed(3)
  prep <- make_prep(y = 0.5, dec = 0, vzero = 2.5, vone = 1.6, boundary = 0.5,
                    bias = 0.2, ndt = 0.4, poutlier = 0.9, n_draws = 500)
  out <- posterior_predict_cogmod_rdm(1, prep)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(500L, 2L))
  # poutlier is 0.9, so with the outlier component on nearly every draw would
  # be below ndt; with it off, none can be
  expect_true(all(out[, 1] > 0.4))
  expect_true(all(out[, 2] %in% c(0, 1)))

  with_out <- posterior_predict_cogmod_rdm(1, prep, predict_outliers = TRUE)
  expect_true(any(with_out[, 1] < 0.4))
})


test_that("the family flag drives predictions when no argument is given", {
  set.seed(4)
  prep <- make_prep(y = 0.5, dec = 0, vzero = 2.5, vone = 1.6, boundary = 0.5,
                    bias = 0.2, ndt = 0.4, poutlier = 0.9, n_draws = 500,
                    family = cogmod_rdm(predict_outliers = TRUE))
  expect_true(any(posterior_predict_cogmod_rdm(1, prep)[, 1] < 0.4))
  # an explicit argument still wins
  expect_true(all(
    posterior_predict_cogmod_rdm(1, prep, predict_outliers = FALSE)[, 1] > 0.4
  ))
})


test_that("posterior_epred_cogmod_rdm refuses rather than guessing", {
  expect_error(posterior_epred_cogmod_rdm(NULL), "posterior_predict")
})


# family ------------------------------------------------------------------

test_that("cogmod_rdm() builds a valid brms custom family", {
  fam <- cogmod_rdm()
  expect_s3_class(fam, "customfamily")
  # Order matters: brms passes the dpars to cogmod_rdm_lpdf in exactly this one.
  expect_equal(fam$dpars,
               c("mu", "driftone", "sigmabias", "boundary", "ndt", "poutlier"))
  expect_equal(unname(cogmod:::.family_links(fam)),
               c("softplus", "softplus", "softplus", "softplus", "log",
                 "logit"))
  expect_equal(fam$vars, "dec[n]")
  expect_false(fam$predict_outliers)
  # no tau / minrt dpars survive the migration
  expect_false(any(c("tau", "minrt") %in% fam$dpars))
})


test_that("cogmod_rdm is one of the outlier-mixture families", {
  expect_true("cogmod_rdm" %in% cogmod:::.OUTLIER_FAMILIES)
  expect_true("cogmod_rdm" %in% cogmod:::.cogmod_families())
  expect_length(cogmod:::.cogmod_families(),
                length(unique(cogmod:::.cogmod_families())))
})


test_that("predict_outliers survives the validation brm() performs", {
  fam <- brms:::validate_family(cogmod_rdm(predict_outliers = TRUE))
  expect_true(fam$predict_outliers)
})


# Stan code ---------------------------------------------------------------

test_that("stanvars carry the likelihood with the outlier component", {
  code <- cogmod_rdm_stanvars()[[1]]$scode
  expect_true(grepl("real cogmod_rdm_lpdf", code))
  expect_true(grepl("int dec", code))
  expect_true(grepl("log_mix\\(poutlier", code))
  expect_true(grepl("half Normal with scale 0.2", code))
  # the decision helpers travel with it
  expect_true(grepl("real cogmod_rdm_wald_ldens", code))
  expect_true(grepl("real cogmod_rdm_wald_lsurv", code))
  # the old parameterization is gone
  expect_false(grepl("real tau", code))
  expect_false(grepl("real minrt", code))
})


test_that("the 1/K is folded into the generated outlier constant", {
  # The constant is log(2 * dnorm(0, 0, 0.2) / K): the outlier
  # log-density at Y = 0, thinned by the K response options. Getting the 1/K
  # wrong is the silent failure the density tests above exist to catch.
  code <- cogmod_rdm_stanvars()[[1]]$scode
  const <- as.numeric(sub(
    ".*real lp_out = ([-0-9.e+]+) -.*", "\\1",
    gsub("\n", " ", code)
  ))
  expect_equal(const, log(2 * stats::dnorm(0, 0, 0.2) / 2),
               tolerance = 1e-12)
})


test_that("Stan cogmod_rdm_lpdf matches dcogmod_rdm", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_rdm")
  grid <- covering_grid(
    Y = c(0.02, 0.25, 0.8, 3, 8),
    mu = c(0, 0.5, 2, 6),
    driftone = c(0.3, 1.5, 5),
    sigmabias = c(0.02, 0.3, 1.0),
    boundary = c(0.2, 0.8, 2.0),
    ndt = c(0, 0.15, 0.4),
    poutlier = c(0, 0.001, 0.4),
    dec = 0:1,
    # a drift of exactly zero reaches the driftless branch on both sides, so
    # that slice is swept in full rather than left to the covering subset
    always = function(g) {
      g$mu == 0 & g$sigmabias == 0.3 & g$boundary == 0.8 & g$poutlier == 0.001
    }
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    stan <- lpdf(g$Y, g$mu, g$driftone, g$sigmabias, g$boundary, g$ndt,
                 g$poutlier, as.integer(g$dec))
    r <- dcogmod_rdm(g$Y, g$mu, g$driftone, g$boundary, g$sigmabias, g$ndt,
                     response = g$dec, poutlier = g$poutlier, log = TRUE)
    if (is.infinite(r)) {
      expect_equal(stan, r)
    } else {
      # Relative, not absolute: these log-densities reach large magnitudes in
      # the far tails, where one ulp is already ~1e-7 absolute. The floor is
      # 1e-6 rather than the 1e-10 the simpler families manage, because the
      # decision density is a six-term signed log-sum that cancels: R and Stan
      # accumulate it in the same order but not in the same arithmetic.
      expect_lt(abs(stan - r) / max(1, abs(r)), 1e-6)
    }
  }

  # a zero start-point range is a model on both sides, not an invalid argument
  expect_equal(
    lpdf(0.5, 2, 1, 0, 0.5, 0.2, 0.02, 0L),
    dcogmod_rdm(0.5, vzero = 2, vone = 1, boundary = 0.5, bias = 0,
                ndt = 0.2, response = 0, poutlier = 0.02, log = TRUE),
    tolerance = 1e-8
  )

  # invalid arguments are rejected on both sides
  expect_equal(lpdf(0.5, 2, 1, -0.1, 0.5, 0.2, 0.02, 0L), -Inf) # sigmabias
  expect_equal(lpdf(0.5, 2, 1, 0.2, 0, 0.2, 0.02, 0L), -Inf)   # boundary
  expect_equal(lpdf(0.5, -1, 1, 0.2, 0.5, 0.2, 0.02, 0L), -Inf) # mu
  expect_equal(lpdf(0.5, 2, 1, 0.2, 0.5, -0.1, 0.02, 0L), -Inf) # ndt
  expect_equal(lpdf(0.5, 2, 1, 0.2, 0.5, 0.2, 1.5, 0L), -Inf)  # poutlier
  expect_equal(lpdf(0.5, 2, 1, 0.2, 0.5, 0.2, 0.02, 2L), -Inf) # dec
  expect_equal(lpdf(0, 2, 1, 0.2, 0.5, 0.2, 0.02, 0L), -Inf)   # Y
})


test_that("Stan cogmod_rdm_lpdf is finite in the tails and for fast responses", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_rdm")

  # Slow responses with a fast losing accumulator: the survival term is what
  # collapses here if it is computed as 1 - CDF.
  for (prm in list(c(5, 10), c(8, 10), c(20, 6))) {
    ll <- lpdf(prm[1], prm[2], prm[2], 0.2, 0.5, 0.1, 0, 0L)
    expect_true(is.finite(ll),
                label = sprintf("finite at Y=%g, drift=%g", prm[1], prm[2]))
  }

  # Responses just above the non-decision time, where (boundary - drift*t)/sqrt(t)
  # grows past 8 and Stan's own std_normal_lccdf would return -Inf.
  for (dy in c(1e-2, 1e-3, 1e-4)) {
    ll <- lpdf(0.1 + dy, 2, 1.5, 0.3, 1.0, 0.1, 0, 0L)
    expect_true(is.finite(ll),
                label = sprintf("finite at %g above the ndt", dy))
  }

  # Below the non-decision time only the outlier component can have produced the
  # response - and at poutlier = 0 there is nothing left.
  expect_true(is.finite(lpdf(0.05, 2, 1, 0.3, 0.5, 0.2, 0.02, 0L)))
  expect_equal(lpdf(0.05, 2, 1, 0.3, 0.5, 0.2, 0, 0L), -Inf)
})


# priors and inits --------------------------------------------------------

test_that("brms leaves the ndt and poutlier intercepts flat", {
  # The reason cogmod_priors() is not a convenience for this family.
  set.seed(5)
  sim <- rcogmod_rdm(100, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  f <- brms::bf(RT | dec(Error) ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_rdm())
  p <- brms::get_prior(f, data = d, family = cogmod_rdm())
  expect_true(all(!nzchar(p$prior[p$dpar %in% c("ndt", "poutlier") &
                                    p$class == "Intercept"])))
})


test_that("cogmod_priors fills ndt and poutlier for cogmod_rdm", {
  set.seed(6)
  sim <- rcogmod_rdm(150, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150),
                  id = rep(1:10, length.out = 150))

  modelled <- brms::bf(RT | dec(Error) ~ 1, driftone ~ 1, sigmabias ~ 1,
                       boundary ~ 1, ndt ~ 1, poutlier ~ 1,
                       family = cogmod_rdm())
  p <- cogmod_priors(modelled, d)
  expect_true(any(p$dpar == "ndt" & p$class == "Intercept" &
                    p$prior == "normal(-1.2, 0.2)"))
  expect_true(any(p$dpar == "poutlier" & p$class == "Intercept" &
                    p$prior == "normal(-5, 1)"))

  # Omitted dpars live on the natural scale and need different priors. brms
  # recognises the name `ndt` and supplies uniform(0, min_Y) - precisely the
  # min-RT bound this parameterization exists to remove - so that row arrives
  # non-empty and has to be overridden rather than filled.
  omitted <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_rdm())
  expect_true(any(grepl("uniform",
                        brms::get_prior(omitted, data = d,
                                        family = cogmod_rdm())$prior)))
  p2 <- cogmod_priors(omitted, d)
  expect_false(any(grepl("uniform", p2$prior)))
  expect_true(any(p2$class == "ndt" & p2$prior == "lognormal(-1.2, 0.2)"))
  expect_true(any(p2$class == "poutlier" & p2$prior == "exponential(100)"))

  # and a mixed formula with group-level terms still builds a Stan program
  mixed <- brms::bf(RT | dec(Error) ~ Condition + (1 | id), ndt ~ Condition,
                    family = cogmod_rdm())
  p3 <- cogmod_priors(mixed, d)
  expect_s3_class(p3, "brmsprior")
  code <- brms::make_stancode(mixed, data = d, prior = p3,
                              stanvars = cogmod_stanvars(mixed))
  expect_false(grepl("uniform\\(", code))
  expect_true(grepl("cogmod_rdm_lpdf", code))
})


test_that("cogmod_priors fences off the sigmabias / boundary ridge", {
  # The two enter the threshold only through their sum, b = boundary + sigmabias,
  # so the likelihood is nearly flat along the ridge that holds the sum fixed -
  # and a softplus link reaches zero only at minus infinity. Flat prior plus
  # unbounded flat region is the improper posterior cogmod_priors() exists to
  # prevent, the same failure as cogmod_lba1().
  set.seed(11)
  sim <- rcogmod_rdm(2000, vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2,
                     ndt = 0.25, poutlier = 0.02)
  # Profiled, not sliced: `boundary` has to be free to slide, since it is the
  # sliding that makes the direction flat. Holding it fixed costs 60-odd log
  # units and would hide the ridge entirely.
  nll <- function(par, bias) {
    if (par[1] <= 0 || par[2] < 0 || par[3] < 0) return(1e10)
    -sum(dcogmod_rdm(sim$rt, vzero = par[2], vone = par[3], boundary = par[1],
                     bias = bias, ndt = 0.25, response = sim$response,
                     poutlier = 0.02, log = TRUE))
  }
  prof <- vapply(c(1e-4, 0.2, 0.35), function(a) {
    o <- stats::optim(c(0.6, 2.5, 1.6), nll, bias = a,
                      control = list(reltol = 1e-10))
    c(boundary = o$par[1], nll = o$value)
  }, numeric(2))
  # the whole ridge, from no start-point variability to half the threshold, is
  # worth only a few log units on 2000 trials
  expect_lt(diff(range(prof["nll", ])), 15)
  # ... and `boundary` slides down to pay for it, which is what makes it a ridge
  expect_true(all(diff(prof["boundary", ]) < 0))

  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = nrow(sim)))
  f <- brms::bf(RT | dec(Error) ~ Condition, driftone ~ Condition,
                sigmabias ~ Condition, boundary ~ Condition, ndt ~ 1,
                family = cogmod_rdm())

  # brms leaves both flat on its own
  raw <- brms::get_prior(f, data = d, family = cogmod_rdm())
  expect_true(all(raw$prior[raw$dpar %in% c("sigmabias", "boundary")] == ""))

  p <- cogmod_priors(f, d)
  pick <- function(dp, cls) p$prior[p$dpar == dp & p$class == cls & nzchar(p$prior)]
  expect_equal(pick("sigmabias", "Intercept"), "normal(0, 1)")
  expect_equal(pick("boundary", "Intercept"), "normal(0, 1)")
  for (dp in c("sigmabias", "boundary")) {
    expect_equal(pick(dp, "b"), "normal(0, 0.5)", label = paste(dp, "slope"))
  }

  # omitted dpars live on their natural scale, which is positive here
  f2 <- brms::bf(RT | dec(Error) ~ 1, family = cogmod_rdm())
  p2 <- cogmod_priors(f2, d)
  expect_equal(p2$prior[p2$class == "sigmabias"], "lognormal(-0.7, 0.75)")
  expect_equal(p2$prior[p2$class == "boundary"], "lognormal(-0.7, 0.75)")

  for (form in list(f, f2)) {
    expect_silent(
      brms::make_stancode(form, data = d, family = cogmod_rdm(),
                          prior = cogmod_priors(form, d),
                          stanvars = cogmod_stanvars(form))
    )
  }
})


test_that("cogmod_inits covers the declared parameters", {
  set.seed(8)
  sim <- rcogmod_rdm(150, ndt = 0.25, poutlier = 0.03)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150))
  f <- brms::bf(RT | dec(Error) ~ Condition, driftone ~ 1, sigmabias ~ 1,
                boundary ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_rdm())
  vals <- cogmod_inits(f, d)(1)
  code <- suppressWarnings(
    brms::make_stancode(f, data = d, family = cogmod_rdm())
  )
  declared <- vapply(cogmod:::.stan_param_decls(code), `[[`, character(1),
                     "name")
  expect_setequal(names(vals), declared)
  expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))))
  # ndt starts small rather than at exp(0) = 1s, which is the point
  v0 <- cogmod_inits(f, d, jitter = 0)(1)
  expect_lt(exp(v0$Intercept_ndt), 0.3)
  # and the race parameters start somewhere a race could plausibly be
  expect_equal(log1p(exp(v0$Intercept_sigmabias)), 0.3, tolerance = 1e-8)
  expect_equal(log1p(exp(v0$Intercept_boundary)), 0.5, tolerance = 1e-8)
  expect_equal(log1p(exp(v0$Intercept_driftone)), 3, tolerance = 1e-8)
})


# a real fit --------------------------------------------------------------

test_that("cogmod_rdm recovers ndt above the fastest observed response", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n <- 2000
  sim <- rcogmod_rdm(n, vzero = 2.5, vone = 1.6, boundary = 0.5, bias = 0.2,
                     ndt = 0.25, poutlier = 0.04)
  d <- data.frame(RT = sim$rt, Error = sim$response)
  # the case the old tau * minrt parameterization could not express
  expect_lt(min(d$RT), 0.25)

  f <- brms::bf(RT | dec(Error) ~ 1, driftone ~ 1, sigmabias ~ 1,
                boundary ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_rdm())
  fit <- brms::brm(f, data = d, prior = cogmod_priors(f, d),
                   init = cogmod_inits(f, d), stanvars = cogmod_stanvars(f),
                   algorithm = "pathfinder", refresh = 0,
                   backend = "cmdstanr")

  sp <- function(x) log1p(exp(x))
  s <- brms::posterior_summary(fit)[, "Estimate"]
  expect_equal(sp(s[["b_Intercept"]]), 2.5, tolerance = 0.3)
  expect_equal(sp(s[["b_driftone_Intercept"]]), 1.6, tolerance = 0.3)
  # sigmabias and boundary trade off along a flat ridge; their sum does not
  expect_equal(sp(s[["b_boundary_Intercept"]]) +
                 sp(s[["b_sigmabias_Intercept"]]), 0.7, tolerance = 0.15)
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
