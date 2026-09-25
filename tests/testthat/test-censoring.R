context("Censoring: cens() on the RT-only families")

# What brms does with `bf(rt | cens(x) ~ ...)` on a custom family is call
# `<family>_lccdf()` / `<family>_lcdf()` in the Stan program and the family's own
# log_lik() method in R - and it does the censoring itself in neither place
# beyond that. So three things have to agree: the R p*() function, the generated
# Stan CDF pair, and the log_lik() branch. The Stan checks go through
# stan_fun(), the one shared compilation.
#
# The two places bmm's cswald broke - `t == ndt` and the far upper tail - are
# tested explicitly throughout.

# Every censorable family in the shifted registry, with the defaults of its own
# r*() function as the decision parameters. Derived, not listed, so a family
# given a `stan_lccdf` slot is covered here without anyone remembering to add it.
cens_shifted <- intersect(cogmod:::.CENS_FAMILIES, names(cogmod:::.SHIFTED))
dec_defaults <- function(nm) {
  fm <- formals(get(paste0("r", nm)))
  fm <- fm[setdiff(names(fm), c("n", "ndt", "poutlier"))]
  lapply(fm, eval)
}

# Equality that is relative where the numbers are ordinary and absolute where a
# log-probability sits next to zero (a CDF within rounding of 1), where a
# relative criterion measures nothing but the two implementations' rounding.
expect_close <- function(a, b, tol = 1e-8, info = NULL) {
  if (is.infinite(a) && is.infinite(b) && sign(a) == sign(b)) {
    return(expect_true(TRUE))
  }
  expect_true(is.finite(a) && is.finite(b) && abs(a - b) <= tol * max(1, abs(b)),
              info = paste0(info, ": ", format(a, digits = 12), " vs ",
                            format(b, digits = 12)))
}


test_that("the censorable set is derived from the registry and excludes the races", {
  expect_true(all(c("cogmod_invgaussian", "cogmod_lognormal", "cogmod_bisa",
                    "cogmod_exgaussian", "cogmod_geg") %in% cogmod:::.CENS_FAMILIES))
  # No choice family: the outlier's 1 / K keeps their defective densities
  # summing to one, and a bare survival term would break that.
  expect_false(any(cogmod:::.CHOICE_FAMILIES %in% cogmod:::.CENS_FAMILIES))
  # No closed-form CDF, no censoring - deliberately, until one is written.
  expect_false(any(c("cogmod_loggamma", "cogmod_exwald", "cogmod_lba1") %in%
                     cogmod:::.CENS_FAMILIES))
  for (nm in cens_shifted) {
    expect_true(exists(paste0("p", nm), mode = "function"), info = nm)
    dp <- cogmod:::.SHIFTED[[nm]]$dpars
    expect_equal(names(dec_defaults(nm)) |> length(), length(dp), info = nm)
  }
})


test_that("p*() is a CDF: bounded, monotone, complementary, and the integral of d*()", {
  q <- c(0, 0.05, 0.2, 0.3, 0.5, 1, 2, 3, Inf)
  for (nm in cens_shifted) {
    pf <- get(paste0("p", nm))
    df <- get(paste0("d", nm))
    for (pout in c(0, 0.05)) {
      args <- c(dec_defaults(nm), list(ndt = 0.2, poutlier = pout))
      lo <- do.call(pf, c(list(q), args))
      up <- do.call(pf, c(list(q), args, list(lower.tail = FALSE)))
      expect_true(all(diff(lo) >= -1e-12), info = nm)
      expect_equal(lo + up, rep(1, length(q)), tolerance = 1e-10, info = nm)
      expect_equal(lo[1], 0, info = nm)
      expect_equal(lo[length(q)], 1, info = nm)
      # Below the shift only the outlier component has any mass...
      expect_equal(do.call(pf, c(list(0.15), args)),
                   pout * cogmod:::.pcontam(0.15), tolerance = 1e-12, info = nm)
      # ...and exactly at it too: the decision process has not started.
      expect_equal(do.call(pf, c(list(0.2), args)),
                   pout * cogmod:::.pcontam(0.2), tolerance = 1e-12, info = nm)
      expect_equal(do.call(pf, c(list(q[2:8]), args, list(log.p = TRUE))),
                   log(lo[2:8]), tolerance = 1e-12, info = nm)

      # The CDF is the integral of the density. Split at the shift, where the
      # integrand is discontinuous. The log-Student-t's density is UNBOUNDED
      # there - integrable, but integrate() declares it divergent - so that one
      # is integrated from just above the shift, at the cost of the sliver it
      # leaves out, and held to a looser tolerance.
      spike <- nm == "cogmod_logstudent"
      tol <- if (spike) 1e-3 else 1e-6
      from <- if (spike) 0.2 + 1e-6 else 0.2
      dens <- function(x) do.call(df, c(list(x), args))
      for (qq in c(0.3, 0.6)) {
        num <- stats::integrate(dens, 0, 0.2, rel.tol = 1e-10)$value +
          stats::integrate(dens, from, qq, rel.tol = 1e-10,
                           subdivisions = 500L)$value
        expect_equal(do.call(pf, c(list(qq), args)), num, tolerance = tol,
                     info = paste(nm, qq))
      }
    }
    # NA in, NA out, as for every p*() in base R.
    expect_true(is.na(do.call(pf, c(list(NA_real_), dec_defaults(nm),
                                    list(ndt = 0.2, poutlier = 0)))))
  }
})


test_that("the survival keeps its digits at the shift and in the far tail", {
  # With poutlier > 0 the outlier survival is 2 * Phi(-t / 0.2): exp(-112) at
  # 3 s, a number rather than -Inf. Written as log1p(-Phi(t / 0.2)) it would
  # have been -Inf from about 1.66 s on - the first of bmm's two bug sites.
  for (nm in cens_shifted) {
    pf <- get(paste0("p", nm))
    args <- c(dec_defaults(nm), list(ndt = 0.2, poutlier = 0.05))
    at_ndt <- do.call(pf, c(list(0.2), args, list(lower.tail = FALSE, log.p = TRUE)))
    expect_equal(at_ndt, log(0.05 * 2 * stats::pnorm(-1) + 0.95),
                 tolerance = 1e-12, info = nm)
    far <- do.call(pf, c(list(c(3, 10)), args, list(lower.tail = FALSE, log.p = TRUE)))
    expect_true(all(is.finite(far)), info = nm)
    # Decreasing, and below the outlier component alone - which at 3 s is
    # already exp(-112) and has nothing to do with how far down these sit: the
    # heavy-tailed families (LogNormal, inverse Gamma...) keep a survival of
    # order 1e-3 to 1e-6 at 10 s, the Wald and the Gamma of order 1e-20.
    expect_true(diff(far) < 0, info = nm)
    expect_true(far[1] < log(0.5), info = nm)
  }
  # And with no outlier component at all the Wald's own survival, far out, is
  # the difference of two terms that 1 - CDF would have returned as 0.
  ls <- pcogmod_invgaussian(10, drift = 3, boundary = 0.5, ndt = 0.2,
                            poutlier = 0, lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(ls) && ls < -40)
  expect_equal(pcogmod_invgaussian(10, drift = 3, boundary = 0.5, ndt = 0.2,
                                   poutlier = 0, lower.tail = FALSE), 0)
})


test_that("the Wald CDF and survival agree with each other and with .pwald_raw() under drift variability", {
  t <- c(0.05, 0.1, 0.3, 0.6, 1)
  for (sv in c(0.3, 1, 2)) {
    lS <- cogmod:::.lswald_raw(t, 3, 0.5, sv)
    lF <- cogmod:::.lpwald_raw(t, 3, 0.5, sv)
    expect_equal(exp(lS) + exp(lF), rep(1, 5), tolerance = 1e-10, info = sv)
    expect_equal(exp(lF), cogmod:::.pwald_raw(t, 3, 0.5, sv), tolerance = 1e-10,
                 info = sv)
  }
  # Continuous into sigmadrift = 0.
  expect_equal(cogmod:::.lswald_raw(0.3, 3, 0.5, 1e-6),
               cogmod:::.lswald_raw(0.3, 3, 0.5, 0), tolerance = 1e-5)
  expect_equal(cogmod:::.lpwald_raw(0.3, 3, 0.5, 1e-6),
               cogmod:::.lpwald_raw(0.3, 3, 0.5, 0), tolerance = 1e-5)
})


test_that("the ex-Gaussian survival is the complement of its CDF", {
  x <- c(-0.3, 0, 0.3, 0.6, 1, 2, 4)
  lo <- pcogmod_exgaussian(x, mu = 0.4, sigma = 0.1, tau = 0.2)
  up <- pcogmod_exgaussian(x, mu = 0.4, sigma = 0.1, tau = 0.2, lower.tail = FALSE)
  expect_equal(lo + up, rep(1, length(x)), tolerance = 1e-10)
  num <- vapply(x, function(q) {
    stats::integrate(function(v) dcogmod_exgaussian(v, 0.4, 0.1, 0.2), -Inf, q,
                     rel.tol = 1e-10)$value
  }, numeric(1))
  expect_equal(lo, num, tolerance = 1e-7)
  # The far tail, as a log survival, is the two positive terms and stays finite.
  ls <- pcogmod_exgaussian(8, 0.4, 0.1, 0.2, lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(ls) && ls < -30)
  expect_equal(pcogmod_exgaussian(8, 0.4, 0.1, 0.2, lower.tail = FALSE), 0)
})


# Stan -------------------------------------------------------------------

test_that("Stan <family>_lcdf and _lccdf match the R p*() functions", {
  skip_if_not_installed("cmdstanr")
  # Y == ndt and Y just above it are the first bmm bug site; 3 s the second.
  grid <- covering_grid(
    y = c(0.1, 0.2, 0.2001, 0.35, 0.8, 3),
    ndt = c(0, 0.2),
    pout = c(0, 0.05, 0.4),
    always = function(g) g$y %in% c(0.2, 3)
  )
  for (nm in cens_shifted) {
    lcdf <- stan_fun(nm, "_lcdf")
    lccdf <- stan_fun(nm, "_lccdf")
    pf <- get(paste0("p", nm))
    dec <- dec_defaults(nm)
    for (r in seq_len(nrow(grid))) {
      y <- grid$y[r]; ndt <- grid$ndt[r]; pout <- grid$pout[r]
      info <- sprintf("%s: Y=%.4f ndt=%.1f poutlier=%.2f", nm, y, ndt, pout)
      rF <- do.call(pf, c(list(y), dec, list(ndt = ndt, poutlier = pout,
                                              log.p = TRUE)))
      rS <- do.call(pf, c(list(y), dec, list(ndt = ndt, poutlier = pout,
                                              lower.tail = FALSE, log.p = TRUE)))
      sF <- do.call(lcdf, c(list(y), unname(dec), list(ndt, pout)))
      sS <- do.call(lccdf, c(list(y), unname(dec), list(ndt, pout)))
      expect_close(sF, rF, info = paste(info, "lcdf"))
      expect_close(sS, rS, info = paste(info, "lccdf"))
      # The two Stan functions are complementary in their own right.
      if (is.finite(sF) && is.finite(sS)) {
        expect_equal(exp(sF) + exp(sS), 1, tolerance = 1e-8, info = info)
      }
    }
    # Invalid parameters are refused the way the lpdf refuses them.
    bad <- do.call(lccdf, c(list(0.5), unname(dec), list(-0.1, 0.02)))
    expect_equal(bad, -Inf, info = nm)
  }
})


test_that("Stan and R take the drift-variability Wald CDF over the same nodes", {
  skip_if_not_installed("cmdstanr")
  lcdf <- stan_fun("cogmod_invgaussian", "_lcdf")
  lccdf <- stan_fun("cogmod_invgaussian", "_lccdf")
  grid <- covering_grid(
    y = c(0.25, 0.4, 0.8, 2, 4),
    sv = c(0.05, 0.5, 2),
    pout = c(0, 0.05)
  )
  for (r in seq_len(nrow(grid))) {
    y <- grid$y[r]; sv <- grid$sv[r]; pout <- grid$pout[r]
    info <- sprintf("Y=%.2f sigmadrift=%.2f poutlier=%.2f", y, sv, pout)
    rF <- pcogmod_invgaussian(y, 3, 0.5, 0.2, sv, poutlier = pout, log.p = TRUE)
    rS <- pcogmod_invgaussian(y, 3, 0.5, 0.2, sv, poutlier = pout,
                              lower.tail = FALSE, log.p = TRUE)
    expect_close(lcdf(y, 3, 0.5, sv, 0, 0.2, pout), rF, info = paste(info, "lcdf"))
    expect_close(lccdf(y, 3, 0.5, sv, 0, 0.2, pout), rS, info = paste(info, "lccdf"))
  }
})


test_that("Stan and R agree on the Wald CDF and survival with a non-decision range", {
  skip_if_not_installed("cmdstanr")
  # sigmandt spreads the shift over [ndt, ndt + sigmandt]. The closed forms
  # (cogmod_wald_st0 and .lwald_st0_fixed()) branch three ways - interval
  # reaching back past the shift, CDF-side difference, survival-side
  # difference - and the smallest range below sits on the midpoint branch. Each
  # is visited, with and without the drift quadrature around it.
  lpdf <- stan_fun("cogmod_invgaussian")
  lcdf <- stan_fun("cogmod_invgaussian", "_lcdf")
  lccdf <- stan_fun("cogmod_invgaussian", "_lccdf")
  grid <- covering_grid(
    y = c(0.2, 0.205, 0.23, 0.3, 0.5, 1, 3),
    sv = c(0, 0.5),
    st0 = c(1e-7, 0.02, 0.15),
    pout = c(0, 0.05),
    always = function(g) g$y %in% c(0.2, 0.205, 3)
  )
  for (r in seq_len(nrow(grid))) {
    y <- grid$y[r]; sv <- grid$sv[r]; st0 <- grid$st0[r]; pout <- grid$pout[r]
    info <- sprintf("Y=%.3f sigmadrift=%.2f sigmandt=%g poutlier=%.2f",
                    y, sv, st0, pout)
    rD <- dcogmod_invgaussian(y, 3, 0.5, 0.2, sv, st0, pout, log = TRUE)
    rF <- pcogmod_invgaussian(y, 3, 0.5, 0.2, sv, st0, pout, log.p = TRUE)
    rS <- pcogmod_invgaussian(y, 3, 0.5, 0.2, sv, st0, pout,
                              lower.tail = FALSE, log.p = TRUE)
    expect_close(lpdf(y, 3, 0.5, sv, st0, 0.2, pout), rD, info = paste(info, "lpdf"))
    expect_close(lcdf(y, 3, 0.5, sv, st0, 0.2, pout), rF, info = paste(info, "lcdf"))
    expect_close(lccdf(y, 3, 0.5, sv, st0, 0.2, pout), rS, info = paste(info, "lccdf"))
  }
})


test_that("Stan cogmod_exgaussian and cogmod_geg CDFs match their R counterparts", {
  skip_if_not_installed("cmdstanr")
  ex_lcdf <- stan_fun("cogmod_exgaussian", "_lcdf")
  ex_lccdf <- stan_fun("cogmod_exgaussian", "_lccdf")
  geg_lcdf <- stan_fun("cogmod_geg", "_lcdf")
  geg_lccdf <- stan_fun("cogmod_geg", "_lccdf")
  for (y in c(-0.2, 0.3, 0.5, 1, 3, 6)) {
    # y = -0.2 is six SDs into the LEFT tail, where the two terms of F_EG are
    # near-equal numbers of order 1e-9. Stan's exp_mod_normal_lcdf() subtracts
    # them before taking the log and kept about seven digits there; the Stan
    # side now does the subtraction in log space, the same way
    # .lcdf_exgaussian() does in R, so the one criterion holds everywhere.
    ltol <- 1e-8
    expect_close(ex_lcdf(y, 0.4, 0.1, 0.2),
                 pcogmod_exgaussian(y, 0.4, 0.1, 0.2, log.p = TRUE),
                 tol = ltol, info = paste("exgaussian lcdf", y))
    expect_close(ex_lccdf(y, 0.4, 0.1, 0.2),
                 pcogmod_exgaussian(y, 0.4, 0.1, 0.2, lower.tail = FALSE, log.p = TRUE),
                 info = paste("exgaussian lccdf", y))
    for (shape in c(0.5, 1, 3)) {
      expect_close(geg_lcdf(y, 0.4, 0.1, 0.2, shape),
                   pcogmod_geg(y, 0.4, 0.1, 0.2, shape, log.p = TRUE),
                   tol = ltol, info = paste("geg lcdf", y, shape))
      # log1m_exp of a CDF within rounding of 1 is where the two sides' rounding
      # shows, so the far tail is held to an absolute criterion only.
      expect_close(geg_lccdf(y, 0.4, 0.1, 0.2, shape),
                   pcogmod_geg(y, 0.4, 0.1, 0.2, shape, lower.tail = FALSE, log.p = TRUE),
                   tol = if (y >= 3) 1e-3 else 1e-8,
                   info = paste("geg lccdf", y, shape))
    }
  }
  # shape = 1 is the ex-Gaussian exactly, in the CDF as in the density.
  expect_equal(geg_lccdf(0.7, 0.4, 0.1, 0.2, 1), ex_lccdf(0.7, 0.4, 0.1, 0.2),
               tolerance = 1e-10)
})


# log_lik ----------------------------------------------------------------

# A brmsprep by hand: brms::get_dpar() only needs the class and a draws x
# observations matrix per dpar, and building one this way is what lets the
# branch be tested without a fit.
fake_prep <- function(family, dpars, Y, cens = NULL, rcens = NULL) {
  nd <- nrow(dpars[[1]])
  structure(list(
    ndraws = nd, nobs = length(Y),
    data = Filter(Negate(is.null), list(Y = Y, cens = cens, rcens = rcens)),
    dpars = dpars,
    family = family
  ), class = "brmsprep")
}


test_that("log_lik() honours cens(): brms leaves that to the family", {
  set.seed(3)
  nd <- 5
  Y <- c(0.5, 0.7, 0.9, 0.4)
  mu <- matrix(stats::rnorm(nd * 4, -0.7, 0.05), nd)
  prep <- fake_prep(
    cogmod_lognormal(),
    list(mu = mu, sigma = matrix(0.5, nd, 4), sigmabias = matrix(0, nd, 4),
         ndt = matrix(0.2, nd, 4), poutlier = matrix(0.02, nd, 4)),
    Y, cens = c(0L, 1L, -1L, 2L), rcens = c(NA, NA, NA, 0.8)
  )
  # observed
  expect_equal(log_lik_cogmod_lognormal(1, prep),
               dcogmod_lognormal(0.5, mu[, 1], 0.5, 0.2, poutlier = 0.02, log = TRUE))
  # right-censored: the survival
  expect_equal(log_lik_cogmod_lognormal(2, prep),
               pcogmod_lognormal(0.7, mu[, 2], 0.5, 0.2, poutlier = 0.02, lower.tail = FALSE,
                                 log.p = TRUE))
  # left-censored: the CDF
  expect_equal(log_lik_cogmod_lognormal(3, prep),
               pcogmod_lognormal(0.9, mu[, 3], 0.5, 0.2, poutlier = 0.02, log.p = TRUE))
  # interval-censored between Y and rcens
  expect_equal(log_lik_cogmod_lognormal(4, prep),
               log(pcogmod_lognormal(0.8, mu[, 4], 0.5, 0.2, poutlier = 0.02) -
                     pcogmod_lognormal(0.4, mu[, 4], 0.5, 0.2, poutlier = 0.02)))
  # No cens() on the formula: the density, as before.
  prep$data$cens <- NULL
  expect_equal(log_lik_cogmod_lognormal(2, prep),
               dcogmod_lognormal(0.7, mu[, 2], 0.5, 0.2, poutlier = 0.02, log = TRUE))

  # The two families outside the registry go through the same branch.
  ex <- fake_prep(
    cogmod_exgaussian(),
    list(mu = matrix(0.4, nd, 2), sigma = matrix(0.1, nd, 2), tau = matrix(0.2, nd, 2)),
    c(0.5, 0.7), cens = c(0L, 1L)
  )
  expect_equal(log_lik_cogmod_exgaussian(2, ex),
               rep(pcogmod_exgaussian(0.7, 0.4, 0.1, 0.2, lower.tail = FALSE,
                                      log.p = TRUE), nd))
  expect_equal(log_lik_cogmod_exgaussian(1, ex),
               rep(dcogmod_exgaussian(0.5, 0.4, 0.1, 0.2, log = TRUE), nd))
  geg <- fake_prep(
    cogmod_geg(),
    list(mu = matrix(0.4, nd, 2), sigma = matrix(0.1, nd, 2),
         tau = matrix(0.2, nd, 2), shape = matrix(2, nd, 2)),
    c(0.5, 0.7), cens = c(0L, 1L)
  )
  expect_equal(log_lik_cogmod_geg(2, geg),
               rep(pcogmod_geg(0.7, 0.4, 0.1, 0.2, 2, lower.tail = FALSE,
                               log.p = TRUE), nd))
})


# brms and the data check -------------------------------------------------

set.seed(5)
d_cens <- data.frame(
  RT = rcogmod_invgaussian(100, ndt = 0.2),
  err = rep(c(0L, 1L), c(90, 10)),
  Condition = rep(c("a", "b"), 50)
)
f_cens <- brms::bf(RT | cens(err) ~ Condition, boundary ~ 1, ndt ~ 1,
                   sigmadrift = 0, family = cogmod_invgaussian())


test_that("brms writes the censored likelihood against the functions the stanvars define", {
  sv <- cogmod_stanvars(f_cens)
  scode <- paste(vapply(sv, function(s) s$scode, character(1)), collapse = "\n")
  expect_match(scode, "real cogmod_invgaussian_lccdf(", fixed = TRUE)
  expect_match(scode, "real cogmod_invgaussian_lcdf(", fixed = TRUE)
  expect_match(scode, "real cogmod_invgaussian_lpdf(", fixed = TRUE)

  code <- suppressWarnings(brms::make_stancode(
    f_cens, data = d_cens, prior = cogmod_priors(f_cens, d_cens), stanvars = sv
  ))
  expect_match(code, "cogmod_invgaussian_lccdf(Y[n] |", fixed = TRUE)
  expect_match(code, "cogmod_invgaussian_lcdf(Y[n] |", fixed = TRUE)
  # And the other two generics are unbothered by the addition term.
  expect_no_error(cogmod_inits(f_cens, d_cens)(1))
})


test_that("cens() is refused where it cannot go and questioned where it is doubtful", {
  # A choice family already models the errors, through dec().
  expect_error(
    cogmod:::.cogmod_checkdata(
      brms::bf(RT | dec(err) + cens(err) ~ 1, family = cogmod_lnr()), d_cens),
    "already models the errors"
  )
  # No closed-form CDF, no survival to score a censored trial with.
  f_lg <- brms::bf(RT | cens(err) ~ 1, family = cogmod_loggamma())
  expect_error(cogmod:::.cogmod_checkdata(f_lg, d_cens), "no closed-form CDF")
  expect_error(cogmod_stanvars(f_lg), "no closed-form CDF")
  expect_error(cogmod_priors(f_lg, d_cens), "no closed-form CDF")

  expect_silent(cogmod:::.cogmod_checkdata(f_cens, d_cens))
  # brms's three encodings are all read.
  expect_silent(cogmod:::.cogmod_checkdata(
    f_cens, transform(d_cens, err = ifelse(err == 1, "right", "none"))))
  expect_silent(cogmod:::.cogmod_checkdata(f_cens, transform(d_cens, err = err == 1)))
  expect_error(cogmod:::.cogmod_checkdata(f_cens, transform(d_cens, err = 3)),
               "censoring indicator")
  # A fifth of the trials censored is past what the model is argued for.
  expect_warning(
    cogmod:::.cogmod_checkdata(f_cens, transform(d_cens, err = rep(c(0L, 1L), 50))),
    "right-censored"
  )
})


test_that("a censored shifted Wald recovers its parameters and loo() runs", {
  skip_if_not_slow()
  skip_if_not_installed("loo")
  set.seed(11)
  n <- 600
  # Correct responses from the model; an independent competing process that
  # gets in first on a minority of trials, which is exactly the non-informative
  # censoring the likelihood assumes.
  latent <- rcogmod_invgaussian(n, drift = 3, boundary = 0.6, ndt = 0.25)
  competitor <- stats::runif(n, 0.35, 2.5)
  d <- data.frame(RT = pmin(latent, competitor),
                  err = as.integer(competitor < latent))
  f <- brms::bf(RT | cens(err) ~ 1, boundary ~ 1, ndt ~ 1, sigmadrift = 0,
                family = cogmod_invgaussian())
  m <- brms::brm(f, data = d, prior = cogmod_priors(f, d),
                 init = cogmod_inits(f, d), stanvars = cogmod_stanvars(f),
                 chains = 2, iter = 1000, refresh = 0, backend = "cmdstanr")
  s <- brms::posterior_summary(m)
  softplus <- function(x) log1p(exp(x))
  expect_equal(softplus(s["b_Intercept", "Estimate"]), 3, tolerance = 0.15)
  expect_equal(softplus(s["b_boundary_Intercept", "Estimate"]), 0.6, tolerance = 0.15)
  expect_equal(exp(s["b_ndt_Intercept", "Estimate"]), 0.25, tolerance = 0.1)
  ll <- brms::log_lik(m)
  expect_true(all(is.finite(ll)))
  # The censored rows are scored with the survival, which is what makes them
  # LARGER than the density of the same time would be for the bulk of them.
  expect_s3_class(suppressWarnings(brms::loo(m)), "loo")
})
