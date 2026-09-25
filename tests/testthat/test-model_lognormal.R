context("Shifted LogNormal - brms")

make_prep <- function(y, mu, sigma, ndt, poutlier, n_draws = 10) {
  structure(
    list(
      data = list(Y = y),
      dpars = list(
        mu = rep(mu, n_draws),
        sigma = rep(sigma, n_draws),
        sigmabias = rep(0, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}

# Reference implementation of the mixture, written out longhand
ref_dens <- function(y, mu, sigma, ndt, poutlier) {
  d_out <- 2 * dnorm(y, 0, 0.2)
  d_dec <- if (y > ndt) dlnorm(y - ndt, mu, sigma) else 0
  poutlier * d_out + (1 - poutlier) * d_dec
}


# dcogmod_lognormal -----------------------------------------------------------

test_that("dcogmod_lognormal matches the mixture density", {
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.02

  for (y in c(0.05, 0.31, 0.5, 0.9, 1.5, 3)) {
    expect_equal(
      dcogmod_lognormal(y, mu, sigma, ndt, poutlier = poutlier),
      ref_dens(y, mu, sigma, ndt, poutlier),
      tolerance = 1e-12,
      label = sprintf("density at y = %.2f", y)
    )
    expect_equal(
      dcogmod_lognormal(y, mu, sigma, ndt, poutlier = poutlier, log = TRUE),
      log(ref_dens(y, mu, sigma, ndt, poutlier)),
      tolerance = 1e-12,
      label = sprintf("log-density at y = %.2f", y)
    )
  }
})

test_that("responses faster than ndt keep positive density (the whole point)", {
  # Without the outlier component these are exactly 0, which is what forces ndt
  # below the fastest observed RT.
  d <- dcogmod_lognormal(0.1, ndt = 0.4, poutlier = 0.02)
  expect_gt(d, 0)
  expect_equal(d, 0.02 * 2 * dnorm(0.1, 0, 0.2), tolerance = 1e-12)

  # With no outlier component they are impossible again
  expect_equal(dcogmod_lognormal(0.1, ndt = 0.4, poutlier = 0), 0)
  expect_true(
    is.infinite(dcogmod_lognormal(0.1, ndt = 0.4, poutlier = 0, log = TRUE))
  )
})

test_that("poutlier = 0 recovers the plain shifted LogNormal", {
  expect_equal(
    dcogmod_lognormal(0.9, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0),
    dlnorm(0.6, -0.7, 0.5),
    tolerance = 1e-12
  )
})

test_that("dcogmod_lognormal is vectorized and integrates to one", {
  d <- dcogmod_lognormal(c(0.5, 0.9), mu = c(-0.7, -0.5), sigma = 0.5, ndt = 0.3)
  expect_length(d, 2)
  expect_equal(
    d,
    c(dlnorm(0.2, -0.7, 0.5), dlnorm(0.6, -0.5, 0.5)),
    tolerance = 1e-12
  )

  total <- integrate(
    function(z) dcogmod_lognormal(z, ndt = 0.3, poutlier = 0.05),
    lower = 0, upper = Inf
  )$value
  expect_equal(total, 1, tolerance = 1e-4)
})

test_that("dcogmod_lognormal returns 0 density for invalid parameters", {
  expect_warning(d <- dcogmod_lognormal(0.5, sigma = -0.5))
  expect_equal(d, 0)
  expect_warning(d <- dcogmod_lognormal(0.5, ndt = -0.3))
  expect_equal(d, 0)
  expect_warning(d <- dcogmod_lognormal(0.5, poutlier = 1.5))
  expect_equal(d, 0)
})


# rcogmod_lognormal -----------------------------------------------------------

test_that("rcogmod_lognormal recovers the mixture mean", {
  set.seed(123)
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  rts <- rcogmod_lognormal(2e4, mu, sigma, ndt, poutlier = poutlier)
  theo <- (1 - poutlier) * (exp(mu + sigma^2 / 2) + ndt) +
    poutlier * (0.2 * sqrt(2 / pi))

  expect_equal(mean(rts), theo, tolerance = 0.05)
  expect_true(all(rts > 0))
  expect_length(rts, 2e4)
})

test_that("rcogmod_lognormal produces responses below ndt only via outliers", {
  set.seed(1)
  expect_true(all(rcogmod_lognormal(5000, ndt = 0.3, poutlier = 0) > 0.3))
  expect_true(any(rcogmod_lognormal(5000, ndt = 0.3, poutlier = 0.5) < 0.3))
})

test_that("rcogmod_lognormal errors on invalid parameters", {
  expect_error(rcogmod_lognormal(10, sigma = -1), "sigma")
  expect_error(rcogmod_lognormal(10, ndt = -1), "ndt")
  expect_error(rcogmod_lognormal(10, poutlier = 2), "poutlier")
  expect_error(rcogmod_lognormal(-5), "non-negative integer")
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_lognormal matches dcogmod_lognormal", {
  for (y in c(0.1, 0.35, 0.9, 2)) {
    prep <- make_prep(y, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02)
    expect_equal(
      log_lik_cogmod_lognormal(1, prep),
      rep(log(ref_dens(y, -0.7, 0.5, 0.3, 0.02)), 10),
      tolerance = 1e-10,
      label = sprintf("log_lik at y = %.2f", y)
    )
  }
})

test_that("log_lik_cogmod_lognormal returns -Inf for invalid parameters", {
  prep <- structure(
    list(
      data = list(Y = 0.5),
      dpars = list(
        mu = c(-0.7, -0.7, -0.7),
        sigma = c(-0.5, 0.5, 0.5),
        sigmabias = c(0, 0, 0),
        ndt = c(0.3, -0.3, 0.3),
        poutlier = c(0.02, 0.02, 1.5)
      )
    ),
    class = "brmsprep"
  )
  expect_true(all(suppressWarnings(log_lik_cogmod_lognormal(1, prep)) == -Inf))
})

test_that("posterior_predict_cogmod_lognormal excludes outliers by default", {
  set.seed(123)
  n_draws <- 1e5
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  prep <- structure(
    list(dpars = list(
      mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
      sigmabias = rep(0, n_draws),
      ndt = rep(ndt, n_draws), poutlier = rep(poutlier, n_draws)
    )),
    class = "brmsprep"
  )

  # No flag on the family: the decision process alone, so nothing below ndt
  rts <- posterior_predict_cogmod_lognormal(1, prep)
  expect_equal(mean(rts), exp(mu + sigma^2 / 2) + ndt, tolerance = 0.02)
  expect_true(all(rts > ndt))

  # Asking for the mixture recovers the mixture mean
  rts_mix <- posterior_predict_cogmod_lognormal(1, prep, predict_outliers = TRUE)
  theo <- (1 - poutlier) * (exp(mu + sigma^2 / 2) + ndt) +
    poutlier * (0.2 * sqrt(2 / pi))

  expect_equal(mean(rts_mix), theo, tolerance = 0.05)
  expect_true(all(rts_mix > 0))
  expect_true(any(rts_mix < ndt))
})

test_that("posterior_epred_cogmod_lognormal excludes outliers by default", {
  mu <- matrix(c(-0.7, -0.5), nrow = 2, ncol = 2)
  sigma <- matrix(0.5, nrow = 2, ncol = 2)
  ndt <- matrix(0.3, nrow = 2, ncol = 2)
  poutlier <- matrix(0.02, nrow = 2, ncol = 2)

  prep <- structure(
    list(dpars = list(mu = mu, sigma = sigma, sigmabias = 0 * mu, ndt = ndt,
                      poutlier = poutlier)),
    class = "brmsprep"
  )

  # Default: the decision-component mean exactly
  expect_equal(
    posterior_epred_cogmod_lognormal(prep),
    exp(mu + sigma^2 / 2) + ndt
  )

  # With the outlier component: the mixture of component means
  expect_equal(
    posterior_epred_cogmod_lognormal(prep, predict_outliers = TRUE),
    0.98 * (exp(mu + sigma^2 / 2) + ndt) + 0.02 * (0.2 * sqrt(2 / pi))
  )

  # and the two differ
  expect_false(isTRUE(all.equal(
    posterior_epred_cogmod_lognormal(prep),
    posterior_epred_cogmod_lognormal(prep, predict_outliers = TRUE)
  )))
})

test_that("the family flag drives predictions when no argument is given", {
  set.seed(123)
  n_draws <- 5e4
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.2 # exaggerated so the two differ unmistakably

  dpars <- list(
    mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
    sigmabias = rep(0, n_draws),
    ndt = rep(ndt, n_draws), poutlier = rep(poutlier, n_draws)
  )
  on <- structure(
    list(dpars = dpars,
         family = list(name = "cogmod_lognormal", predict_outliers = TRUE)),
    class = "brmsprep"
  )
  off <- structure(
    list(dpars = dpars,
         family = list(name = "cogmod_lognormal", predict_outliers = FALSE)),
    class = "brmsprep"
  )

  expect_true(any(posterior_predict_cogmod_lognormal(1, on) < ndt))
  expect_true(all(posterior_predict_cogmod_lognormal(1, off) > ndt))

  # an explicit argument overrides the flag either way
  expect_true(all(
    posterior_predict_cogmod_lognormal(1, on, predict_outliers = FALSE) > ndt
  ))
})

test_that(".predict_outliers resolves the argument, then the family flag", {
  bare <- structure(list(family = list(name = "cogmod_lognormal")),
                    class = "brmsprep")
  on <- structure(list(family = list(name = "cogmod_lognormal",
                                     predict_outliers = TRUE)),
                  class = "brmsprep")

  # no flag anywhere: exclude them
  expect_false(cogmod:::.predict_outliers(NULL, bare))
  # flag on the family is honoured
  expect_true(cogmod:::.predict_outliers(NULL, on))
  # an explicit argument overrides the flag either way
  expect_false(cogmod:::.predict_outliers(FALSE, on))
  expect_true(cogmod:::.predict_outliers(TRUE, bare))
})

test_that("with_outliers sets the flag where prepare_predictions can see it", {
  # brms rebuilds prep$family from object$formula$family, so the flag must live
  # there; object$family alone is dropped. Verified against brms 2.23.1. If this
  # changes upstream, with_outliers() and its docs need revisiting.
  fake <- structure(
    list(family = list(name = "cogmod_lognormal"),
         formula = list(family = list(name = "cogmod_lognormal"))),
    class = "brmsfit"
  )
  m <- with_outliers(fake)
  expect_true(m$formula$family$predict_outliers)
  expect_true(m$family$predict_outliers)

  m2 <- without_outliers(m)
  expect_false(m2$formula$family$predict_outliers)
  expect_false(m2$family$predict_outliers)

  for (f in list(with_outliers, without_outliers)) {
    expect_error(f(list()), "must be a brmsfit")
  }
})

test_that("brms does not forward ... to custom family prediction methods", {
  # This is why the flag rides on the model instead of being an argument.
  flat <- function(f) gsub("\\s+", " ", paste(deparse(body(f)), collapse = " "))

  # posterior_predict.brmsfit sends ... to prepare_predictions(), then calls the
  # brmsprep method with named arguments only
  fit_body <- flat(brms:::posterior_predict.brmsfit)
  expect_true(grepl("prepare_predictions(object", fit_body, fixed = TRUE))
  expect_false(grepl("...", sub(".*posterior_predict\\(prep", "", fit_body),
                     fixed = TRUE))

  # epred and log_lik have no ... at all in their custom dispatch
  expect_false("..." %in% names(formals(brms:::posterior_epred_custom)))
  expect_false("..." %in% names(formals(brms:::log_lik_custom)))
})


# p_outlier ---------------------------------------------------------------

test_that("p_outlier returns the mixture responsibility", {
  y <- c(0.05, 0.4, 0.9, 3.0) # below ndt / just after / bulk / far tail
  nd <- 4
  prep <- list(
    data = list(Y = y),
    dpars = list(
      mu = matrix(-0.9, nd, length(y)), sigma = matrix(0.5, nd, length(y)),
      sigmabias = matrix(0, nd, length(y)),
      ndt = matrix(0.25, nd, length(y)), poutlier = matrix(0.02, nd, length(y))
    )
  )

  out <- with_mocked_bindings(
    p_outlier(structure(list(), class = "brmsfit")),
    prepare_predictions = function(object, ...) prep,
    get_dpar = function(prep, dpar, ...) prep$dpars[[dpar]],
    .package = "brms"
  )

  # Reference responsibility, written out longhand
  ref <- vapply(y, function(v) {
    g <- 2 * dnorm(v, 0, 0.2)
    f <- if (v > 0.25) dlnorm(v - 0.25, -0.9, 0.5) else 0
    0.02 * g / (0.02 * g + 0.98 * f)
  }, numeric(1))

  expect_equal(out$p_outlier, ref, tolerance = 1e-12)
  expect_equal(out$rt, y)
  expect_named(out, c("rt", "p_outlier"))

  # The shape the article plots: certainty below ndt, then collapsing and
  # staying collapsed. Before 0.2.1 the curve rose again in the far slow tail,
  # because the half-t outlier component had heavier tails than the LogNormal
  # and eventually explained slow responses better than the model did. The half
  # Normal cannot, which is why it replaced it.
  expect_equal(out$p_outlier[1], 1)
  expect_lt(out$p_outlier[3], out$p_outlier[2])
  expect_lt(out$p_outlier[4], out$p_outlier[3])
})

test_that("p_outlier stays finite where the natural-scale ratio would be 0/0", {
  # Responsibility is computed in logs. On the natural scale both densities
  # underflow to exactly 0 in the far tails, turning the ratio into 0/0 before
  # either component is genuinely negligible.
  # sigma = 0.05 puts the decision process at 0.62-0.70 s, so 0.66 is the bulk
  # and everything else here is a far tail on one side or the other.
  y <- c(1e-8, 0.05, 0.66, 50, 1e6)
  nd <- 3
  prep <- list(
    data = list(Y = y),
    dpars = list(
      mu = matrix(-0.9, nd, length(y)),
      sigma = matrix(0.05, nd, length(y)), # tight enough to underflow dlnorm
      sigmabias = matrix(0, nd, length(y)),
      ndt = matrix(0.25, nd, length(y)),
      poutlier = matrix(1e-8, nd, length(y)) # tiny rate, as in early warmup
    )
  )

  out <- with_mocked_bindings(
    p_outlier(structure(list(), class = "brmsfit")),
    prepare_predictions = function(object, ...) prep,
    get_dpar = function(prep, dpar, ...) prep$dpars[[dpar]],
    .package = "brms"
  )

  expect_true(all(is.finite(out$p_outlier)))
  expect_true(all(out$p_outlier >= 0 & out$p_outlier <= 1))
  expect_equal(out$p_outlier[1], 1) # below ndt: only the outlier could produce it
  expect_lt(out$p_outlier[3], 1e-3) # in the bulk: the decision process owns it
})

test_that("p_outlier rejects non-brmsfit input", {
  expect_error(p_outlier(list()), "must be a brmsfit")
})


# family ------------------------------------------------------------------

test_that("cogmod_lognormal() builds a valid brms custom family", {
  fam <- cogmod_lognormal()

  expect_s3_class(fam, "customfamily")
  expect_identical(fam$dpars, c("mu", "sigma", "sigmabias", "ndt", "poutlier"))
  expect_identical(
    unname(c(fam$link, fam$link_sigma, fam$link_sigmabias, fam$link_ndt,
             fam$link_poutlier)),
    c("identity", "softplus", "softplus", "log", "logit")
  )
  # the outlier scale is a package constant, never a dpar: brms would estimate it
  expect_false("minrt" %in% fam$dpars)
})

test_that("predict_outliers survives the family validation brm() performs", {
  # The flag has to be an ordinary field on the family and reach
  # object$formula$family untouched, or prepare_predictions() never sees it.
  # Verified against brms 2.23.1.
  expect_false(cogmod_lognormal()$predict_outliers)
  expect_true(cogmod_lognormal(predict_outliers = TRUE)$predict_outliers)

  fam <- cogmod_lognormal(predict_outliers = TRUE)
  expect_true(brms:::validate_family(fam)$predict_outliers)

  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  expect_true(f$family$predict_outliers)

  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  expect_true(brms:::validate_formula(f, data = d)$family$predict_outliers)
})

test_that("stanvars carry the likelihood with the outlier component", {
  sv <- cogmod_lognormal_stanvars()
  code <- paste(vapply(sv, function(x) x$scode, character(1)), collapse = "\n")

  expect_true(grepl("cogmod_lognormal_lpdf", code, fixed = TRUE))
  expect_true(grepl("log_mix", code, fixed = TRUE))
  expect_true(grepl("12.5 * square(Y)", code, fixed = TRUE))
  expect_true(grepl("half Normal with scale 0.2", code, fixed = TRUE))
})


# priors ------------------------------------------------------------------

test_that("brms leaves the ndt and poutlier intercepts flat", {
  # Both are improper, and both directions of each are flat in the likelihood
  # (poutlier -> 1, ndt -> 0), so a model fitted without priors has an improper
  # posterior. The article says so; this pins the upstream behaviour it rests
  # on. Verified against brms 2.23.1.
  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  p <- brms::get_prior(
    brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1),
    family = cogmod_lognormal(), data = d
  )
  flat <- p$prior[p$class == "Intercept" & p$dpar %in% c("ndt", "poutlier")]
  expect_true(all(flat == ""))
})

# The two locations cogmod_priors() sets are pinned here once. Every other
# assertion about them is expressed relative to these, so changing a location in
# the source leaves exactly one test to update rather than four.
.NDT_INTERCEPT_PRIOR <- "normal(-1.2, 0.5)"
.SLOPE_PRIOR <- "normal(0, 0.2)"

.normal_loc <- function(x) as.numeric(sub("^normal\\(([^,]+),.*$", "\\1", x))
.normal_scale <- function(x) sub("^normal\\([^,]+, *(.*)\\)$", "\\1", x)


test_that("cogmod_priors survives arbitrary formula shapes", {
  # The whole point of reading get_prior() rather than guessing: every row comes
  # from the model brms is actually going to build, so a prior matching no
  # parameter is impossible by construction.
  set.seed(1)
  d <- data.frame(
    RT = rcogmod_lognormal(200, ndt = 0.3, poutlier = 0.02),
    x = rnorm(200), g = factor(sample(letters[1:5], 200, TRUE)),
    cond = factor(sample(c("A", "B"), 200, TRUE))
  )
  fam <- cogmod_lognormal()
  forms <- list(
    plain = brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam),
    no_intercept = brms::bf(RT ~ 0 + Intercept + x, sigma ~ 1,
                            ndt ~ 0 + Intercept + x, poutlier ~ 1, family = fam),
    interaction = brms::bf(RT ~ cond * x, sigma ~ 1, ndt ~ cond,
                           poutlier ~ 1 + (1 | g), family = fam),
    ranef_both = brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1 + (1 | g),
                          poutlier ~ 1 + (1 + x | g), family = fam),
    smooth = brms::bf(RT ~ s(x), sigma ~ 1, ndt ~ s(x), poutlier ~ 1,
                      family = fam)
  )

  for (nm in names(forms)) {
    p <- cogmod_priors(forms[[nm]], d)
    expect_s3_class(p, "brmsprior")
    # No ndt or poutlier row is left flat unless a blanket row covers it -
    # brms reports a blanket row plus one per coefficient, and filling the
    # blanket is what makes the per-coefficient ones proper.
    np <- p[p$dpar %in% c("ndt", "poutlier") & !nzchar(p$prior), ]
    covered <- vapply(seq_len(nrow(np)), function(i) {
      any(!nzchar(p$coef) & nzchar(p$prior) & p$class == np$class[i] &
            p$dpar == np$dpar[i] & p$group == np$group[i])
    }, logical(1))
    expect_true(all(covered), label = nm)
    expect_gt(sum(p$dpar %in% c("ndt", "poutlier") & nzchar(p$prior)), 0)
    # brms accepts it, silently: no unmatched prior, no unused blanket row
    expect_silent(
      brms::make_stancode(forms[[nm]], data = d, family = cogmod_lognormal(),
                          prior = p, stanvars = cogmod_lognormal_stanvars())
    )
  }

  # `0 + Intercept` has no Intercept class, so the location has to land on the
  # b coefficient named "Intercept" instead
  p <- cogmod_priors(forms$no_intercept, d)
  b <- p[p$dpar == "ndt" & p$class == "b", ]
  expect_equal(b$prior[b$coef == "Intercept"], .NDT_INTERCEPT_PRIOR)
  expect_equal(b$prior[!nzchar(b$coef)], .SLOPE_PRIOR)
})

test_that("cogmod_priors passes other families straight through", {
  d <- data.frame(y = rnorm(30), x = rnorm(30))

  f <- brms::bf(y ~ x, family = stats::gaussian())
  # the message names the family it found, not "<none found>"
  expect_message(p <- cogmod_priors(f, d), "gaussian")
  expect_message(cogmod_priors(f, d), "cogmod_lognormal")
  # exactly the brms defaults, untouched
  expect_equal(p, brms::validate_prior(brms::empty_prior(), f, d,
                                       family = stats::gaussian()))

  # a formula carrying no family at all is a message, not an error
  expect_message(q <- cogmod_priors(brms::bf(y ~ x), d), "none found")
  expect_s3_class(q, "brmsprior")
})

test_that("returned priors compose and replace with c()", {
  d <- data.frame(RT = rcogmod_lognormal(50, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_lognormal())

  p <- c(
    cogmod_priors(f, d),
    brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
    replace = TRUE
  )
  expect_s3_class(p, "brmsprior")
  code <- brms::make_stancode(f, data = d, family = cogmod_lognormal(), prior = p,
                              stanvars = cogmod_lognormal_stanvars())
  expect_true(grepl("| -2, 0.1", code, fixed = TRUE))
})

test_that("the likelihood is flat in both directions a prior has to cover", {
  y <- rcogmod_lognormal(400, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02)
  ll <- function(nd, p) sum(dcogmod_lognormal(y, -0.7, 0.5, nd, poutlier = p, log = TRUE))

  # poutlier -> 1: the plateau is exactly flat in poutlier itself once it
  # saturates, which is the direction the prior has to cover.
  expect_equal(ll(0.3, plogis(40)), ll(0.3, plogis(60)), tolerance = 1e-8)

  # The other parameters no longer drop out there, though. Under the half-t
  # outlier component of 0.2.0 they did - it could explain any response,
  # however slow - so the density collapsed onto it entirely. A half Normal
  # explains none of the slow ones, so the slowest observations keep a pull on
  # `ndt` even at poutlier = 1 - 1e-12, and the degenerate mode is thousands of
  # log-likelihood units worse rather than hundreds.
  expect_gt(abs(ll(1e-9, 1 - 1e-12) - ll(0.29, 1 - 1e-12)), 1)
  expect_lt(ll(0.3, 1 - 1e-12), ll(0.3, 0.02) - 1000)

  # ndt -> 0 on a log link: flat in log(ndt), which has nothing to do with the
  # mixture and is why a prior on poutlier alone is not enough
  expect_equal(ll(1e-8, 0.02), ll(1e-12, 0.02), tolerance = 1e-6)
})


# omitted dpars -----------------------------------------------------------

test_that("omitting ndt/poutlier from bf() still yields proper priors", {
  # A dpar left out of bf() becomes a plain auxiliary parameter on the natural
  # scale, not a linear predictor on the link scale, so it needs a different
  # prior. brms's own defaults there are uniform(0, min_Y) for ndt - the very
  # min-RT bound this parameterization removes - and flat over [0, 1] for
  # poutlier, which puts half its mass above 0.5.
  set.seed(1)
  d <- data.frame(RT = rcogmod_lognormal(100, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, family = cogmod_lognormal())

  raw <- brms::get_prior(f, data = d, family = cogmod_lognormal())
  expect_match(raw$prior[raw$class == "ndt"], "min_Y")
  expect_equal(raw$prior[raw$class == "poutlier"], "")

  p <- cogmod_priors(f, d)
  expect_equal(p$prior[p$class == "ndt"], "lognormal(-1.2, 0.5)")
  expect_equal(p$prior[p$class == "poutlier"], "exponential(100)")

  code <- brms::make_stancode(f, data = d, family = cogmod_lognormal(), prior = p,
                              stanvars = cogmod_lognormal_stanvars())
  expect_true(grepl("lognormal_lpdf(ndt | -1.2, 0.5)", code, fixed = TRUE))
  expect_true(grepl("exponential_lpdf(poutlier | 100)", code, fixed = TRUE))
  expect_false(grepl("uniform_lpdf(ndt", code, fixed = TRUE))
})

test_that("the natural-scale priors describe the same belief as the link ones", {
  # lognormal(m, s) on ndt is exactly normal(m, s) on log(ndt)
  set.seed(1)
  x <- rlnorm(2e4, -1.2, 0.5)
  expect_equal(mean(log(x)), -1.2, tolerance = 0.01)
  expect_equal(sd(log(x)), 0.5, tolerance = 0.02)

  # exponential(100) keeps the centre of logit-normal(-5, 1) but moves the mode
  # to zero: omitting poutlier from the formula says you expect no outliers.
  # both medians sit at about 0.7%: 0.00693 vs 0.00669, i.e. within 0.03 of a
  # percentage point of each other (a relative check would be the wrong test
  # here - what matters is that neither says "expect a percent or two")
  expect_lt(abs(qexp(0.5, 100) - plogis(-5)), 5e-4)
  expect_lt(qexp(0.95, 100), 0.05)
  # essentially all of it lies inside the [0, 1] support, so truncation is moot
  expect_gt(pexp(1, 100), 1 - 1e-12)
})


# the start-point range ---------------------------------------------------

# For a fixed distance D the decision time is D * exp(mu + sigma Z), i.e.
# LogNormal(mu + log D, sigma), so the density, CDF and survival with a
# start-point range are averages over D ~ Uniform(1, 1 + A) of the plain
# lognormal functions. One-dimensional quadrature of those is the reference for
# the closed-form kernels, series branches included.
ref_acc <- function(t, mu, sigma, A, what = c("dens", "cdf", "surv")) {
  what <- match.arg(what)
  g <- switch(what,
    dens = function(d) stats::dlnorm(t, mu + log(d), sigma),
    cdf = function(d) stats::plnorm(t, mu + log(d), sigma),
    surv = function(d) stats::plnorm(t, mu + log(d), sigma, lower.tail = FALSE)
  )
  stats::integrate(g, 1, 1 + A, rel.tol = 1e-11, abs.tol = 0,
                   subdivisions = 500)$value / A
}

test_that("sigmabias = 0 is the plain shifted LogNormal, bit for bit", {
  t <- c(0.05, 0.3, 1, 4)
  expect_identical(cogmod:::.lognormal_acc_ldens(t, -0.7, 0.6, 0),
                   stats::dlnorm(t, -0.7, 0.6, log = TRUE))
  expect_identical(cogmod:::.lognormal_acc_lcdf(t, -0.7, 0.6, 0),
                   stats::plnorm(t, -0.7, 0.6, log.p = TRUE))
  expect_identical(cogmod:::.lognormal_acc_lccdf(t, -0.7, 0.6, 0),
                   stats::plnorm(t, -0.7, 0.6, lower.tail = FALSE,
                                 log.p = TRUE))
  # the general kernel is continuous with it through the series branch
  expect_equal(cogmod:::.lognormal_acc_ldens(t, -0.7, 0.6, 1e-9),
               stats::dlnorm(t, -0.7, 0.6, log = TRUE), tolerance = 1e-8)
  expect_equal(cogmod:::.lognormal_acc_lcdf(t, -0.7, 0.6, 1e-9),
               stats::plnorm(t, -0.7, 0.6, log.p = TRUE), tolerance = 1e-8)
  # the default argument is that zero, in all three user-facing functions
  expect_identical(dcogmod_lognormal(0.5), dcogmod_lognormal(0.5, sigmabias = 0))
  expect_identical(pcogmod_lognormal(0.5), pcogmod_lognormal(0.5, sigmabias = 0))
  set.seed(5); a <- rcogmod_lognormal(20)
  set.seed(5); b <- rcogmod_lognormal(20, sigmabias = 0)
  expect_identical(a, b)
})

test_that("the start-point kernels agree with quadrature over the start point", {
  grid <- covering_grid(
    t = c(0.02, 0.1, 0.4, 1, 2.5, 8, 30),
    mu = c(-1.5, -0.2, 1),
    sigma = c(0.3, 0.8, 1.5),
    # 1e-7 and 5e-5 take the series branch, the rest the general kernel
    A = c(1e-7, 5e-5, 1e-3, 0.05, 0.5, 3, 20),
    # both sides of the series switch, at the far ends of the time axis
    always = function(g) {
      (g$t == 30 & g$mu == 1 & g$sigma == 0.3 & g$A == 1e-7) |
        (g$t == 0.02 & g$mu == -1.5 & g$sigma == 0.3 & g$A == 1e-3) |
        (g$t == 30 & g$mu == 1 & g$sigma == 1.5 & g$A == 20)
    }
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    info <- sprintf("t = %g, mu = %g, sigma = %g, A = %g", g$t, g$mu, g$sigma,
                    g$A)
    ld <- cogmod:::.lognormal_acc_ldens(g$t, g$mu, g$sigma, g$A)
    lF <- cogmod:::.lognormal_acc_lcdf(g$t, g$mu, g$sigma, g$A)
    lS <- cogmod:::.lognormal_acc_lccdf(g$t, g$mu, g$sigma, g$A)
    expect_false(is.nan(ld) || is.nan(lF) || is.nan(lS), info = info)
    expect_lte(lF, 0)
    expect_lte(lS, 0)
    rd <- ref_acc(g$t, g$mu, g$sigma, g$A, "dens")
    rF <- ref_acc(g$t, g$mu, g$sigma, g$A, "cdf")
    rS <- ref_acc(g$t, g$mu, g$sigma, g$A, "surv")
    if (rd > 1e-250) expect_equal(exp(ld), rd, tolerance = 1e-7, info = info)
    if (rF > 1e-250) expect_equal(exp(lF), rF, tolerance = 1e-7, info = info)
    if (rS > 1e-250) expect_equal(exp(lS), rS, tolerance = 1e-7, info = info)
    # and the two tails are complementary in their own right
    if (is.finite(lF) && is.finite(lS)) {
      expect_equal(exp(lF) + exp(lS), 1, tolerance = 1e-9, info = info)
    }
  }
})

test_that("with a start-point range the mixture still integrates to one", {
  for (A in c(0.5, 3)) {
    for (pout in c(0, 0.05)) {
      total <- stats::integrate(
        function(t) dcogmod_lognormal(t, -0.7, 0.5, ndt = 0.2, poutlier = pout,
                                      sigmabias = A),
        lower = 0, upper = Inf, subdivisions = 2000
      )$value
      expect_equal(total, 1, tolerance = 1e-5,
                   info = sprintf("sigmabias %g, poutlier %g", A, pout))
    }
  }
})

test_that("pcogmod_lognormal integrates dcogmod_lognormal with a start-point range", {
  for (q in c(0.25, 0.5, 1, 2.5)) {
    num <- stats::integrate(
      function(t) dcogmod_lognormal(t, -0.7, 0.5, ndt = 0.2, poutlier = 0.03,
                                    sigmabias = 1),
      0, q, rel.tol = 1e-10
    )$value
    expect_equal(pcogmod_lognormal(q, -0.7, 0.5, ndt = 0.2, poutlier = 0.03,
                                   sigmabias = 1), num, tolerance = 1e-7)
  }
})

test_that("the mean and the sampler carry the start-point range", {
  set.seed(21)
  n <- 200000
  x <- rcogmod_lognormal(n, mu = -0.7, sigma = 0.5, ndt = 0.2, sigmabias = 1)
  # E[T] = ndt + (1 + A / 2) exp(mu + sigma^2 / 2)
  expect_equal(mean(x), 0.2 + 1.5 * exp(-0.7 + 0.125), tolerance = 0.005)
  # and the sampler's quantiles match the CDF
  for (p in c(0.1, 0.5, 0.9)) {
    q <- unname(stats::quantile(x, p))
    expect_equal(pcogmod_lognormal(q, -0.7, 0.5, ndt = 0.2, sigmabias = 1), p,
                 tolerance = 0.01)
  }
  # posterior_epred() gets the same factor through the registry's mean
  spec <- cogmod:::.shifted_spec("cogmod_lognormal")
  expect_equal(spec$mean(list(mu = -0.7, sigma = 0.5, sigmabias = 1)),
               1.5 * exp(-0.7 + 0.125))
})

test_that("a negative start-point range is rejected", {
  expect_warning(d <- dcogmod_lognormal(0.5, sigmabias = -0.1))
  expect_equal(d, 0)
  expect_error(rcogmod_lognormal(5, sigmabias = -0.1), "sigmabias")
})

test_that("Stan cogmod_lognormal_lpdf, _lcdf and _lccdf match R with a start-point range", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  lpdf <- stan_fun("cogmod_lognormal")
  lcdf <- stan_fun("cogmod_lognormal", "_lcdf")
  lccdf <- stan_fun("cogmod_lognormal", "_lccdf")
  grid <- covering_grid(
    Y = c(0.21, 0.25, 0.4, 0.8, 1.5, 6),
    mu = c(-1.2, -0.5, 0.3),
    sigma = c(0.3, 0.7, 1.2),
    # 0 takes the plain lognormal branch, 1e-6 the series, the rest the general
    # kernel on both sides of a = 0
    sigmabias = c(0, 1e-6, 0.3, 2),
    poutlier = c(0, 0.03),
    always = function(g) g$sigmabias %in% c(0, 1e-6) & g$Y %in% c(0.21, 6)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    info <- sprintf("Y = %g, mu = %g, sigma = %g, sigmabias = %g, poutlier = %g",
                    g$Y, g$mu, g$sigma, g$sigmabias, g$poutlier)
    r <- dcogmod_lognormal(g$Y, g$mu, g$sigma, 0.2, g$sigmabias, g$poutlier,
                           log = TRUE)
    s <- lpdf(g$Y, g$mu, g$sigma, g$sigmabias, 0.2, g$poutlier)
    expect_lt(abs(s - r) / max(1, abs(r)), 1e-10, label = info)
    rF <- pcogmod_lognormal(g$Y, g$mu, g$sigma, 0.2, g$sigmabias, g$poutlier,
                            log.p = TRUE)
    rS <- pcogmod_lognormal(g$Y, g$mu, g$sigma, 0.2, g$sigmabias, g$poutlier,
                            lower.tail = FALSE, log.p = TRUE)
    sF <- lcdf(g$Y, g$mu, g$sigma, g$sigmabias, 0.2, g$poutlier)
    sS <- lccdf(g$Y, g$mu, g$sigma, g$sigmabias, 0.2, g$poutlier)
    expect_lt(abs(sF - rF) / max(1, abs(rF)), 1e-9, label = paste(info, "lcdf"))
    expect_lt(abs(sS - rS) / max(1, abs(rS)), 1e-9, label = paste(info, "lccdf"))
  }
  expect_equal(lpdf(0.5, -0.7, 0.5, -0.1, 0.2, 0.02), -Inf)
})

test_that("cogmod_priors fences the start-point range for cogmod_lognormal", {
  set.seed(22)
  d <- data.frame(RT = rcogmod_lognormal(150, ndt = 0.25, poutlier = 0.03),
                  Condition = rep(c("a", "b"), length.out = 150))
  modelled <- brms::bf(RT ~ 1, sigmabias ~ Condition, family = cogmod_lognormal())
  p <- cogmod_priors(modelled, d)
  expect_true(any(p$dpar == "sigmabias" & p$class == "Intercept" &
                    p$prior == "normal(0, 1)"))
  expect_true(any(p$dpar == "sigmabias" & p$class == "b" &
                    p$prior == "normal(0, 0.5)"))
  omitted <- brms::bf(RT ~ 1, family = cogmod_lognormal())
  p2 <- cogmod_priors(omitted, d)
  expect_true(any(p2$class == "sigmabias" &
                    p2$prior == "lognormal(-0.35, 0.75)"))
  # pinned in the formula it is not a parameter, so there is nothing to fence
  pinned <- brms::bf(RT ~ 1, sigmabias = 0, family = cogmod_lognormal())
  p3 <- cogmod_priors(pinned, d)
  expect_false(any(p3$dpar == "sigmabias" | p3$class == "sigmabias"))
  code <- brms::make_stancode(pinned, data = d, prior = p3,
                              stanvars = cogmod_stanvars(pinned))
  expect_true(grepl("real sigmabias = 0;", code, fixed = TRUE))
  # and cens() builds with the range modelled (half the trials censored, which
  # draws the heavy-censoring warning: not what this test is about)
  d$cens <- rep(c(0L, 1L), length.out = 150)
  censored <- brms::bf(RT | cens(cens) ~ 1, sigmabias ~ 1,
                       family = cogmod_lognormal())
  code2 <- brms::make_stancode(censored, data = d,
                               prior = suppressWarnings(cogmod_priors(censored, d)),
                               stanvars = cogmod_stanvars(censored))
  expect_true(grepl("cogmod_lognormal_lccdf", code2, fixed = TRUE))
})
