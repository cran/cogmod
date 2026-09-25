# Shared machinery for the shifted RT families
# =============================================
#
# Every family in this group has the same structure: a decision-time
# distribution shifted by `ndt`, mixed with a half Normal outlier component
# of weight `poutlier` and a fixed scale. Only the decision distribution
# differs. Keeping the mixture in one place is what stops the nine families
# drifting apart - the Stan code, the R density, the RNG, the likelihood, the
# predictions, the parameter checks and the priors are all generated from the
# one registry below, so a family cannot end up with, say, an outlier component
# in Stan but not in `log_lik()`, or an R check that rejects what Stan accepts.
#
# See ?rcogmod_lognormal for the account of why `ndt` is expressed directly rather
# than as a fraction of an observed minimum, what the outlier component is for,
# and why the outlier component's scale is a constant rather than a dpar.


# The registry ------------------------------------------------------------

# Each entry describes the *decision* component only.
#
#  dpars/links/lb/ub : the distributional parameters, before ndt and poutlier
#  stan_check        : Stan expression that is TRUE for invalid parameters
#  stan_dens         : Stan expression for the decision log-density at `t_adj`
#  ldens             : R log-density of the decision component, at t > 0
#  stan_lcdf         : optional, Stan expression for the decision log-CDF at
#                      `t_adj` > 0
#  stan_lccdf        : optional, the log-survival, same footing. A family with
#                      both gets a `<name>_lcdf` and a `<name>_lccdf` beside its
#                      `<name>_lpdf`, which is what brms's `cens()` addition
#                      term calls - so these two slots are what makes a family
#                      censorable. Write the survival directly, never as
#                      log(1 - exp(lcdf)): the far tail is exactly where a
#                      censored slow error lands, and that form cancels there
#  lcdf / lccdf      : the R counterparts, at t > 0, as functions of (t, p)
#  rng               : R sampler for the decision component
#  mean              : E[decision time]; Inf where it does not exist
#  prior             : optional, per-dpar priors for cogmod_priors() to fill,
#                      where the likelihood has a flat direction brms would
#                      otherwise leave improper. Each entry is a named character
#                      vector with any of `link` (modelled intercept, on the
#                      dpar's link scale), `nat` (dpar omitted from bf(), so on
#                      the natural scale) and `slope`.
#  label             : human-readable name, used in the generated Stan comments
#' @keywords internal
.SHIFTED <- list(
  cogmod_lognormal = list(
    # `sigmabias` is a start-point range: the decision time is a LogNormal
    # multiplied by a Uniform(1, 1 + sigmabias) distance, which is the
    # single-accumulator LBA with a LogNormal rate and its threshold offset
    # pinned at 1 (sigma cannot pin the scale of a LogNormal rate). At
    # sigmabias = 0 the family is the shifted LogNormal exactly as it was
    # before the parameter existed, and cogmod_lnr() is a race of two of these
    # accumulators. The derivation and the numerics live with the kernels
    # below (.lognormal_acc_ldens() and friends).
    dpars = c("mu", "sigma", "sigmabias"),
    links = c("identity", "softplus", "softplus"),
    lb = c(NA, 0, 0), ub = c(NA, NA, NA),
    # A start-point range of exactly zero is a model, not an invalid parameter
    # - the plain LogNormal - so its bound is closed, as in cogmod_lba1() and
    # cogmod_lnr(). The kernels take the plain lognormal branch there, so the
    # nested model is reached exactly and at its own cost.
    lb_open = c(TRUE, TRUE, FALSE),
    stan_check = "sigma <= 0 || sigmabias < 0",
    stan_dens = "cogmod_lognormal_acc_ldens(t_adj, mu, sigma, sigmabias)",
    ldens = function(t, p) .lognormal_acc_ldens(t, p$mu, p$sigma, p$sigmabias),
    stan_lcdf = "cogmod_lognormal_acc_logcdf(t_adj, mu, sigma, sigmabias)",
    stan_lccdf = "cogmod_lognormal_acc_logsurv(t_adj, mu, sigma, sigmabias)",
    lcdf = function(t, p) .lognormal_acc_lcdf(t, p$mu, p$sigma, p$sigmabias),
    lccdf = function(t, p) .lognormal_acc_lccdf(t, p$mu, p$sigma, p$sigmabias),
    prelude = ".LOGNORMAL_STAN_PRELUDE",
    # The LogNormal is drawn first, so at sigmabias = 0 the stream is exactly
    # what rcogmod_lognormal() produced before the parameter existed.
    rng = function(n, p) {
      r <- stats::rlnorm(n, p$mu, p$sigma)
      (1 + stats::runif(n, min = 0, max = p$sigmabias)) * r
    },
    # E[D] E[exp(mu + sigma Z)], the distance and the rate being independent
    mean = function(p) (1 + p$sigmabias / 2) * exp(p$mu + p$sigma^2 / 2),
    init = list(mu = -0.7, sigma = 0.5, sigmabias = 0.5),
    # `sigmabias` has the flat direction cogmod_lba1()'s has: as the range
    # approaches zero the model converges to the LogNormal and the likelihood
    # stops changing, while the softplus link reaches zero only at minus
    # infinity. The rows are cogmod_lnr()'s verbatim - same unit (the threshold
    # offset of 1), same job - so that a fitted LogNormal can warm-start one
    # lane of the race. Fixing `sigmabias = 0` in bf() sidesteps all of it, and
    # is the recommendation.
    prior = list(
      sigmabias = c(link = "normal(0, 1)", nat = "lognormal(-0.35, 0.75)",
                    slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "mu: mean of the decision time on the log scale (meanlog).",
      "sigma: SD of the decision time on the log scale (> 0).",
      paste("sigmabias: start-point range, in units of the threshold offset",
            "(>= 0); 0 is the plain LogNormal.")
    ),
    label = "LogNormal"
  ),
  cogmod_logstudent = list(
    # log(decision time) ~ scaled Student-t, i.e. a robust LogNormal. `dof` is
    # the degrees of freedom - what brms::student() calls `nu`, renamed because
    # cogmod_lnr() already spends `nuzero`/`nuone` on drift rates, and because
    # brms recognises the name `nu` and injects defaults of its own for it
    # (normal(2.7, 0.8) modelled, gamma(2, 0.1) omitted). `dof` arrives flat like
    # every other dpar this package defines, so cogmod_priors() simply fills it.
    dpars = c("mu", "sigma", "dof"),
    links = c("identity", "softplus", "log"),
    lb = c(NA, 0, 0), ub = c(NA, NA, NA),
    stan_check = "sigma <= 0 || dof <= 0",
    stan_dens = "student_t_lpdf(log(t_adj) | dof, mu, sigma) - log(t_adj)",
    ldens = function(t, p) {
      stats::dt((log(t) - p$mu) / p$sigma, df = p$dof, log = TRUE) -
        log(p$sigma) - log(t)
    },
    stan_lcdf = "student_t_lcdf(log(t_adj) | dof, mu, sigma)",
    stan_lccdf = "student_t_lccdf(log(t_adj) | dof, mu, sigma)",
    lcdf = function(t, p) {
      stats::pt((log(t) - p$mu) / p$sigma, df = p$dof, log.p = TRUE)
    },
    lccdf = function(t, p) {
      stats::pt((log(t) - p$mu) / p$sigma, df = p$dof, lower.tail = FALSE,
                log.p = TRUE)
    },
    rng = function(n, p) exp(p$mu + p$sigma * stats::rt(n, p$dof)),
    # E[exp(sigma * T)] with T a Student-t diverges for EVERY finite dof: the t
    # has polynomial tails and exp() outruns them, so the moment generating
    # function does not exist anywhere. Unlike cogmod_logweibull(), where the
    # mean exists below sigma = 1, there is no region of the parameter space
    # where this one has an expectation, so posterior_epred() has nothing to
    # return at all - see posterior_epred_cogmod_logstudent(). The median is
    # exact: ndt + exp(mu).
    mean = NULL,
    init = list(mu = -0.7, sigma = 0.4, dof = 5),
    # Two flat-ish directions, and the prior is what fences both.
    #
    # Upward: the LogNormal is only reached as dof -> Inf, and it is approached
    # slowly - at dof = 100 the density still differs from the LogNormal by 66%
    # in the tail, and it takes dof ~ 1e4 to get within 1%. So above about
    # dof = 30 the likelihood has almost stopped moving, and a log link puts
    # that region at plus infinity.
    #
    # Downward: a t is symmetric ON THE LOG SCALE, so a small dof does not only
    # buy the slow tail this family exists for - it also piles mass just above
    # `ndt`, which is exactly what poutlier is for. The two then compete. See
    # the `note` below for the numbers.
    #
    # normal(1.8, 0.7) on the log link is centred at dof = 6.0 with 95% between
    # 1.5 and 24: heavy enough for the slow tail to be worth having, light
    # enough that the fast-side spike stays under a tenth of a percent.
    prior = list(
      dof = c(link = "normal(1.8, 0.7)", nat = "lognormal(1.8, 0.7)",
              slope = "normal(0, 0.3)")
    ),
    dpar_doc = c(
      "mu: location of the decision time on the log scale.",
      "sigma: scale of the decision time on the log scale (> 0).",
      paste("dof: degrees of freedom of the Student-t on the log-RT scale",
            "(> 0). Smaller is"),
      "   heavier-tailed; dof -> Inf is the LogNormal."
    ),
    note = c(
      paste("The decision density is UNBOUNDED at `ndt`: as t_adj -> 0 it grows",
            "like 1 / (t_adj *"),
      paste("|log t_adj|^(dof + 1)), where a LogNormal decays to zero. The spike",
            "is integrable for"),
      paste("every dof > 0, so the posterior stays proper, but the likelihood",
            "has no maximum -"),
      "it is the prior on `ndt` that keeps the sampler off min(RT).",
      "",
      paste("A Student-t is symmetric on the log scale, so a small `dof` fattens",
            "BOTH tails. At"),
      paste("dof = 2, 1.5% of the decision distribution falls below 0.05 s",
            "against a LogNormal's"),
      paste("5e-9, which is territory `poutlier` also claims; at dof = 5 it is",
            "0.1%, while the"),
      paste("slow tail is still five orders of magnitude heavier than a",
            "LogNormal's at 5 s."),
      "Expect `dof` and `poutlier` to trade off if the prior on either is widened."
    ),
    label = "Log-Student-t"
  ),
  cogmod_loggamma = list(
    dpars = c("mu", "sigma", "shape"),
    links = c("identity", "softplus", "identity"),
    lb = c(NA, 0, NA), ub = c(NA, NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = paste0(
      "-log(sigma) - log(t_adj) + cogmod_loggamma_lconst(shape)",
      " + cogmod_loggamma_lkernel(shape, (log(t_adj) - mu) / sigma)"
    ),
    prelude = ".LOGGAMMA_STAN_PRELUDE",
    ldens = function(t, p) .dloggamma(t, p$mu, p$sigma, p$shape, log = TRUE),
    rng = function(n, p) .rloggamma(n, p$mu, p$sigma, p$shape),
    mean = function(p) .mloggamma(p$mu, p$sigma, p$shape),
    init = list(mu = -0.7, sigma = 0.5, shape = 0),
    dpar_doc = c(
      "mu: location of the decision time on the log scale.",
      "sigma: scale of the decision time on the log scale (> 0).",
      paste("shape: shape of the log-gamma on the log-RT scale",
            "(unconstrained). shape = 0 is the"),
      "   LogNormal, shape = sigma the Gamma, shape = 1 the Weibull."
    ),
    note = c(
      paste("Note that when sigma * shape >= 1 the decision density is",
            "unbounded at `ndt`, so"),
      paste("the likelihood is unbounded too; the prior on shape is what",
            "keeps the sampler out"),
      "of that region. See ?rcogmod_loggamma."
    ),
    label = "Log-Gamma"
  ),
  cogmod_invgaussian = list(
    # mu is the drift rate and boundary the decision threshold, so the underlying
    # inverse Gaussian has mean boundary / mu and shape boundary^2.
    #
    # `sigmadrift` is the across-trial SD of the drift: with sigmadrift > 0 the
    # drift is drawn once per trial from a Normal(mu, sigmadrift) *truncated at
    # zero*, and the decision time is the Wald first passage given that draw.
    # Marginalising over the draw is a Gaussian integral, so the density stays
    # closed form - see .dwald_raw(). Truncating is what cogmod_lba1() does with
    # the same quantity, and it is needed here for the same reason: a
    # single-boundary accumulator with a negative drift never terminates, so an
    # untruncated Normal would leave the density defective (only 0.99 of its
    # mass at mu = 3, boundary = 0.5, sigmadrift = 1.5, and 0.69 at mu = 0.5,
    # boundary = 1, sigmadrift = 2). cogmod_ddm()'s `sigmadrift` needs no
    # truncation because a diffusion between two boundaries always absorbs at
    # one of them.
    #
    # `sigmandt` is the between-trial RANGE of the non-decision time, st0 in
    # the usual notation: each trial's non-decision time is drawn from
    # Uniform(ndt, ndt + sigmandt), so `ndt` is the lower bound, exactly as for
    # cogmod_ddm()'s `sigmandt` and rtdists's `t0` / `st0`. Smearing the shift
    # over an interval turns the density into a difference of two CDFs and the
    # CDF into a difference of two integrated CDFs, and for the Wald both are
    # closed form at a fixed drift - see .lwald_st0_fixed(). The decision
    # component still has no mass below `ndt`, which is why the shared mixture
    # template needs no change to carry it.
    #
    # It is on a `log` link where the other three are `softplus`: it is the same
    # quantity as cogmod_ddm()'s `sigmandt`, in the same unit, and it lives at
    # 0.01-0.2 s, where log and softplus agree to within a few percent anyway.
    # What log buys there is a prior that means the same thing whether or not
    # the parameter is written in bf(): a normal on the link scale IS a
    # lognormal on the natural one, so the two rows of the prior table below
    # describe one distribution. Softplus has no such counterpart, and its
    # advantage - an additive rather than multiplicative scale for large values
    # - never arises for a range of tens of milliseconds.
    dpars = c("mu", "boundary", "sigmadrift", "sigmandt"),
    links = c("softplus", "softplus", "softplus", "log"),
    lb = c(0, 0, 0, 0), ub = c(NA, NA, NA, NA),
    # A zero drift SD is the plain Wald and a zero range a fixed non-decision
    # time, so those bounds are closed - unlike the drift and the threshold,
    # which a Wald needs strictly positive.
    lb_open = c(TRUE, TRUE, FALSE, FALSE),
    stan_check = "mu <= 0 || boundary <= 0 || sigmadrift < 0 || sigmandt < 0",
    stan_dens = "cogmod_invgaussian_decision_lpdf(t_adj | mu, boundary, sigmadrift, sigmandt)",
    prelude = ".WALD_STAN_PRELUDE",
    ldens = function(t, p) {
      .dwald_raw(t, p$mu, p$boundary, p$sigmadrift, p$sigmandt)
    },
    # Closed form at sigmadrift = 0, with or without sigmandt; 64-point
    # Gauss-Legendre over the drift above it, in Stan and in R alike, off the
    # same nodes - see .lwald_sv().
    stan_lcdf = "cogmod_invgaussian_decision_lcdf(t_adj | mu, boundary, sigmadrift, sigmandt)",
    stan_lccdf = "cogmod_invgaussian_decision_lccdf(t_adj | mu, boundary, sigmadrift, sigmandt)",
    lcdf = function(t, p) {
      .lpwald_raw(t, p$mu, p$boundary, p$sigmadrift, p$sigmandt)
    },
    lccdf = function(t, p) {
      .lswald_raw(t, p$mu, p$boundary, p$sigmadrift, p$sigmandt)
    },
    rng = function(n, p) {
      .rwald_raw(n, p$mu, p$boundary, p$sigmadrift, p$sigmandt)
    },
    # E[T] = boundary / mu + sigmandt / 2 holds only while the drift is fixed.
    # Once it varies the density decays as t^-2 - drifts arbitrarily close to
    # zero take arbitrarily long - so the mean diverges, exactly as it does for
    # cogmod_lba1(). Inf is what posterior_epred() then returns.
    mean = function(p) {
      ifelse(p$sigmadrift > 0, Inf, p$boundary / p$mu + p$sigmandt / 2)
    },
    init = list(mu = 3, boundary = 0.5, sigmadrift = 0.5, sigmandt = 0.05),
    # Same flat direction as cogmod_lba1()'s `sigmabias` and cogmod_ddm()'s own
    # `sigmadrift`: the floor at zero is only reached by the softplus link at
    # minus infinity, and the likelihood stops changing well before then, which
    # is an improper posterior wherever brms would leave it flat. There is a
    # second reason here - scaling (mu, boundary, sigmadrift) up by a common
    # factor sends the Wald to the reciprocal-normal (LATER) limit, so large
    # values of all three describe very nearly the same distribution. The prior
    # is what keeps the sampler out of both.
    #
    # sigmandt takes cogmod_ddm()'s prior for the same parameter, verbatim:
    # the same flat direction (a log link reaches zero only at minus infinity),
    # the same unit, and the same reputation as the variability a model is
    # least able to recover.
    prior = list(
      sigmadrift = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                     slope = "normal(0, 0.5)"),
      sigmandt = c(link = "normal(-3, 1)", nat = "lognormal(-3, 1)",
                   slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "mu: drift rate, the average speed of evidence accumulation (> 0).",
      "boundary: decision threshold, the evidence needed to respond (> 0).",
      paste("sigmadrift: between-trial SD of the drift rate (>= 0), which is",
            "drawn from a"),
      "   Normal(mu, sigmadrift) truncated at zero. 0 is the plain Wald.",
      paste("sigmandt: between-trial range of the non-decision time (st0,",
            ">= 0), in the same"),
      "   unit as Y, with `ndt` its lower bound. 0 is a fixed non-decision time."
    ),
    note = c(
      paste("sigmadrift and poutlier both fatten the right tail, and they are",
            "only weakly"),
      paste("distinguishable: on 2000 simulated trials at mu = 3, boundary =",
            "0.5, sigmadrift ="),
      paste("0.8, estimating sigmadrift buys about 2 log-likelihood units over",
            "fixing it at"),
      "zero. Fix it (`sigmadrift = 0` in bf()) unless the design can identify it.",
      paste("sigmandt is harder still: it shares the leading edge with ndt and",
            "poutlier, and"),
      paste("should be fixed (`sigmandt = 0` in bf()) unless there is a lot of",
            "data, a strong"),
      "prior, or both."
    ),
    label = "Wald (inverse Gaussian)"
  ),
  cogmod_exwald = list(
    # Schwarz (2001): the decision time is a Wald convolved with an Exponential
    # of mean `tau` - a diffusive decision stage followed by an exponentially
    # distributed residual stage. The mechanistic counterpart of
    # cogmod_exgaussian(), whose first stage is a descriptive Gaussian instead,
    # and `tau` means the same thing in both.
    #
    # No `sigmadrift` here, deliberately. It and `tau` both fatten the right
    # tail and would be near-impossible to tell apart; cogmod_invgaussian() is
    # where the drift-variability route lives.
    dpars = c("mu", "boundary", "tau"),
    links = c("softplus", "softplus", "softplus"),
    lb = c(0, 0, 0), ub = c(NA, NA, NA),
    stan_check = "mu <= 0 || boundary <= 0 || tau <= 0",
    stan_dens = "cogmod_exwald_decision_lpdf(t_adj | mu, boundary, tau)",
    prelude = ".EXWALD_STAN_PRELUDE",
    ldens = function(t, p) .dexwald_raw(t, p$mu, p$boundary, p$tau),
    rng = function(n, p) {
      .rwald_raw(n, p$mu, p$boundary) + stats::rexp(n, 1 / p$tau)
    },
    # Finite, unlike the drift-variability Wald: the exponential stage has a
    # mean and the fixed-drift Wald has one too, so they simply add.
    mean = function(p) p$boundary / p$mu + p$tau,
    init = list(mu = 3, boundary = 0.5, tau = 0.15),
    # `tau` is a length of time in seconds, and brms leaves it flat. Same belief,
    # and the same numbers, as the `tau` of cogmod_exgaussian(): it is the same
    # quantity, the mean of an exponential residual stage. See cogmod_priors().
    prior = list(
      tau = c(link = "normal(-1.5, 0.7)", nat = "lognormal(-1.5, 0.7)",
              slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "mu: drift rate, the average speed of evidence accumulation (> 0).",
      "boundary: decision threshold, the evidence needed to respond (> 0).",
      paste("tau: mean of the exponential residual stage, in seconds (> 0).",
            "The decision time is")
    ),
    note = c(
      paste("`ndt` and `tau` both delay the response and are separated only by",
            "shape - `ndt` is a"),
      paste("hard floor, `tau` a variable stage - so they sit on a ridge. The",
            "prior on `ndt` is"),
      "what holds it; widen it and expect the two to trade off.",
      "",
      paste("The density has two branches, both exact. Where mu^2 > 2 / tau the",
            "convolution collapses"),
      paste("to a closed form in the Wald CDF; below that the same expression",
            "continues"),
      paste("analytically through the Faddeeva function. They meet exactly at",
            "mu^2 = 2 / tau.")
    ),
    label = "ex-Wald"
  ),
  cogmod_bisa = list(
    # Birnbaum-Saunders (fatigue life), carried in the same (drift, threshold)
    # currency as cogmod_invgaussian() - and that is the point of it. The Wald
    # accumulates by diffusion, so evidence can move either way at any instant;
    # here it arrives in discrete cycles and only ever TOWARDS the boundary.
    # What is random is the SIZE of each increment, never its sign.
    #
    # Sum n such increments and apply the CLT: S_n ~ Normal(n mu, n), the
    # increment SD being fixed at 1 per cycle - the same convention that fixes
    # the Wald's diffusion coefficient, which is what keeps `mu` and `boundary`
    # comparable across the two families. Then
    #
    #   P(T <= n) = P(S_n >= boundary) = Phi((mu n - boundary) / sqrt(n)),
    #
    # and treating the cycle count n as continuous is the first-crossing
    # approximation. So the standardised time is
    #
    #   z = (mu t - boundary) / sqrt(t)   ~   Normal(0, 1)   exactly,
    #
    # which is the usual (1 / a)(sqrt(t / b) - sqrt(b / t)) written in these
    # parameters: b = boundary / mu is the scale and a = 1 / sqrt(mu * boundary)
    # the shape. That map is a bijection - boundary = sqrt(b) / a and
    # mu = 1 / (a sqrt(b)) - so the mechanistic parameterization gives up
    # nothing, and cogmod_invgaussian(sigmadrift = 0) then differs from this
    # family ONLY in how the evidence arrives.
    #
    # Note there is no third parameter and there cannot be one: the shape is
    # pinned by mu * boundary once the increment SD is fixed, exactly as the
    # Wald's own shape is pinned at boundary^2.
    dpars = c("mu", "boundary"),
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || boundary <= 0",
    stan_dens = "cogmod_bisa_decision_lpdf(t_adj | mu, boundary)",
    prelude = ".BISA_STAN_PRELUDE",
    ldens = function(t, p) .dbisa_raw(t, p$mu, p$boundary),
    # (mu t - boundary) / sqrt(t) is exactly standard normal, so the CDF is one
    # Phi and the survival is the same Phi of the negated argument - no
    # log(1 - .) anywhere.
    stan_lcdf = "std_normal_lcdf((mu * t_adj - boundary) / sqrt(t_adj))",
    stan_lccdf = "std_normal_lcdf((boundary - mu * t_adj) / sqrt(t_adj))",
    lcdf = function(t, p) {
      stats::pnorm((p$mu * t - p$boundary) / sqrt(t), log.p = TRUE)
    },
    lccdf = function(t, p) {
      stats::pnorm((p$boundary - p$mu * t) / sqrt(t), log.p = TRUE)
    },
    rng = function(n, p) .rbisa_raw(n, p$mu, p$boundary),
    # The Wald's mean plus 1 / (2 mu^2), and always finite - nothing here varies
    # across trials, so posterior_epred() always has a number to return. The
    # variance is boundary / mu^3 + 5 / (4 mu^4), against the Wald's
    # boundary / mu^3.
    mean = function(p) p$boundary / p$mu + 1 / (2 * p$mu^2),
    # The same starting point as cogmod_invgaussian(), the parameters meaning
    # the same thing there.
    init = list(mu = 3, boundary = 0.5),
    dpar_doc = c(
      paste("mu: drift rate, the average size of the per-cycle evidence",
            "increment (> 0)."),
      "boundary: decision threshold, the evidence needed to respond (> 0)."
    ),
    note = c(
      paste("The density is the Wald's, tilted: f_BS(t) = f_Wald(t) *",
            "(mu t + boundary) / (2 boundary)."),
      paste("Equivalently it is an EQUAL MIXTURE of the Wald with the same",
            "parameters and that"),
      paste("Wald's length-biased version, so at a given (mu, boundary) it is",
            "both slower and more"),
      paste("spread out: at mu = 3, boundary = 0.5 the mean is 0.222 s against",
            "the Wald's 0.167 and"),
      paste("the SD 0.184 against 0.136. The right tail is still",
            "exponential-order, exp(-mu^2 t / 2),"),
      "like the Wald's and unlike a LogNormal's.",
      "",
      paste("The median is exactly boundary / mu, which is the Wald's MEAN.",
            "Do not read the two"),
      "families' parameters as describing the same central tendency."
    ),
    label = "Birnbaum-Saunders (fatigue life)"
  ),
  cogmod_gamma = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "gamma_lpdf(t_adj | mu, inv(sigma))",
    ldens = function(t, p) {
      stats::dgamma(t, shape = p$mu, scale = p$sigma, log = TRUE)
    },
    stan_lcdf = "gamma_lcdf(t_adj | mu, inv(sigma))",
    stan_lccdf = "gamma_lccdf(t_adj | mu, inv(sigma))",
    lcdf = function(t, p) {
      stats::pgamma(t, shape = p$mu, scale = p$sigma, log.p = TRUE)
    },
    lccdf = function(t, p) {
      stats::pgamma(t, shape = p$mu, scale = p$sigma, lower.tail = FALSE,
                    log.p = TRUE)
    },
    rng = function(n, p) stats::rgamma(n, shape = p$mu, scale = p$sigma),
    mean = function(p) p$mu * p$sigma,
    init = list(mu = 2, sigma = 0.2),
    label = "Gamma"
  ),
  cogmod_invgamma = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "inv_gamma_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      p$mu * log(p$sigma) - lgamma(p$mu) - (p$mu + 1) * log(t) - p$sigma / t
    },
    stan_lcdf = "inv_gamma_lcdf(t_adj | mu, sigma)",
    stan_lccdf = "inv_gamma_lccdf(t_adj | mu, sigma)",
    # T = 1 / G with G ~ Gamma(shape mu, rate sigma), so P(T <= t) is the UPPER
    # tail of the Gamma at 1 / t, and the survival its lower tail.
    lcdf = function(t, p) {
      stats::pgamma(1 / t, shape = p$mu, rate = p$sigma, lower.tail = FALSE,
                    log.p = TRUE)
    },
    lccdf = function(t, p) {
      stats::pgamma(1 / t, shape = p$mu, rate = p$sigma, log.p = TRUE)
    },
    rng = function(n, p) 1 / stats::rgamma(n, shape = p$mu, rate = p$sigma),
    # The mean of an inverse Gamma exists only for shape > 1.
    mean = function(p) ifelse(p$mu > 1, p$sigma / (p$mu - 1), Inf),
    init = list(mu = 4, sigma = 1.5),
    label = "inverse Gamma"
  ),
  cogmod_weibull = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "weibull_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      stats::dweibull(t, shape = p$mu, scale = p$sigma, log = TRUE)
    },
    stan_lcdf = "weibull_lcdf(t_adj | mu, sigma)",
    stan_lccdf = "weibull_lccdf(t_adj | mu, sigma)",
    lcdf = function(t, p) {
      stats::pweibull(t, shape = p$mu, scale = p$sigma, log.p = TRUE)
    },
    lccdf = function(t, p) {
      stats::pweibull(t, shape = p$mu, scale = p$sigma, lower.tail = FALSE,
                      log.p = TRUE)
    },
    rng = function(n, p) stats::rweibull(n, shape = p$mu, scale = p$sigma),
    mean = function(p) p$sigma * gamma(1 + 1 / p$mu),
    init = list(mu = 2, sigma = 0.5),
    label = "Weibull"
  ),
  cogmod_invweibull = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "frechet_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      log(p$mu) - log(p$sigma) - (1 + p$mu) * (log(t) - log(p$sigma)) -
        (t / p$sigma)^(-p$mu)
    },
    stan_lcdf = "frechet_lcdf(t_adj | mu, sigma)",
    stan_lccdf = "frechet_lccdf(t_adj | mu, sigma)",
    # log F = -(t / sigma)^-mu exactly; the survival is its log1m_exp, which is
    # the stable direction here because log F -> 0 only as t -> Inf.
    lcdf = function(t, p) -(t / p$sigma)^(-p$mu),
    lccdf = function(t, p) .log1m_exp(-(t / p$sigma)^(-p$mu)),
    rng = function(n, p) p$sigma * (-log(stats::runif(n)))^(-1 / p$mu),
    # E[Frechet] is finite only for shape > 1.
    mean = function(p) ifelse(p$mu > 1, p$sigma * gamma(1 - 1 / p$mu), Inf),
    init = list(mu = 3, sigma = 0.4),
    label = "inverse Weibull (Frechet)"
  ),
  cogmod_lba1 = list(
    # Single-accumulator LBA: mu is the drift rate, sigma its across-trial SD
    # (fixed to 1 by convention - the evidence scale is arbitrary), sigmabias
    # the start-point range A, and boundary the threshold offset, so b = sigmabias + boundary.
    dpars = c("mu", "sigma", "sigmabias", "boundary"),
    links = c("softplus", "softplus", "softplus", "softplus"),
    lb = c(NA, 0, 0, 0), ub = c(NA, NA, NA, NA),
    # A start-point range of exactly zero is a model, not an invalid parameter:
    # the accumulator then starts at 0 every trial, the finishing time is b / v,
    # and 1 / (RT - ndt) is normal - the recinormal, i.e. the LATER model of
    # Carpenter & Williams (1995). Its bound is therefore closed, for the same
    # reason cogmod_invgaussian() closes `sigmadrift` and cogmod_ddm() closes
    # its three between-trial variabilities: the nested simpler model is
    # reachable by pinning the extension parameter at its floor.
    #
    # Nothing else has to change for that to be exact. `delta = sigmabias / st`
    # is then 0, which takes the Taylor branch of .lba_dens_over_A(), and there
    # the series is phi(z1) * (drift + sigma * z1) / st with drift + sigma * z1
    # = b / t identically - that is phi(z1) * b / (sigma * t^2), the recinormal
    # density itself rather than an approximation to it. .lba_lsurv_trunc()
    # likewise collapses to the survival of a deterministic start point, which
    # is what cogmod_lba2() needs.
    lb_open = c(TRUE, TRUE, FALSE, TRUE),
    stan_check = "sigma <= 0 || sigmabias < 0 || boundary <= 0",
    stan_dens = "cogmod_lba1_decision_lpdf(t_adj | mu, sigma, sigmabias, boundary)",
    prelude = ".LBA_STAN_PRELUDE",
    ldens = function(t, p) .dlba1_raw(t, p$mu, p$sigma, p$sigmabias, p$boundary),
    rng = function(n, p) .rlba1_raw(n, p$mu, p$sigma, p$sigmabias, p$boundary),
    # E[(b - U(0, A)) / v] with v a normal truncated at 0 diverges: the drift
    # density is positive at v = 0, so E[1/v] does not exist. posterior_epred()
    # therefore has nothing to return - see posterior_epred_cogmod_lba1().
    mean = NULL,
    init = list(mu = 3, sigma = 1, sigmabias = 0.5, boundary = 0.5),
    # Every one of these four is a length on the evidence scale, and that scale
    # has no unit: multiply them all by any c > 0 and (b - z) / v is unchanged,
    # so the likelihood is EXACTLY constant along that ray. Fixing any one of
    # them in bf() at a NON-ZERO value pins it; zero is the one value the
    # scaling leaves alone, so `sigmabias = 0` does not. .warn_scale_ray() says
    # so when nothing pins it.
    scale_ray = c("mu", "sigma", "sigmabias", "boundary"),
    # `sigmabias` is estimable but treacherous, and the two are worth keeping
    # apart. Estimating it, the likelihood becomes *flat* as it approaches zero
    # - the LBA is converging to the recinormal, and once the start-point range
    # is small enough making it smaller stops changing the density - while a
    # softplus link reaches zero only at minus infinity. Flat prior plus
    # infinite flat region is the improper posterior cogmod_priors() exists to
    # prevent; without these rows `sigmabias` runs off to -11 with Rhat near 1.7.
    # `boundary` is fenced for the same reason: b = sigmabias + boundary, so the
    # two share the ridge. See ?rcogmod_lba1.
    #
    # Pinning `sigmabias = 0` in bf() is the other option and has none of that
    # trouble: it is not a parameter, so there is no flat direction to sample
    # along and cogmod_priors() emits no row for it. That is the recinormal, and
    # it is the honest thing to do when the design cannot identify a start-point
    # range rather than letting a prior invent one.
    #
    # `link` is the intercept prior when the dpar is modelled in bf() (softplus
    # scale here), `nat` the one for a dpar omitted from bf(), which brms
    # declares as a plain auxiliary parameter on the natural scale instead.
    # Both say the same thing: a start-point range and a threshold offset of
    # order 0.5, and neither of them zero.
    #
    # `slope` widens the blanket normal(0, 0.2) the other dpars get. A condition
    # difference in the start-point range is a real effect of ordinary size -
    # 0.74 on the softplus scale between speed and accuracy instructions on the
    # data in vignette("rt_models") - which normal(0, 0.2) would shrink most of
    # the way to zero. The job here is to fence off the flat direction at
    # sigmabias = 0, not to shrink effects.
    prior = list(
      sigmabias = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                    slope = "normal(0, 0.5)"),
      boundary = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                   slope = "normal(0, 0.5)")
    ),
    label = "single-accumulator LBA"
  ),
  cogmod_logweibull = list(
    # log(decision time) ~ Gumbel(mu, sigma), i.e. the log-Weibull.
    dpars = c("mu", "sigma"),
    links = c("identity", "softplus"),
    lb = c(NA, 0), ub = c(NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = "gumbel_lpdf(log(t_adj) | mu, sigma) - log(t_adj)",
    ldens = function(t, p) {
      lt <- log(t)
      z <- (lt - p$mu) / p$sigma
      -log(p$sigma) - z - exp(-z) - lt
    },
    stan_lcdf = "gumbel_lcdf(log(t_adj) | mu, sigma)",
    stan_lccdf = "gumbel_lccdf(log(t_adj) | mu, sigma)",
    # Gumbel: log F = -exp(-z), and the survival is its log1m_exp.
    lcdf = function(t, p) -exp(-(log(t) - p$mu) / p$sigma),
    lccdf = function(t, p) .log1m_exp(-exp(-(log(t) - p$mu) / p$sigma)),
    rng = function(n, p) exp(p$mu - p$sigma * log(-log(stats::runif(n)))),
    # E[exp(Gumbel)] = exp(mu) * Gamma(1 - sigma), and only exists for
    # sigma < 1. Note this is NOT exp(mu + sigma * gamma_euler), which is the
    # geometric mean rather than the mean.
    mean = function(p) ifelse(p$sigma < 1, exp(p$mu) * gamma(1 - p$sigma), Inf),
    init = list(mu = -0.8, sigma = 0.3),
    label = "Log-Weibull (Gumbel)"
  )
)


# The families built on the outlier mixture. with_outliers(),
# without_outliers(), p_outlier(), cogmod_priors() and cogmod_inits() all work
# on every one of them. The choice families in core_choice.R share the mixture
# but not the univariate density, so they have a registry of their own and are
# folded in here.
#' @keywords internal
.OUTLIER_FAMILIES <- c(names(.SHIFTED), .CHOICE_FAMILIES)


# The families that can take brms's `cens()` addition term: every registry
# entry carrying a survival, plus the two unshifted RT families, whose CDFs Stan
# ships. Derived rather than listed so that giving an entry its `stan_lccdf`
# is all it takes.
#
# The choice families are deliberately NOT here. Their likelihood is a set of
# defective densities summing to one over the response options, and the
# outlier's 1 / K factor exists to preserve that identity; a censored likelihood
# - a density for the winner, a bare survival for the loser - breaks it by
# construction, and brms would generate the `_lccdf` call regardless. Errors go
# through `dec()` there, and through `cens()` here.
#' @keywords internal
.CENS_FAMILIES <- c(
  names(Filter(function(e) !is.null(e$stan_lccdf), .SHIFTED)),
  "cogmod_exgaussian", "cogmod_geg"
)


#' @keywords internal
.shifted_spec <- function(name) {
  spec <- .SHIFTED[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the shifted RT families. Supported: ",
      paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Family constructor ------------------------------------------------------

# Builds the brms custom family for `name`, appending `ndt` and `poutlier` to
# whatever decision dpars the registry lists.
#' @keywords internal
.shifted_family <- function(name, links = NULL, predict_outliers = FALSE) {
  spec <- .shifted_spec(name)
  dec_links <- if (is.null(links)) spec$links else links
  if (length(dec_links) != length(spec$dpars)) {
    stop("Expected ", length(spec$dpars), " link functions for ", name, "().",
         call. = FALSE)
  }
  fam <- brms::custom_family(
    name = name,
    dpars = c(spec$dpars, "ndt", "poutlier"),
    links = c(dec_links, "log", "logit"),
    lb = c(spec$lb, 0, 0),
    ub = c(spec$ub, NA, 1),
    type = "real"
  )
  # It rides on the family because that is the only thing brms carries down to
  # a custom family's prediction methods. See ?rcogmod_lognormal.
  fam$predict_outliers <- isTRUE(predict_outliers)
  fam
}


# Stan code ---------------------------------------------------------------

# Generates the family's Stan code: the `<name>_lpdf` function - the same
# mixture skeleton for every family, with only the decision density swapped in -
# followed, for the families that have a CDF to offer, by the `<name>_lcdf` and
# `<name>_lccdf` that `cens()` needs (see .shifted_lcdfs()). The outlier component's
# scale is a literal because Stan functions cannot see the data block, and
# because a dpar would be estimated whenever the user left it out of the
# formula.
#' @keywords internal
.shifted_lpdf <- function(name, prelude = "") {
  spec <- .shifted_spec(name)
  # Families whose decision density needs helper functions declare them here.
  if (!is.null(spec$prelude)) prelude <- get(spec$prelude)
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(.POUTLIER_SCALE, format = "g", digits = 15, width = 1)

  # The outlier term log(2) + normal_lpdf(Y | 0, .POUTLIER_SCALE) has *no*
  # parameter in it - the location and the scale are both constants, and Y is
  # data. Called as written, Stan nevertheless recomputes its normalising
  # constant for every observation on every leapfrog step. Folding that constant
  # into a literal here and keeping only the Y-dependent part measured ~1.4x
  # faster per gradient evaluation on a 4000-observation LogNormal fit, with the
  # posterior unchanged.
  lc <- log(2) - 0.5 * log(2 * pi * .POUTLIER_SCALE^2)
  lp_out <- sprintf(
    "%s - %s * square(Y)",
    formatC(lc, format = "g", digits = 17, width = 1),
    formatC(1 / (2 * .POUTLIER_SCALE^2), format = "g", digits = 15, width = 1)
  )
  args <- paste(
    sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")),
    collapse = ", "
  )
  # Families that have something of their own to say about their parameters, or
  # about a boundary in their parameter space, carry it in the registry so the
  # generated code keeps it.
  dpar_doc <- if (is.null(spec$dpar_doc)) "" else {
    paste0(paste0("// ", spec$dpar_doc, collapse = "\n"), "\n")
  }
  note <- if (is.null(spec$note)) "" else {
    paste0("//\n", paste0("// ", spec$note, collapse = "\n"), "\n")
  }
  lpdf <- sprintf(
    "%s
// Log-likelihood for one observation from the shifted %s model.
// Y: observed reaction time.
%s// ndt: non-decision time, same unit as Y (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is a half Normal with scale %s s. It keeps the density
// strictly positive below `ndt`, where the shifted decision component has none.
// That is what removes the hard min-RT boundary and lets `ndt` be estimated
// directly rather than as a fraction of an observed minimum. The scale is a
// constant in SECONDS. This family expects reaction times in seconds; give it
// another unit and the component contributes nothing anywhere in the data,
// which silently reinstates the min-RT boundary it exists to remove.
//
// It is written out rather than called as normal_lpdf() because both of its
// parameters are constant: written that way, Stan recomputes the normalising
// constant for every observation on every leapfrog step.
%sreal %s_lpdf(%s) {
    // Parameter checks
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();

    // The leading constant includes the log(2) that folds the symmetric
    // Normal onto [0, Inf).
    real lp_out = %s;
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    return log_mix(poutlier, lp_out, %s);
}
",
    prelude, spec$label, dpar_doc, scale, note, name, args, spec$stan_check,
    lp_out, spec$stan_dens
  )
  # The CDF and the survival ride along whenever the family has them, so the
  # family's Stan code is complete for `cens()` without anyone asking for it.
  paste0(lpdf, .shifted_lcdfs(name))
}


# Generates `<name>_lcdf` and `<name>_lccdf`, the two functions brms's `cens()`
# addition term calls on a custom family:
#
#   if (cens[n] == 0)       target += <name>_lpdf(Y[n] | ...);
#   else if (cens[n] == 1)  target += <name>_lccdf(Y[n] | ...);   right-censored
#   else if (cens[n] == -1) target += <name>_lcdf(Y[n] | ...);    left-censored
#   else if (cens[n] == 2)  target += log_diff_exp(<name>_lcdf(rcens[n] | ...),
#                                                  <name>_lcdf(Y[n] | ...));
#
# Both are the mixture's own: p * G_out(Y) + (1 - p) * G_dec(Y - ndt), with G
# the CDF or the survival as appropriate. Below `ndt` the decision component
# has not started, so its CDF is 0 and its survival 1 - which is why a
# right-censored observation faster than `ndt` is not impossible: it says the
# decision process had not finished, and it had not.
#
# The half Normal's survival is 2 * Phi(-Y / scale), taken through
# std_normal_lcdf(-z). std_normal_lccdf(z) is the obvious call and the wrong
# one: it collapses to -inf once z passes about 8.3, which at a 0.2 s scale is
# Y = 1.66 s - the middle of the range a censored slow error lands in. Its CDF
# is 2 * Phi(Y / scale) - 1 = erf(Y / (scale * sqrt(2))), which has no such
# problem in either direction.
#' @keywords internal
.shifted_lcdfs <- function(name) {
  spec <- .shifted_spec(name)
  if (is.null(spec$stan_lcdf) || is.null(spec$stan_lccdf)) return("")
  num <- function(v) formatC(v, format = "g", digits = 17, width = 1)
  args <- paste(
    sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")),
    collapse = ", "
  )
  sprintf("
// Log CDF and log survival of the same mixture, for brms's cens() addition
// term. See ?rcogmod_invgaussian for what censoring a reaction time means and
// when it is the right model.
real %s_lcdf(%s) {
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();
    // Outlier CDF: 2 Phi(Y / s) - 1 = erf(Y / (s sqrt(2)))
    real lF_out = log(erf(Y * %s));
    real t_adj  = Y - ndt;
    if (t_adj <= 0) return log(poutlier) + lF_out;
    return log_mix(poutlier, lF_out, %s);
}

real %s_lccdf(%s) {
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return 0;
    // Outlier survival: 2 Phi(-Y / s), through the lower tail (see above)
    real lS_out = 0.69314718055994529 + std_normal_lcdf(-Y * %s);
    real t_adj  = Y - ndt;
    // Not yet past the non-decision time: the decision process cannot have
    // finished, so its survival is exactly 1.
    if (t_adj <= 0) return log_mix(poutlier, lS_out, 0);
    return log_mix(poutlier, lS_out, %s);
}
",
    name, args, spec$stan_check, num(1 / (.POUTLIER_SCALE * sqrt(2))),
    spec$stan_lcdf,
    name, args, spec$stan_check, num(1 / .POUTLIER_SCALE),
    spec$stan_lccdf
  )
}


# Density, RNG and mean ---------------------------------------------------

# Recycle the decision dpars plus ndt/poutlier to a common length. `...` holds
# the decision dpars, by name.
#' @keywords internal
.prepare_shifted <- function(name, x = NULL, n = NULL, ndt, poutlier, ...) {
  spec <- .shifted_spec(name)
  dec <- list(...)
  missing <- setdiff(spec$dpars, names(dec))
  if (length(missing)) {
    stop("Missing parameter(s) for ", name, "(): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  dec <- dec[spec$dpars]

  # The decision parameters are checked against the registry's own bounds, so
  # the R functions reject exactly what `stan_check` rejects and the two cannot
  # disagree. Lower bounds are open unless the family says otherwise - `sigma <=
  # 0` rather than `sigma < 0` - which is why the default comparison is `<=`;
  # cogmod_invgaussian() closes the bound on `sigmadrift`, a zero drift SD being
  # the plain Wald rather than an invalid model. Same convention, and same code,
  # as .prepare_choice() in core_choice.R.
  open <- if (is.null(spec$lb_open)) {
    rep(TRUE, length(spec$dpars))
  } else {
    spec$lb_open
  }
  for (j in seq_along(spec$dpars)) {
    d <- spec$dpars[j]
    if (!is.na(spec$lb[j])) {
      bad <- if (open[j]) dec[[d]] <= spec$lb[j] else dec[[d]] < spec$lb[j]
      if (any(bad, na.rm = TRUE)) {
        stop("`", d, "` must be ",
             if (open[j]) "greater than " else "at least ", spec$lb[j], ".",
             call. = FALSE)
      }
    }
    if (!is.na(spec$ub[j]) && any(dec[[d]] >= spec$ub[j], na.rm = TRUE)) {
      stop("`", d, "` must be less than ", spec$ub[j], ".", call. = FALSE)
    }
  }
  # ndt and poutlier are not in the registry: every family shares them, and
  # both admit their boundary values (ndt = 0 is an unshifted model, and
  # poutlier = 0 switches the outlier component off).
  if (any(ndt < 0, na.rm = TRUE)) stop("`ndt` must be non-negative.")
  if (any(poutlier < 0 | poutlier > 1, na.rm = TRUE)) {
    stop("`poutlier` must be in [0, 1].")
  }

  lens <- c(vapply(dec, length, integer(1)), length(ndt), length(poutlier))
  if (!is.null(x)) {
    # A zero-length quantile gives a zero-length answer, as it does for every
    # d/p/q function in base R. Taking `max()` over the parameter lengths
    # instead turned `numeric(0)` into a length-1 NA - `rep_len(numeric(0), 1)`
    # is NA - which then reached the density as a missing value.
    if (length(x) == 0L) {
      m <- 0L
    } else if (any(lens == 0L)) {
      stop("Parameters must not be zero-length when `x` is not.", call. = FALSE)
    } else {
      m <- max(length(x), lens)
    }
  } else if (!is.null(n)) {
    if (length(n) > 1) n <- length(n)
    if (length(n) != 1 || n < 0 || n != floor(n)) {
      stop("n must be a single non-negative integer.")
    }
    m <- max(n, lens)
    if (n == 0) m <- 0
  } else {
    stop("Either 'x' or 'n' must be provided.")
  }

  params <- lapply(dec, rep_len, m)
  params$ndt <- rep_len(ndt, m)
  params$poutlier <- rep_len(poutlier, m)
  if (!is.null(x)) params$x <- rep_len(x, m)
  params$ndraws <- m
  params
}


# Log-density of the decision component alone, at t (which may be <= 0, giving
# -Inf). Shape is preserved, so this works on draws x observations matrices.
#' @keywords internal
.ldec <- function(name, t, p) {
  spec <- .shifted_spec(name)
  # `.dens_mask()` substitutes a value the density will accept wherever the time
  # or a parameter is not one - non-positive, infinite or missing - and says
  # which entries those were; they are overwritten with -Inf immediately below.
  # A bare `pmax(t, 1e-300)` would cover the non-positive case but let NA
  # through to a branch that cannot take it. See `.dens_mask()`.
  msk <- .dens_mask(spec, t, p)
  ld <- spec$ldens(msk$t, msk$p)
  ld[!msk$ok] <- -Inf
  ld[is.na(ld)] <- -Inf
  dim(ld) <- dim(t)
  ld
}


# Mixture log-density, shared by every family's d*() function.
#' @keywords internal
.dshifted <- function(name, x, ndt, poutlier, log = FALSE, ...) {
  params <- tryCatch(
    .prepare_shifted(name, x = x, ndt = ndt, poutlier = poutlier, ...),
    error = function(e) {
      warning(conditionMessage(e), ". Returning 0 density / -Inf log-density.")
      list(ndraws = length(x), error = TRUE)
    }
  )
  if (!is.null(params$error)) {
    return(rep(ifelse(log, -Inf, 0), params$ndraws))
  }

  lp_out <- .dcontam(params$x, log = TRUE)
  lp_dec <- .ldec(name, params$x - params$ndt, params)

  ld <- .log_mix(params$poutlier, lp_out, lp_dec)
  ld[is.na(ld)] <- -Inf
  if (log) ld else exp(ld)
}


# Log CDF (or log survival) of the decision component alone, at t. The
# counterpart of .ldec(): the same mask, and the same reason for it. Outside the
# support the answer is known without evaluating anything - nothing has happened
# by a non-positive time, everything has by an infinite one - and a missing
# time stays missing, as it does in every p*() function in base R.
#' @keywords internal
.lcdf_dec <- function(name, t, p, lower.tail = TRUE) {
  spec <- .shifted_spec(name)
  if (is.null(spec$lcdf) || is.null(spec$lccdf)) {
    stop(name, "() has no closed-form CDF, so it cannot be censored. ",
         "Families that can: ", paste0(.CENS_FAMILIES, "()", collapse = ", "),
         ".", call. = FALSE)
  }
  msk <- .dens_mask(spec, t, p)
  out <- if (lower.tail) spec$lcdf(msk$t, msk$p) else spec$lccdf(msk$t, msk$p)
  # Wherever the time itself was the problem the answer is one of the two
  # boundary values; wherever a parameter was, it is -Inf like the density.
  below <- !is.na(t) & t <= 0
  above <- !is.na(t) & is.infinite(t) & t > 0
  out[!msk$ok] <- -Inf
  out[below] <- if (lower.tail) -Inf else 0
  out[above] <- if (lower.tail) 0 else -Inf
  out[is.na(t)] <- NA_real_
  dim(out) <- dim(t)
  out
}


# Mixture CDF, shared by every family's p*() function, and the R side of the
# Stan `<name>_lcdf` / `<name>_lccdf` pair (.shifted_lcdfs()). Everything is
# done on the log scale and in the tail that was asked for: the survival is
# never 1 - CDF, because a right-censored slow response sits exactly where that
# subtraction loses every digit.
#' @keywords internal
.pshifted <- function(name, q, ndt, poutlier, lower.tail = TRUE,
                      log.p = FALSE, ...) {
  params <- .prepare_shifted(name, x = q, ndt = ndt, poutlier = poutlier, ...)
  lp_out <- if (lower.tail) .lpcontam(params$x) else .lscontam(params$x)
  lp_dec <- .lcdf_dec(name, params$x - params$ndt, params, lower.tail)
  lp <- .log_mix(params$poutlier, lp_out, lp_dec)
  lp[is.na(params$x)] <- NA_real_
  if (log.p) lp else exp(lp)
}


# Mixture RNG, shared by every family's r*() function.
#' @keywords internal
.rshifted <- function(name, n, ndt, poutlier, ...) {
  spec <- .shifted_spec(name)
  params <- .prepare_shifted(name, n = n, ndt = ndt, poutlier = poutlier,
                                ...)
  m <- params$ndraws

  is_out <- stats::runif(m) < params$poutlier
  out <- numeric(m)

  if (any(is_out)) {
    out[is_out] <- abs(stats::rnorm(sum(is_out), 0, .POUTLIER_SCALE))
  }
  if (any(!is_out)) {
    keep <- !is_out
    p <- lapply(params[spec$dpars], function(v) v[keep])
    out[keep] <- spec$rng(sum(keep), p) + params$ndt[keep]
  }
  out
}


# brms methods ------------------------------------------------------------

# Pull the decision dpars for observation i (or all of them) off a brmsprep.
#' @keywords internal
.dpars_from_prep <- function(name, prep, i = NULL) {
  spec <- .shifted_spec(name)
  out <- lapply(spec$dpars, function(d) {
    if (is.null(i)) brms::get_dpar(prep, d) else brms::get_dpar(prep, d, i = i)
  })
  stats::setNames(out, spec$dpars)
}


#' @keywords internal
.log_lik_shifted <- function(name, i, prep) {
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  dec <- .dpars_from_prep(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)

  n_draws <- max(vapply(c(dec, list(ndt, poutlier)), length, integer(1)))
  if (n_draws == 0) return(numeric(0))

  # brms applies `cens()` to its own families inside log_lik() and leaves a
  # custom family's method to do it for itself - so without this branch loo()
  # would score every censored trial with the density of a response that was
  # never observed, and say nothing.
  ll <- .censor_ll(
    prep, i, y,
    ldens = function(y) do.call(.dshifted, c(
      list(name = name, x = rep(y, length.out = n_draws), ndt = ndt,
           poutlier = poutlier, log = TRUE),
      dec
    )),
    lcdf = function(y, lower.tail) do.call(.pshifted, c(
      list(name = name, q = rep(y, length.out = n_draws), ndt = ndt,
           poutlier = poutlier, lower.tail = lower.tail, log.p = TRUE),
      dec
    ))
  )
  ll[is.na(ll)] <- -Inf
  ll
}


# The log-likelihood of observation `i` under whatever `cens()` says about it,
# mirroring brms:::log_lik_censor() for the built-in families: 0 (or no `cens`
# at all) is an observed response, 1 right-censored, -1 left-censored and 2
# interval-censored between Y and `rcens`. `ldens(y)` and `lcdf(y, lower.tail)`
# return one log value per draw.
#' @keywords internal
.censor_ll <- function(prep, i, y, ldens, lcdf) {
  cens <- prep$data$cens[i]
  if (is.null(cens) || is.na(cens) || cens == 0) return(ldens(y))
  if (cens == 1) return(lcdf(y, lower.tail = FALSE))
  if (cens == -1) return(lcdf(y, lower.tail = TRUE))
  if (cens == 2) {
    return(.log_sub_exp(lcdf(prep$data$rcens[i], lower.tail = TRUE),
                        lcdf(y, lower.tail = TRUE)))
  }
  stop("Unknown censoring code ", cens, " in `cens()`.", call. = FALSE)
}


#' @keywords internal
.posterior_predict_shifted <- function(name, i, prep, predict_outliers = NULL) {
  dec <- .dpars_from_prep(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }
  n_draws <- max(vapply(c(dec, list(ndt)), length, integer(1)))

  out <- do.call(.rshifted, c(
    list(name = name, n = n_draws, ndt = ndt, poutlier = poutlier),
    dec
  ))
  as.matrix(out)
}


#' @keywords internal
.posterior_epred_shifted <- function(name, prep, predict_outliers = NULL) {
  spec <- .shifted_spec(name)
  dec <- .dpars_from_prep(name, prep)
  ndt <- brms::get_dpar(prep, "ndt")

  # get_dpar() returns a bare scalar for a dpar that is constant across draws
  # and observations, so the draws x observations shape is taken from whichever
  # one is actually a matrix.
  dims <- NULL
  for (d in c(dec, list(ndt))) {
    if (is.matrix(d)) {
      dims <- dim(d)
      break
    }
  }
  n <- if (is.null(dims)) max(vapply(dec, length, integer(1))) else prod(dims)
  dec <- lapply(dec, rep_len, n)

  if (is.null(spec$mean)) {
    stop(
      "The decision time of ", name, "() has no finite mean, so there is no ",
      "expectation to return. Use posterior_predict() and summarise the draws ",
      "instead.",
      call. = FALSE
    )
  }
  mean_dec <- spec$mean(dec)
  if (!is.null(dims) && length(mean_dec) == prod(dims)) dim(mean_dec) <- dims
  mean_dec <- mean_dec + ndt

  if (!.predict_outliers(predict_outliers, prep)) {
    return(mean_dec)
  }
  poutlier <- brms::get_dpar(prep, "poutlier")
  (1 - poutlier) * mean_dec + poutlier * .mcontam()
}


# log(Phi(x)) for Stan, shared by every prelude that takes a normal tail:
# .WALD_STAN_PRELUDE and .LOGNORMAL_STAN_PRELUDE below (and so the ex-Wald's
# and the LNR's), .RDM_STAN_PRELUDE in model_rdm.R and .EXGAUSSIAN_STAN_PRELUDE
# in model_exgaussian.R (and so the GEG's) paste it in front of their own
# functions. It is its own object so that helper-stan.R can strip the duplicate
# copies when it concatenates every family into one test program - a function
# defined twice will not compile - and so that the next family to need a normal
# tail takes this one rather than a Stan built-in. The reasoning is in the Stan
# comment.
#
# It sits here, ahead of every prelude, because those preludes paste it in at
# load time: an object defined further down the file does not exist yet.
#' @keywords internal
.LOG_PHI_STAN_PRELUDE <- "
// log(Phi(x)), the one piece of arithmetic every normal tail in this package is
// built from. Neither of Stan's two routes to it is good enough for both jobs
// it has here, which is to be right in the far tail *and* to hand back a
// usable derivative there. All three claims below were measured against
// central differences of the log probability over 20000 responses.
//
// The erfc route - std_normal_lcdf() is not it, but lognormal_lcdf(),
// lognormal_lccdf() and the log(u1) + log1m(u2 / u1) the LogNormal used to
// write are - has good partials, to about 4e-6 on a summed gradient of order
// 1e3. But erfc underflows near x = -38, and then the value is log(0) and the
// partials are inf or 0/0. That is not a harmless -inf in a mixture:
// log_mix() in the lpdf stays finite when the decision component is -inf, but
// reverse mode multiplies the (zero) adjoint into the stored partial, and
// 0 * inf is NaN, so a single response turns the gradient of the whole model
// to NaN - 'Gradient evaluated at the initial value is not finite' at the
// start of a fit, divergent transitions afterwards.
//
// std_normal_lcdf() has the range: its value is exact against R's
// pnorm(log.p = TRUE) as far as x = -1e7. Its partials are not - they sat
// 1.7e-3 from central differences where the erfc route sat 4e-6 on the LNR,
// and 2e-4 to 7e-2 on the RDM, which took every tail through it - so it is not
// a drop-in for the tails of a race, where those partials are the gradient of
// the drifts and the scales.
//
// So: erfc in the body of the distribution, and below x = -25 the asymptotic
// expansion of the tail,
//
//   Phi(x) = phi(x) / (-x) * (1 - 1/x^2 + 3/x^4 - 15/x^6 + 105/x^8 - 945/x^10)
//
// whose leading term is the exponent itself. Nothing underflows, the result
// stays finite and differentiable as far as x = -1e150, and the six terms
// agree with pnorm(log.p = TRUE) to 4e-16 relative from x = -25 down - the
// last bit of a double - so the two branches meet with no step in the density.
real cogmod_log_Phi(real x) {
  if (x < -25) {
    real z = inv_square(x);
    real series = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 - 945 * z))));
    return -0.5 * square(x) - log(-x) - 0.91893853320467274 + log(series);
  }
  if (x > 0) return log1p(-0.5 * erfc(x * 0.7071067811865476));
  return log(0.5 * erfc(-x * 0.7071067811865476));
}
"


# Unshifted Wald ----------------------------------------------------------

# The decision component of cogmod_invgaussian(), written out here rather than
# delegating to dcogmod_invgaussian()/rcogmod_invgaussian(): those are the *mixture*
# functions and route back through .dshifted(), so calling them from the
# registry would recurse.
#
# With `sigmadrift` = 0 this is the textbook Wald: an inverse Gaussian with mean
# boundary / drift and shape boundary^2. With sigmadrift > 0 the drift is drawn
# once per trial from Normal(drift, sigmadrift) truncated at zero, and the
# marginal density is
#
#   f(t) = boundary / sqrt(2 pi t^3 D) * exp(-(boundary - drift t)^2 / (2 t D))
#            * Phi(znum) / Phi(drift / sigmadrift),     D = 1 + sigmadrift^2 t
#
# with znum = (boundary sigmadrift^2 + drift) / (sigmadrift sqrt(D)). The first
# two factors are the integral of the Wald over an *untruncated* Normal drift,
# which is a Gaussian integral in the drift and so stays elementary; the ratio
# of normal CDFs is what the truncation at zero contributes - the numerator is
# the mass above zero after seeing t, the denominator the mass above zero
# before.
#
# D = 1 at sigmadrift = 0, and the two CDFs then both have an infinite argument
# and cancel, so the expression degenerates to the plain Wald exactly rather
# than approaching it. It is still written as a branch, because at sigmadrift =
# 0 the CDF arguments are 0/0 rather than infinite.
#
# The tail is worth knowing about: as t grows the exponent tends to a constant
# and the prefactor to t^-2, so the density decays as t^-2 whenever sigmadrift
# is positive, and E[T] does not exist. That is the registry's `mean` returning
# Inf.
#
# `sigmandt` smears the shift over Uniform(0, sigmandt), and the density is then
# a difference of CDFs rather than anything in this closed form - see
# .lwald_st0_fixed(). It is routed through .lwald_raw() for exactly the elements
# that need it, so a zero range costs nothing and returns bit for bit what it
# did before the parameter existed.
#' @keywords internal
.dwald_raw <- function(t, drift, boundary, sigmadrift = 0, sigmandt = 0) {
  s2 <- sigmadrift^2
  D <- 1 + s2 * t
  ld <- log(boundary) - 0.5 * (log(2 * pi) + 3 * log(t) + log(D)) -
    (boundary - drift * t)^2 / (2 * t * D)

  # The division below is by `s`, which stands in as 1 wherever sigmadrift is 0
  # so that the truncation term stays finite; multiplying by the indicator then
  # drops it. Written this way rather than by subsetting because `t` and the
  # parameters may be draws x observations matrices, whose shape has to survive.
  s <- sigmadrift + (sigmadrift <= 0)
  ld <- ld + (sigmadrift > 0) *
    (.lpnorm_upper((boundary * s2 + drift) / (s * sqrt(D))) -
      .lpnorm_upper(drift / s))

  if (!any(sigmandt > 0, na.rm = TRUE)) return(ld)
  # Logical indexing keeps whatever shape `ld` has.
  n <- length(ld)
  st0 <- rep_len(sigmandt, n)
  use <- !is.na(st0) & st0 > 0
  ld[use] <- .lwald_raw(rep_len(t, n)[use], rep_len(drift, n)[use],
                        rep_len(boundary, n)[use], rep_len(sigmadrift, n)[use],
                        st0[use], "dens")
  ld
}

# log(Phi(z)), taken through the lower tail as log1p(-Phi(-z)), because Phi(z)
# is the quantity that saturates at 1 for the large positive z a well-separated
# drift produces. Same form, and same reason, as .dlba1_raw()'s normalizer.
#' @keywords internal
.lpnorm_upper <- function(z) {
  log1p(-exp(stats::pnorm(-z, log.p = TRUE)))
}

# Michael-Schucany-Haas two-root method, applied to the drift the trial actually
# got. .rnorm_truncated() returns the mean untouched when sd is 0, but the draw
# is only taken when some trial really has a variable drift, which keeps the
# random stream of the plain Wald exactly as it was. The same goes for the
# non-decision offset: Uniform(0, sigmandt) is drawn only when some trial has a
# positive range. `ndt` itself is added by .rshifted(), so what is returned is
# the decision time plus the trial's excess over the lower bound.
#' @keywords internal
.rwald_raw <- function(n, drift, boundary, sigmadrift = 0, sigmandt = 0) {
  if (any(sigmadrift > 0, na.rm = TRUE)) {
    drift <- .rnorm_truncated(n, mean = drift, sd = sigmadrift, lower = 0)
  }
  ig_mu <- boundary / drift
  lambda <- boundary^2
  y <- stats::rnorm(n)^2
  z <- y * (ig_mu / lambda)
  x1_over_mu <- 1 + z / 2 * (1 - sqrt(1 + 4 / z))
  u <- stats::runif(n)
  out <- ig_mu * ifelse(u < 1 / (1 + x1_over_mu), x1_over_mu, 1 / x1_over_mu)
  if (any(sigmandt > 0, na.rm = TRUE)) {
    # runif() with max == min returns min, so a zero range adds exactly 0.
    out <- out + stats::runif(n, 0, rep_len(sigmandt, n))
  }
  out
}

# CDF of the fixed-drift decision component. `exp(2 * boundary * drift) *
# Phi(.)` is taken as `exp(2 * boundary * drift + lcdf)` for the same reason as
# in .wald_lcdf_direct(): the factor overflows on its own well before the
# product stops being representable.
#' @keywords internal
.pwald_fixed <- function(t, drift, boundary) {
  st <- sqrt(t)
  stats::pnorm((drift * t - boundary) / st) +
    exp(2 * boundary * drift +
      stats::pnorm(-(drift * t + boundary) / st, log.p = TRUE))
}

# The same CDF marginalised over a drift drawn from Normal(drift, sigmadrift)
# truncated at zero. Unlike the density, this integral is *not* elementary: it
# is a bivariate normal probability, so there is nothing to write in closed form
# short of Phi_2. It is taken by 64-point Gauss-Legendre quadrature over the
# drift instead, on the interval carrying the truncated normal's mass. The
# integrand is analytic there, so the rule converges geometrically: measured
# against adaptive integration of the density, the error is at machine precision
# (~5e-15) across drift 0.5-5, boundary 0.2-2 and sigmadrift 0.05-3, and 48
# nodes already suffice everywhere but the widest of those.
#' @keywords internal
.pwald_sv <- function(t, drift, boundary, sigmadrift, nsd = 10) {
  lo <- pmax(drift - nsd * sigmadrift, 0)
  hi <- drift + nsd * sigmadrift
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  # The truncation's normalizer, so that the weights below sum to one.
  lnorm <- .lpnorm_upper(drift / sigmadrift)

  out <- numeric(length(t))
  for (j in seq_along(.GAUSS_LEGENDRE$x)) {
    v <- mid + half * .GAUSS_LEGENDRE$x[j]
    w <- .GAUSS_LEGENDRE$w[j] * half *
      exp(stats::dnorm(v, drift, sigmadrift, log = TRUE) - lnorm)
    out <- out + w * .pwald_fixed(t, v, boundary)
  }
  out
}

# Dispatch on whether the drift varies at all, so a fixed-drift Wald keeps its
# exact closed-form CDF. `t` is expected to be finite and strictly positive; the
# parameters are recycled against it.
#' @keywords internal
.pwald_raw <- function(t, drift, boundary, sigmadrift = 0) {
  n <- length(t)
  drift <- rep_len(drift, n)
  boundary <- rep_len(boundary, n)
  sigmadrift <- rep_len(sigmadrift, n)

  sv <- sigmadrift > 0
  out <- numeric(n)
  if (any(!sv)) {
    out[!sv] <- .pwald_fixed(t[!sv], drift[!sv], boundary[!sv])
  }
  if (any(sv)) {
    out[sv] <- .pwald_sv(t[sv], drift[sv], boundary[sv], sigmadrift[sv])
  }
  pmin(pmax(out, 0), 1)
}

# Gauss-Legendre nodes and weights on [-1, 1], from the Golub-Welsch
# eigendecomposition of the Jacobi matrix. Built once when the namespace loads;
# .pwald_sv() is the only thing that needs them.
#' @keywords internal
.GAUSS_LEGENDRE <- local({
  k <- 64L
  i <- seq_len(k - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, k, k)
  jacobi[cbind(i, i + 1L)] <- b
  jacobi[cbind(i + 1L, i)] <- b
  e <- eigen(jacobi, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
})

# Stan counterpart of .dwald_raw(). The branch is the same one, and so is the
# reason for it: at sigmadrift = 0 the two normal truncation factors are 0/0
# rather than 1.
#
# Every normal tail below is cogmod_log_Phi(), never std_normal_lcdf(): the
# built-in's value is right but its partials are not, and `sigmadrift` is
# differentiated almost entirely through the two truncation factors of the
# closed-form branch. Measured 2026-09-18, this is what the gradient check was
# seeing - 1.2e-3 relative on d/d sigmadrift, 2.7e-4 on d/d mu, against 1e-8
# for the families that were already on cogmod_log_Phi(). The truncation factor
# used to be written log1m_exp(std_normal_lcdf(-z)), the lower tail because
# with a positive drift and threshold the argument is positive and Phi itself
# saturates there; cogmod_log_Phi(z) is that same quantity computed from
# log1p(-erfc(z / sqrt(2)) / 2), so the saturation is gone along with the
# round trip through log1m_exp().
#' @keywords internal
#
# The CDF and survival follow, for cens(). At sigmadrift = 0 both are closed
# form; above it they are the fixed-drift ones marginalised over the truncated
# normal drift by the same 64-point Gauss-Legendre rule as .lwald_sv(), with
# the nodes and log-weights written out from .GAUSS_LEGENDRE so that R and Stan
# integrate over literally the same points. A non-decision range (sigmandt)
# enters through cogmod_wald_st0, the counterpart of .lwald_st0_fixed(): the
# same closed forms, the same tail choice, the same small-range switch, so the
# two sides take identical branches. The helpers that are not Stan
# `_lcdf`/`_lccdf` functions are named to avoid those suffixes: Stan would
# otherwise insist on the `|` call syntax for them.
#' @keywords internal
.WALD_STAN_PRELUDE <- local({
  num <- function(v) formatC(v, format = "g", digits = 17, width = 1)
  k <- length(.GAUSS_LEGENDRE$x)
  paste0(.LOG_PHI_STAN_PRELUDE, "
// Log density of the fixed-drift Wald with drift v > 0 and threshold a > 0,
// at t > 0: an inverse Gaussian with mean a / v and shape a^2.
real cogmod_wald_ldens(real t, real v, real a) {
  return log(a) - 0.5 * (log(2 * pi()) + 3 * log(t)) - square(a - v * t) / (2 * t);
}

// Log CDF of the same Wald. The exp(2 a v) factor overflows on its own long
// before the product it belongs to stops being representable, so it is folded
// into the exponent.
real cogmod_wald_logcdf(real t, real v, real a) {
  if (t <= 0) return negative_infinity();
  real st = sqrt(t);
  return log_sum_exp(
    cogmod_log_Phi((v * t - a) / st),
    2 * a * v + cogmod_log_Phi(-(v * t + a) / st)
  );
}

// Log survival of the same Wald: the DIFFERENCE of the two terms,
//   S(t) = Phi((a - v t) / sqrt(t)) - exp(2 a v) Phi(-(a + v t) / sqrt(t)),
// rather than log1m_exp(logcdf), which has lost every digit by the time a slow
// error is censored. Same form as cogmod_rdm_wald_lsurv() at A = 0.
real cogmod_wald_lsurv(real t, real v, real a) {
  if (t <= 0) return 0;
  real st = sqrt(t);
  real m1 = cogmod_log_Phi((a - v * t) / st);
  real m2 = 2 * a * v + cogmod_log_Phi(-(a + v * t) / st);
  return m1 > m2 ? log_diff_exp(m1, m2) : negative_infinity();
}

// log INT_0^t F(s) ds, the integrated CDF, at t > 0:
//   I(t) = t F(t) - (a / v) [Phi(alpha) - exp(2 a v) Phi(beta)],
// alpha = (v t - a) / sqrt(t), beta = -(v t + a) / sqrt(t). The bracket is
// (v / a) E[T; T <= t] and t F(t) exceeds its multiple, so both differences
// are taken in log space and are positive.
real cogmod_wald_liF(real t, real v, real a) {
  real st = sqrt(t);
  real lPa = cogmod_log_Phi((v * t - a) / st);
  real lPb = 2 * a * v + cogmod_log_Phi(-(v * t + a) / st);
  real lF = log_sum_exp(lPa, lPb);
  real lP = lPa > lPb ? log_diff_exp(lPa, lPb) : negative_infinity();
  real x = log(t) + lF;
  real y = log(a / v) + lP;
  return x > y ? log_diff_exp(x, y) : negative_infinity();
}

// log INT_t^Inf S(s) ds, the integrated survival:
//   R(t) = (a / v - t) Phi(-alpha) + (a / v + t) exp(2 a v) Phi(beta),
// which is a / v - t below the shift (S = 1 there). Above the mean the first
// coefficient is negative and R is a difference, positive because it is the
// integral of a survival.
real cogmod_wald_liS(real t, real v, real a) {
  if (t <= 0) return log(a / v - t);
  real st = sqrt(t);
  real lQa = cogmod_log_Phi(-(v * t - a) / st);
  real lQb = 2 * a * v + cogmod_log_Phi(-(v * t + a) / st);
  real c1 = a / v - t;
  real c2 = a / v + t;
  if (c1 >= 0) return log_sum_exp(log(c1) + lQa, log(c2) + lQb);
  real x = log(c2) + lQb;
  real y = log(-c1) + lQa;
  return x > y ? log_diff_exp(x, y) : negative_infinity();
}

// The fixed-drift Wald with its shift smeared over Uniform(0, st0): the log
// density (what = 0), log CDF (1) or log survival (2) of the decision
// component seen at t = Y - ndt. Each is a difference quotient of the
// functions above,
//   f = [F(t) - F(t - st0)] / st0 = [S(t - st0) - S(t)] / st0
//   G = [I(t) - I(t - st0)] / st0,   1 - G = [R(t - st0) - R(t)] / st0,
// taken between whichever pair is small - the CDFs near the shift, the
// survivals in the tail, decided by the midpoint's CDF - and reduced to the
// one-sided form when the interval reaches back past the shift. Below
// st0 / t = 1e-5 the quotient has lost five digits and the midpoint value,
// with error O(st0^2), is the better answer. Same branches, same order, as
// .lwald_st0_fixed() in R.
real cogmod_wald_st0(real t, real v, real a, real st0, int what) {
  if (st0 <= 0 || st0 < 1e-5 * t) {
    real m = st0 <= 0 ? t : t - 0.5 * st0;
    if (what == 0) return cogmod_wald_ldens(m, v, a);
    if (what == 1) return cogmod_wald_logcdf(m, v, a);
    return cogmod_wald_lsurv(m, v, a);
  }
  real t0 = t - st0;
  real ls = log(st0);
  if (what == 0) {
    if (t0 <= 0) return cogmod_wald_logcdf(t, v, a) - ls;
    if (cogmod_wald_logcdf(t - 0.5 * st0, v, a) < -0.69314718055994529) {
      real x = cogmod_wald_logcdf(t, v, a);
      real y = cogmod_wald_logcdf(t0, v, a);
      return x > y ? log_diff_exp(x, y) - ls : negative_infinity();
    }
    real x = cogmod_wald_lsurv(t0, v, a);
    real y = cogmod_wald_lsurv(t, v, a);
    return x > y ? log_diff_exp(x, y) - ls : negative_infinity();
  }
  if (what == 1) {
    real x = cogmod_wald_liF(t, v, a);
    if (t0 <= 0) return x - ls;
    real y = cogmod_wald_liF(t0, v, a);
    return x > y ? log_diff_exp(x, y) - ls : negative_infinity();
  }
  real x = cogmod_wald_liS(t0, v, a);
  real y = cogmod_wald_liS(t, v, a);
  return x > y ? log_diff_exp(x, y) - ls : negative_infinity();
}

// Any of the three above, marginalised over a drift ~ Normal(mu, sigmadrift)
// truncated at zero, by ", k, "-point Gauss-Legendre quadrature on the interval
// carrying the truncated normal's mass. Assembled with log_sum_exp so the
// survival keeps its digits where every term is tiny.
real cogmod_wald_sv_lquad(real t, real mu, real boundary, real sigmadrift,
                          real st0, int what) {
  vector[", k, "] gx = [", paste(num(.GAUSS_LEGENDRE$x), collapse = ", "), "]';
  vector[", k, "] lgw = [", paste(num(log(.GAUSS_LEGENDRE$w)), collapse = ", "), "]';
  real lo = fmax(mu - 10 * sigmadrift, 0);
  real hi = mu + 10 * sigmadrift;
  real half = 0.5 * (hi - lo);
  real mid = 0.5 * (hi + lo);
  real lnorm = cogmod_log_Phi(mu / sigmadrift);
  vector[", k, "] terms;
  for (j in 1:", k, ") {
    real v = mid + half * gx[j];
    real lw = lgw[j] + log(half) + normal_lpdf(v | mu, sigmadrift) - lnorm;
    terms[j] = lw + cogmod_wald_st0(t, v, boundary, st0, what);
  }
  return log_sum_exp(terms);
}

// Log density of the Wald decision time (no shift), with the drift rate drawn
// once per trial from Normal(mu, sigmadrift) truncated at zero and the
// non-decision time spread over a range sigmandt above ndt. sigmadrift = 0 and
// sigmandt = 0 is the plain Wald: an inverse Gaussian with mean boundary / mu
// and shape boundary^2. With sigmandt = 0 the drift is marginalised in closed
// form - a Gaussian integral - which is the branch a fit usually runs in.
real cogmod_invgaussian_decision_lpdf(real t, real mu, real boundary,
                                      real sigmadrift, real sigmandt) {
  if (sigmandt <= 0) {
    real base = log(boundary) - 0.5 * (log(2 * pi()) + 3 * log(t));
    if (sigmadrift <= 0) {
      return base - square(boundary - mu * t) / (2 * t);
    }
    real s2 = square(sigmadrift);
    real D = 1 + s2 * t;
    return base - 0.5 * log(D) - square(boundary - mu * t) / (2 * t * D)
      + cogmod_log_Phi((boundary * s2 + mu) / (sigmadrift * sqrt(D)))
      - cogmod_log_Phi(mu / sigmadrift);
  }
  if (sigmadrift <= 0) return cogmod_wald_st0(t, mu, boundary, sigmandt, 0);
  return cogmod_wald_sv_lquad(t, mu, boundary, sigmadrift, sigmandt, 0);
}

// Log CDF and log survival of the Wald decision time (no shift), the
// counterparts of cogmod_invgaussian_decision_lpdf().
real cogmod_invgaussian_decision_lcdf(real t, real mu, real boundary,
                                      real sigmadrift, real sigmandt) {
  if (sigmadrift <= 0) return cogmod_wald_st0(t, mu, boundary, sigmandt, 1);
  return cogmod_wald_sv_lquad(t, mu, boundary, sigmadrift, sigmandt, 1);
}

real cogmod_invgaussian_decision_lccdf(real t, real mu, real boundary,
                                       real sigmadrift, real sigmandt) {
  if (sigmadrift <= 0) return cogmod_wald_st0(t, mu, boundary, sigmandt, 2);
  return cogmod_wald_sv_lquad(t, mu, boundary, sigmadrift, sigmandt, 2);
}
")
})


# ex-Wald -----------------------------------------------------------------

# The decision time of cogmod_exwald() is Wald(drift, boundary) + Exp(1 / tau).
# Writing g = 1 / tau, the convolution
#
#   f(t) = g e^{-g t} INT_0^t f_Wald(u) e^{g u} du
#
# has an elementary closed form, because f_Wald(u) e^{g u} is itself a Wald
# density with a different drift:
#
#   f_Wald(u; v, a) e^{g u} = e^{a (v - k)} f_Wald(u; k, a),   k = sqrt(v^2 - 2g)
#
# so the integral is e^{a(v - k)} F_Wald(t; k, a) and
#
#   f(t) = g exp(-g t + a (v - k)) F_Wald(t; k, a).
#
# Exact - it agrees with brute-force convolution to 5e-14 - but only while
# v^2 > 2g, i.e. tau > 2 / drift^2. That is NOT the usual regime: at a drift of
# 3 and a threshold of 0.5 it demands tau > 0.22 s, where a residual stage is
# more often nearer 0.1 s, so the other branch is reached routinely.
#
# Below the line k is imaginary, k = i * kappa with kappa = sqrt(2g - v^2), and
# the same expression continues analytically. Because F_Wald's two terms are
# complex conjugates of one another there,
#
#   INT_0^t f_Wald(u) e^{g u} du = 2 e^{a v} Re[e^{-i a kappa} Phi(-alpha + i beta)]
#
# with alpha = a / sqrt(t) and beta = kappa sqrt(t). Writing Phi through the
# Faddeeva function w(z) = e^{-z^2} erfc(-i z) - which is what makes this
# computable, w being analytic and well conditioned throughout the upper half
# plane - the substitution alpha * beta = a * kappa makes every large factor
# cancel on paper, leaving
#
#   f(t) = g exp(-(a - v t)^2 / (2 t)) Re[w(z)],
#     z = (kappa sqrt(t) + i a / sqrt(t)) / sqrt(2).
#
# The exponent is the Wald's own, so nothing here can overflow, and the whole
# branch costs one real exponential and one w(). At kappa -> 0 it reduces to
# `g exp(a v - v^2 t / 2) erfc(a / sqrt(2t))`, which is the closed form at
# k = 0, so the two branches meet exactly rather than approximately: measured
# either side of the seam the relative step is 5e-8, against 0.74 for the
# 64-point Gauss-Legendre rule this replaced.
#
# Quadrature was tried first and abandoned. The log-integrand of the original
# convolution is bimodal - d/du = 0 is k2 u^2 + 3 u - a^2 = 0, which has two
# interior roots whenever 9 > 4 kappa^2 a^2, one at the Wald bulk near zero and
# one at u = t where exp(g u) is climbing - and the two peak widths vary
# independently over orders of magnitude, so no fixed panel layout resolves
# both. Single-panel, exponentially tilted, split-at-the-mode and
# z = a / sqrt(u) rules all reached only 1e-1 relative error somewhere in the
# region an RT fit actually visits.
#' @keywords internal
.dexwald_raw <- function(t, drift, boundary, tau) {
  g <- 1 / tau
  k2 <- drift^2 - 2 * g

  closed <- function() {
    k <- sqrt(pmax(k2, 0))
    log(g) - g * t + boundary * (drift - k) + .log_pwald_fixed(t, k, boundary)
  }
  faddeeva <- function() {
    st <- sqrt(t)
    rw <- .re_faddeeva(sqrt(pmax(-k2, 0)) * st / sqrt(2),
                       boundary / (st * sqrt(2)))
    # Re w is the Voigt function, strictly positive off the real axis; the guard
    # is for a rounding-level negative where it has underflowed.
    log(g) - (boundary - drift * t)^2 / (2 * t) + log(pmax(rw, 0))
  }

  if (all(k2 >= 0)) return(closed())
  if (all(k2 < 0)) return(faddeeva())

  # Both branches occur in the same call, which happens as soon as `tau` or the
  # drift is modelled. Computed whole and selected rather than subset, because
  # `t` may be a draws x observations matrix whose shape has to survive - the
  # same reason .dwald_raw() masks instead of indexing.
  use <- (k2 >= 0) & rep(TRUE, length(t))
  dim(use) <- dim(t)
  ifelse(use, closed(), faddeeva())
}

# log F_Wald(t; drift, boundary) for a fixed drift. The same expression as
# .pwald_fixed(), carried in logs throughout: the CDF underflows to exactly zero
# well before the log density stops being representable, and `exp(2 a v)`
# overflows on its own long before the product does.
#' @keywords internal
.log_pwald_fixed <- function(t, drift, boundary) {
  st <- sqrt(t)
  l1 <- stats::pnorm((drift * t - boundary) / st, log.p = TRUE)
  l2 <- 2 * boundary * drift +
    stats::pnorm(-(drift * t + boundary) / st, log.p = TRUE)
  m <- pmax(l1, l2)
  out <- m + log(exp(l1 - m) + exp(l2 - m))
  # Both terms vanish; pmax() is -Inf and the difference would be NaN.
  out[!is.finite(m)] <- -Inf
  out
}

# Log survival of the fixed-drift Wald: the DIFFERENCE of the same two terms,
#
#   S(t) = Phi((a - v t) / sqrt(t)) - exp(2 a v) Phi(-(a + v t) / sqrt(t)),
#
# which is what stays accurate as t grows. log(1 - F) does not: F is within
# rounding of 1 by the time a slow error is censored, and the subtraction
# returns 0 or -Inf. Same form as the RDM's cogmod_rdm_wald_lsurv() at a zero
# start-point range, and the same two bmm bugs it exists to avoid - -Inf at
# t = ndt and NaN in the upper tail.
#' @keywords internal
.lswald_fixed <- function(t, drift, boundary) {
  st <- sqrt(t)
  l1 <- stats::pnorm((boundary - drift * t) / st, log.p = TRUE)
  l2 <- 2 * boundary * drift +
    stats::pnorm(-(boundary + drift * t) / st, log.p = TRUE)
  .log_sub_exp(l1, l2)
}

# Log density of the fixed-drift Wald: the sigmadrift = 0 branch of
# .dwald_raw(), on its own so the st0 kernels below can call it at a shifted
# time without going back through the mixture.
#' @keywords internal
.ldwald_fixed <- function(t, drift, boundary) {
  log(boundary) - 0.5 * (log(2 * pi) + 3 * log(t)) -
    (boundary - drift * t)^2 / (2 * t)
}

# The Wald with its shift smeared over Uniform(0, st0): the non-decision time
# of a trial is ndt + U, so the decision component seen at t_adj = Y - ndt is
# the average of the fixed-shift one over U in [0, st0]. Averaging a density
# integrates it, averaging a CDF integrates that, and for the Wald every one of
# those integrals is closed form (F and S are the fixed-drift CDF and survival,
# alpha = (v t - a) / sqrt(t), beta = -(v t + a) / sqrt(t)):
#
#   f_st0(t) = [F(t) - F(t - st0)] / st0 = [S(t - st0) - S(t)] / st0
#   G_st0(t) = [I(t) - I(t - st0)] / st0,
#       I(t) = INT_0^t F = t F(t) - (a / v) [Phi(alpha) - e^{2 a v} Phi(beta)]
#   1 - G_st0(t) = [R(t - st0) - R(t)] / st0,
#       R(t) = INT_t^Inf S = (a / v - t) Phi(-alpha) + (a / v + t) e^{2 a v} Phi(beta)
#
# with F, I = 0 and S = 1, R(s) = a / v - s below the shift, which is what the
# interval straddling it (t < st0) reduces to. Both integrals were checked
# against adaptive integration of F and S before being written down: 4e-16
# absolute for I, 5e-14 relative for R.
#
# Three numerical points, none of them new to this file:
#
#  * `e^{2 a v} Phi(beta)` is kept as one exponent, as in .log_pwald_fixed():
#    the factor overflows on its own long before the product does.
#  * A difference is taken between whichever pair is SMALL. The density is
#    F(t) - F(t - st0) near the shift and S(t - st0) - S(t) in the tail, and
#    the two differ only in which side has lost its digits; the midpoint's CDF
#    decides. Same device as `dPhi` in .lba_dens_over_A(). R inherits the
#    digit loss of .lswald_fixed() in the far tail and no more.
#  * Below st0 / t = 1e-5 the difference quotient has lost five digits and the
#    midpoint value f(t - st0 / 2) - error O(st0^2) - is the better answer, so
#    the kernel switches to it there, as the RDM does below .RDM_EPS_A. The log
#    link never reaches zero, but the sampler visits 1e-4 s routinely, and the
#    two branches have to meet. Measured at drift 3, boundary 0.5 over t from
#    0.21 to 6 s, the step at the switch is below 4e-9 relative, and 1e-5 is
#    where it is smallest: at 1e-4 the midpoint's truncation error has grown
#    to 3e-7, at 1e-6 the quotient's cancellation has.
#
# `what` is "dens", "cdf" or "surv", all on the log scale. Every argument is a
# vector of the same length, and st0 may be zero for some elements: those get
# the plain kernels, so this is the one function the drift quadrature has to
# call. The Stan side (cogmod_wald_st0) is the same branches in the same order.
#' @keywords internal
.lwald_st0_fixed <- function(t, v, a, st0, what) {
  plain <- switch(what, dens = .ldwald_fixed, cdf = .log_pwald_fixed,
                  surv = .lswald_fixed)
  n <- length(t)
  out <- numeric(n)

  # No range, or a range too small to difference across: the point value, at
  # the midpoint of the interval where there is one.
  point <- st0 <= 0 | st0 < 1e-5 * t
  if (any(point)) {
    m <- t[point] - pmax(st0[point], 0) / 2
    out[point] <- plain(m, v[point], a[point])
  }

  g <- !point
  if (!any(g)) return(out)
  t1 <- t[g]
  st <- st0[g]
  t0 <- t1 - st
  vv <- v[g]
  aa <- a[g]
  ls <- log(st)
  left <- t0 <= 0 # the interval reaches back past the shift
  res <- numeric(length(t1))

  if (what == "dens") {
    lF1 <- .log_pwald_fixed(t1, vv, aa)
    res[left] <- lF1[left] - ls[left]
    if (any(!left)) {
      i <- !left
      lower <- .log_pwald_fixed(t1[i] - st[i] / 2, vv[i], aa[i]) < log(0.5)
      d <- numeric(sum(i))
      if (any(lower)) {
        j <- which(i)[lower]
        d[lower] <- .log_sub_exp(lF1[j], .log_pwald_fixed(t0[j], vv[j], aa[j]))
      }
      if (any(!lower)) {
        j <- which(i)[!lower]
        d[!lower] <- .log_sub_exp(.lswald_fixed(t0[j], vv[j], aa[j]),
                                  .lswald_fixed(t1[j], vv[j], aa[j]))
      }
      res[i] <- d - ls[i]
    }
  } else if (what == "cdf") {
    lI1 <- .lwald_liF_fixed(t1, vv, aa)
    res[left] <- lI1[left] - ls[left]
    if (any(!left)) {
      i <- !left
      res[i] <- .log_sub_exp(lI1[i], .lwald_liF_fixed(t0[i], vv[i], aa[i])) -
        ls[i]
    }
  } else {
    lR1 <- .lwald_liS_fixed(t1, vv, aa)
    lR0 <- numeric(length(t1))
    lR0[left] <- log(aa[left] / vv[left] - t0[left])
    if (any(!left)) {
      lR0[!left] <- .lwald_liS_fixed(t0[!left], vv[!left], aa[!left])
    }
    res <- .log_sub_exp(lR0, lR1) - ls
  }
  out[g] <- res
  out
}

# log I(t) = log INT_0^t F(s) ds for the fixed-drift Wald, t > 0. The bracket
# Phi(alpha) - e^{2 a v} Phi(beta) is (v / a) E[T; T <= t], positive for every
# t, and t F(t) exceeds (a / v) times it because F is increasing; both
# differences are taken in log space. Near the shift I is about 2 t / a^2 of
# t F(t), so a digit or two goes there; nothing is lost where anything is
# observable.
#' @keywords internal
.lwald_liF_fixed <- function(t, v, a) {
  st <- sqrt(t)
  lPa <- stats::pnorm((v * t - a) / st, log.p = TRUE)
  lPb <- 2 * a * v + stats::pnorm(-(v * t + a) / st, log.p = TRUE)
  lF <- .log_add_exp(lPa, lPb)
  lP <- .log_sub_exp(lPa, lPb)
  .log_sub_exp(log(t) + lF, log(a / v) + lP)
}

# log R(t) = log INT_t^Inf S(s) ds for the fixed-drift Wald, t > 0. Below the
# mean a / v both coefficients are positive and the terms add; above it the
# first is negative and R is the difference, positive because it is the
# integral of a survival, and losing digits in the far tail at the rate S
# itself does. R(0) is the mean.
#' @keywords internal
.lwald_liS_fixed <- function(t, v, a) {
  st <- sqrt(t)
  lQa <- stats::pnorm(-(v * t - a) / st, log.p = TRUE)
  lQb <- 2 * a * v + stats::pnorm(-(v * t + a) / st, log.p = TRUE)
  c1 <- a / v - t
  c2 <- a / v + t
  out <- numeric(length(t))
  pos <- c1 >= 0
  if (any(pos)) {
    out[pos] <- .log_add_exp(log(c1[pos]) + lQa[pos], log(c2[pos]) + lQb[pos])
  }
  if (any(!pos)) {
    i <- !pos
    out[i] <- .log_sub_exp(log(c2[i]) + lQb[i], log(-c1[i]) + lQa[i])
  }
  out
}

# The fixed-drift log density, log CDF or log survival, marginalised over a
# drift drawn from Normal(drift, sigmadrift) truncated at zero - the same
# 64-point Gauss-Legendre rule as .pwald_sv(), on the same interval, but
# assembled with log_sum_exp so that the survival keeps its digits in the tail.
# The kernel at each node is .lwald_st0_fixed(), so a non-decision range rides
# along for free. The Stan side (cogmod_wald_sv_lquad) is generated from the
# same node table, so the two cannot disagree.
#' @keywords internal
.lwald_sv <- function(t, drift, boundary, sigmadrift, st0, what, nsd = 10) {
  lo <- pmax(drift - nsd * sigmadrift, 0)
  hi <- drift + nsd * sigmadrift
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  lnorm <- .lpnorm_upper(drift / sigmadrift)

  k <- length(.GAUSS_LEGENDRE$x)
  terms <- matrix(NA_real_, length(t), k)
  for (j in seq_len(k)) {
    v <- mid + half * .GAUSS_LEGENDRE$x[j]
    lw <- log(.GAUSS_LEGENDRE$w[j]) + log(half) +
      stats::dnorm(v, drift, sigmadrift, log = TRUE) - lnorm
    terms[, j] <- lw + .lwald_st0_fixed(t, v, boundary, st0, what)
  }
  m <- apply(terms, 1, max)
  out <- m + log(rowSums(exp(terms - m)))
  out[!is.finite(m)] <- -Inf
  out
}

# Log density, log CDF and log survival of the decision component of
# cogmod_invgaussian(), dispatching on whether the drift varies. `t` finite and
# strictly positive; the parameters are recycled against it. A fixed drift goes
# straight to the closed forms; a variable one to the quadrature - except the
# density at a zero range, which .dwald_raw() has in closed form for a variable
# drift too, and which is the common case.
#' @keywords internal
.lwald_raw <- function(t, drift, boundary, sigmadrift, sigmandt, what) {
  n <- length(t)
  drift <- rep_len(drift, n)
  boundary <- rep_len(boundary, n)
  sigmadrift <- rep_len(sigmadrift, n)
  sigmandt <- rep_len(sigmandt, n)
  sv <- sigmadrift > 0
  out <- numeric(n)
  if (any(!sv)) {
    out[!sv] <- .lwald_st0_fixed(t[!sv], drift[!sv], boundary[!sv],
                                 sigmandt[!sv], what)
  }
  if (any(sv)) {
    q <- sv
    if (what == "dens") {
      z <- sv & sigmandt <= 0
      if (any(z)) out[z] <- .dwald_raw(t[z], drift[z], boundary[z], sigmadrift[z])
      q <- sv & sigmandt > 0
    }
    if (any(q)) {
      out[q] <- .lwald_sv(t[q], drift[q], boundary[q], sigmadrift[q],
                          sigmandt[q], what)
    }
  }
  if (what != "dens") out <- pmin(out, 0)
  out
}

#' @keywords internal
.lpwald_raw <- function(t, drift, boundary, sigmadrift = 0, sigmandt = 0) {
  .lwald_raw(t, drift, boundary, sigmadrift, sigmandt, "cdf")
}

#' @keywords internal
.lswald_raw <- function(t, drift, boundary, sigmadrift = 0, sigmandt = 0) {
  .lwald_raw(t, drift, boundary, sigmadrift, sigmandt, "surv")
}

# Coefficients of Weideman's (1994) rational approximation to the Faddeeva
# function, built once when the namespace loads - the same arrangement as
# .GAUSS_LEGENDRE, and for the same reason: the Stan prelude is generated from
# this object, so the two implementations cannot drift apart.
#
# N = 24 terms. The approximation is exact to 8e-11 at w(0), and N = 32 gives
# the ex-Wald density to the last digit N = 24 already reaches, so the extra
# terms buy nothing here.
#' @keywords internal
.WEIDEMAN <- local({
  N <- 24L
  M <- 2L * N
  L <- sqrt(N / sqrt(2))
  k <- seq(-M + 1L, M - 1L)
  tt <- L * tan(k * pi / M / 2)
  f <- c(0, exp(-tt^2) * (L^2 + tt^2))
  n <- length(f)
  # fftshift
  fs <- c(f[(floor(n / 2) + 1):n], f[1:floor(n / 2)])
  a <- Re(stats::fft(fs)) / (2 * M)
  list(N = N, L = L, a = rev(a[2:(N + 1)]))
})

# Re w(x + i y) for y > 0, by Horner on Weideman's rational approximation.
# Written out in real and imaginary parts rather than with R's complex type,
# because the Stan counterpart has to do exactly the same arithmetic in the same
# order, and because it keeps the draws x observations shape of the inputs.
#' @keywords internal
.re_faddeeva <- function(x, y) {
  L <- .WEIDEMAN$L
  a <- .WEIDEMAN$a
  # d = L - i z with z = x + i y, so d = (L + y) - i x
  dr <- L + y
  di <- -x
  den <- dr^2 + di^2
  # Z = (L + i z) / d, numerator (L - y) + i x
  nr <- L - y
  ni <- x
  zr <- (nr * dr + ni * di) / den
  zi <- (ni * dr - nr * di) / den

  pr <- a[1] + 0 * zr
  pim <- 0 * zr
  for (n in seq_along(a)[-1]) {
    tr <- pr * zr - pim * zi + a[n]
    pim <- pr * zi + pim * zr
    pr <- tr
  }

  # w = 2 p / d^2 + (1 / sqrt(pi)) / d, real part only
  d2r <- dr^2 - di^2
  d2i <- 2 * dr * di
  q <- d2r^2 + d2i^2
  2 * (pr * d2r + pim * d2i) / q + (1 / sqrt(pi)) * dr / den
}

# Stan counterpart. The Weideman coefficients are written out from .WEIDEMAN, so
# the two implementations are generated from one table. The Wald prelude comes
# first because the closed-form branch is written in the Wald's log CDF
# (cogmod_wald_logcdf) - one definition, shared, rather than a copy here.
#' @keywords internal
.EXWALD_STAN_PRELUDE <- local({
  num <- function(v) formatC(v, format = "g", digits = 17, width = 1)
  paste0(.WALD_STAN_PRELUDE, sprintf("
// Re w(x + i y) for y > 0, where w is the Faddeeva function
// w(z) = exp(-z^2) erfc(-i z). Weideman (1994), %d-term rational approximation,
// evaluated by Horner in explicit real and imaginary parts.
real cogmod_re_faddeeva(real x, real y) {
  real L = %s;
  array[%d] real a = {%s};
  real dr = L + y;
  real di = -x;
  real den = square(dr) + square(di);
  real zr = ((L - y) * dr + x * di) / den;
  real zi = (x * dr - (L - y) * di) / den;
  real pr = a[1];
  real pim = 0;
  for (n in 2:%d) {
    real tr = pr * zr - pim * zi + a[n];
    pim = pr * zi + pim * zr;
    pr = tr;
  }
  real d2r = square(dr) - square(di);
  real d2i = 2 * dr * di;
  real q = square(d2r) + square(d2i);
  return 2 * (pr * d2r + pim * d2i) / q + 0.56418958354775628 * dr / den;
}

// Log density of the ex-Wald decision time (no shift): a Wald with drift `mu`
// and threshold `boundary`, convolved with an Exponential of mean `tau`.
//
// Where mu^2 > 2 / tau the convolution collapses to a closed form in the Wald
// CDF, with drift k = sqrt(mu^2 - 2 / tau). Below that k is imaginary and the
// same expression continues analytically into
//
//   f(t) = g exp(-(boundary - mu t)^2 / (2 t)) Re[w(z)],
//     z = (kappa sqrt(t) + i boundary / sqrt(t)) / sqrt(2)
//
// with kappa = sqrt(2 / tau - mu^2). The exponent is the Wald's own, so nothing
// overflows. The two branches meet exactly - at kappa = 0 the second reduces to
// the first - and measured either side of the seam the relative step is 5e-8.
real cogmod_exwald_decision_lpdf(real t, real mu, real boundary, real tau) {
  real g = inv(tau);
  real k2 = square(mu) - 2 * g;

  if (k2 >= 0) {
    real k = sqrt(k2);
    return log(g) - g * t + boundary * (mu - k)
      + cogmod_wald_logcdf(t, k, boundary);
  }

  real st = sqrt(t);
  real rw = cogmod_re_faddeeva(sqrt(-k2) * st / sqrt2(), boundary / (st * sqrt2()));
  if (rw <= 0) return negative_infinity();
  return log(g) - square(boundary - mu * t) / (2 * t) + log(rw);
}
",
    .WEIDEMAN$N, num(.WEIDEMAN$L), .WEIDEMAN$N,
    paste(num(.WEIDEMAN$a), collapse = ", "), .WEIDEMAN$N))
})


# Unshifted Birnbaum-Saunders -----------------------------------------------

# The decision component of cogmod_bisa(). Written in the same (drift,
# threshold) parameters as .dwald_raw(), and in those parameters the density is
# the Wald's own with one extra factor:
#
#   f(t) = (mu t + boundary) / (2 sqrt(2 pi) t^{3/2})
#            * exp(-(mu t - boundary)^2 / (2 t))
#        = f_Wald(t; mu, boundary) * (mu t + boundary) / (2 boundary).
#
# That is not a coincidence and it is worth seeing why, because it is the whole
# character of the family. In the textbook (a, b) parameters the exponent is
# -z^2 / 2 with z = (1 / a)(sqrt(t / b) - sqrt(b / t)); substituting
# b = boundary / mu and a = 1 / sqrt(mu * boundary) gives a * sqrt(b) = 1 / mu,
# so z collapses to (mu t - boundary) / sqrt(t) - the Wald's own exponent
# numerator - and the prefactor (sqrt(t / b) + sqrt(b / t)) / (2 a t) collapses
# to (mu t + boundary) / (2 t^{3/2}). One sign is all that separates them.
#
# The tilt factor is (1 + t / b) / 2, which is what makes the family an equal
# mixture of that Wald and its length-biased version: the sum of the two
# first-passage routes to the same boundary, one of which weights long crossings
# in proportion to their length. It costs one log and one square more than
# nothing, so this is the cheapest density in the package after the LogNormal.
#
# `t` is expected finite and strictly positive; the parameters are recycled
# against it, and may arrive as draws x observations matrices.
#' @keywords internal
.dbisa_raw <- function(t, drift, boundary) {
  log(drift * t + boundary) - 0.5 * (log(2 * pi) + 3 * log(t)) - log(2) -
    (drift * t - boundary)^2 / (2 * t)
}

# Exact, and a one-liner, because (1 / a)(sqrt(T / b) - sqrt(b / T)) is exactly
# standard normal: invert that for T and a single Normal draw gives a single
# Birnbaum-Saunders draw. No rejection, no two-root correction of the kind
# .rwald_raw() needs.
#' @keywords internal
.rbisa_raw <- function(n, drift, boundary) {
  b <- boundary / drift
  h <- stats::rnorm(n) / (2 * sqrt(drift * boundary))
  b * (h + sqrt(1 + h^2))^2
}

# Stan counterpart of .dbisa_raw(). No branch and no special function: unlike
# every other first-passage family here, this one is closed form everywhere in
# the parameter space.
#' @keywords internal
.BISA_STAN_PRELUDE <- "
// Log density of the Birnbaum-Saunders decision time (no shift), with evidence
// arriving in discrete cycles of average size `mu` and unit SD, towards a
// threshold `boundary`. (mu * t - boundary) / sqrt(t) is then exactly standard
// normal, which leaves the density elementary:
//
//   f(t) = f_Wald(t; mu, boundary) * (mu t + boundary) / (2 boundary),
//
// the Wald's own density tilted by the length-biasing factor. Note the sign:
// the exponent carries (mu t - boundary), the prefactor (mu t + boundary).
real cogmod_bisa_decision_lpdf(real t, real mu, real boundary) {
  return log(mu * t + boundary) - 0.5 * (log(2 * pi()) + 3 * log(t)) - log(2)
    - square(mu * t - boundary) / (2 * t);
}
"


# Log CDF and log survival of the half Normal outlier component, each computed
# in its own tail. The survival goes through pnorm(-x / scale) rather than
# log1p(-pnorm(x / scale)): the latter is -Inf once x / scale passes about 8,
# which is 1.66 s here - see .shifted_lcdfs() for the Stan side of the same
# point.
#' @keywords internal
.lpcontam <- function(x) {
  out <- log(2 * stats::pnorm(x / .POUTLIER_SCALE) - 1)
  out[!is.na(x) & x <= 0] <- -Inf
  out
}

#' @keywords internal
.lscontam <- function(x) {
  out <- log(2) + stats::pnorm(-x / .POUTLIER_SCALE, log.p = TRUE)
  out[!is.na(x) & x <= 0] <- 0
  out
}

# CDF of the half Normal outlier component.
#' @keywords internal
.pcontam <- function(x) {
  out <- 2 * stats::pnorm(x / .POUTLIER_SCALE) - 1
  out[x <= 0] <- 0
  out
}


# Shared expose helper: compiles the generated lpdf and hands it back as an R
# function, for checking the Stan code against the R density.
#' @keywords internal
.shifted_expose <- function(name) {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {
",
    .shifted_lpdf(name),
    "}"
  )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions[[paste0(name, "_lpdf")]]
}


# Unshifted LBA, one accumulator at a time --------------------------------

# Shared by cogmod_lba1() (RT only, one accumulator) and cogmod_lba2() (choice,
# two racing accumulators), which differ in how many accumulators they run and
# in what they condition on, not in the per-accumulator arithmetic. They live
# here rather than with either family for the same reason as the mixture itself:
# one copy cannot drift out of step with the other.

# The LBA defective density divided by the start-point range A, i.e.
#
#   [ drift * (Phi(z2) - Phi(z1)) + sigma * (phi(z1) - phi(z2)) ] / A
#
# with z2 = z1 + delta and delta = A / (sigma * t).
#
# Both differences vanish linearly in delta, so evaluating them directly and
# then dividing by A loses every significant digit once the start-point range is
# small: at A = 1e-2 the density is already wrong, and by A = 1e-4 it no longer
# integrates to one. (The previous implementation also floored the bracket at
# 1e-10, which turned that underflow into a spurious density floor spread over
# the whole line, so the "density" integrated to far more than one.)
#
# Below delta = 1e-4 the Taylor expansion in delta is used instead. Its
# truncation error there is ~1e-12 relative, while the direct form still has
# ~12 good digits, so the two agree across the switch. Above it, the two
# differences are taken tail-aware: Phi(z2) - Phi(z1) from the upper tail when
# both are positive, and phi(z1) - phi(z2) as phi(z1) * -expm1(...), which stays
# accurate however close the two arguments are.
#' @keywords internal
.lba_dens_over_A <- function(drift, sigma, st, z1, delta) {
  n <- max(length(drift), length(sigma), length(st), length(z1), length(delta))
  drift <- rep_len(drift, n)
  sigma <- rep_len(sigma, n)
  st <- rep_len(st, n)
  z1 <- rep_len(z1, n)
  delta <- rep_len(delta, n)

  phi1 <- stats::dnorm(z1)
  out <- numeric(n)

  small <- delta < 1e-4
  if (any(small)) {
    i <- small
    d <- delta[i]
    z <- z1[i]
    v <- drift[i]
    sg <- sigma[i]
    series <- (v + sg * z) -
      (d / 2) * (v * z + sg * (z^2 - 1)) +
      (d^2 / 6) * (v * (z^2 - 1) + sg * (z^3 - 3 * z))
    out[i] <- phi1[i] * series / st[i]
  }
  if (any(!small)) {
    i <- !small
    z2 <- z1[i] + delta[i]
    dPhi <- ifelse(
      z1[i] > 0,
      stats::pnorm(z1[i], lower.tail = FALSE) -
        stats::pnorm(z2, lower.tail = FALSE),
      stats::pnorm(z2) - stats::pnorm(z1[i])
    )
    dphi <- phi1[i] * -expm1(-delta[i] * (z1[i] + z2) / 2)
    out[i] <- (drift[i] * dPhi + sigma[i] * dphi) / (delta[i] * st[i])
  }
  out
}


# Log-survival of one accumulator at `t` - the probability it has NOT yet
# reached the threshold - given that its drift is positive. Only cogmod_lba2()
# needs it: a race scores the winner's density against the loser's survival.
# It is the same integral over the same start-point range as the density, so it
# belongs beside it.
#
# The drift is a Normal truncated at zero (see ?rcogmod_lba2), so the survival
# wanted is P(v > 0, unfinished at t) / P(v > 0). Writing b - A - v t = z1 * st
# and b - v t = z2 * st, the UNtruncated textbook CDF
#
#   F(t) = 1 + ((b - A - vt)/A) Phi(z1) - ((b - vt)/A) Phi(z2)
#            + (t s / A) (phi(z1) - phi(z2))
#
# becomes `(h(z2) - h(z1)) / delta` with `h(z) = z Phibar(z) - phi(z)`, and the
# untruncated survival `S = 1 - F` becomes `(g(z2) - g(z1)) / delta` with
# `g(z) = z Phi(z) + phi(z)`. A finished accumulator always had a positive
# drift, so P(v > 0, unfinished) = S - q = (1 - q) - F, with q = Phi(-drift /
# sigma) the probability of a negative one. Which form to take is decided by
# the sign of the drift, because that is what decides which quantities are
# small: with a positive drift q is at most a half and S - q cancels only in
# the far tail, where the answer decays like 1 / t and there is nothing to
# lose; with a negative drift both 1 - q and F are tiny, and taking them from
# the upper tail - Phibar rather than 1 - Phi, h rather than g - keeps every
# digit. Both quotients cancel in their own right as delta -> 0, hence the same
# Taylor branch as the density: g has derivatives Phi, phi, -z phi and h has
# Phibar, -phi, z phi.
#
# Returned on the log scale, as log(S - q) - log(1 - q), so that a loser that
# has almost surely finished contributes a large negative number rather than
# log(0).
#' @keywords internal
.lba_lsurv_trunc <- function(drift, sigma, z1, delta) {
  n <- max(length(drift), length(sigma), length(z1), length(delta))
  drift <- rep_len(drift, n)
  sigma <- rep_len(sigma, n)
  z1 <- rep_len(z1, n)
  delta <- rep_len(delta, n)

  phi1 <- stats::dnorm(z1)
  ratio <- drift / sigma
  small <- delta < 1e-4
  neg <- drift < 0
  num <- numeric(n)

  # positive drift: S - q, with S the untruncated survival
  i <- !neg & small
  if (any(i)) {
    d <- delta[i]
    z <- z1[i]
    s <- stats::pnorm(z) + (d / 2) * phi1[i] - (d^2 / 6) * z * phi1[i]
    num[i] <- s - stats::pnorm(-ratio[i])
  }
  i <- !neg & !small
  if (any(i)) {
    z2 <- z1[i] + delta[i]
    g <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
    num[i] <- (g(z2) - g(z1[i])) / delta[i] - stats::pnorm(-ratio[i])
  }
  # negative drift: (1 - q) - F, with F the untruncated CDF, all upper tail
  i <- neg & small
  if (any(i)) {
    d <- delta[i]
    z <- z1[i]
    f <- stats::pnorm(z, lower.tail = FALSE) - (d / 2) * phi1[i] +
      (d^2 / 6) * z * phi1[i]
    num[i] <- stats::pnorm(ratio[i]) - f
  }
  i <- neg & !small
  if (any(i)) {
    z2 <- z1[i] + delta[i]
    h <- function(z) z * stats::pnorm(z, lower.tail = FALSE) - stats::dnorm(z)
    num[i] <- stats::pnorm(ratio[i]) - (h(z2) - h(z1[i])) / delta[i]
  }

  out <- rep(-Inf, n)
  ok <- is.finite(num) & num > 0
  if (any(ok)) {
    out[ok] <- pmin(log(num[ok]) - stats::pnorm(ratio[ok], log.p = TRUE), 0)
  }
  out
}


# Decision component of cogmod_lba1(), written out here for the same reason as the
# Wald: dcogmod_lba1() is the mixture and would recurse.
#' @keywords internal
.dlba1_raw <- function(t, drift, sigma, sigmabias, boundary) {
  st <- pmax(sigma * t, 1e-10)
  z1 <- (boundary - drift * t) / st
  n2 <- .lba_dens_over_A(drift, sigma, st, z1, sigmabias / st)
  log_norm <- log1p(-exp(stats::pnorm(-drift / sigma, log.p = TRUE)))
  ifelse(n2 > 0, log(n2) - log_norm, -Inf)
}

#' @keywords internal
.rlba1_raw <- function(n, drift, sigma, sigmabias, boundary) {
  v <- .rnorm_truncated(n, mean = drift, sd = sigma, lower = 0)
  sp <- stats::runif(n, min = 0, max = sigmabias)
  (sigmabias + boundary - sp) / v
}

# Stan counterpart of .lba_dens_over_A(), .lba_lsurv_trunc() and .dlba1_raw().
# Same branches, and the same reason for them: both differences cancel as the
# start-point range goes to zero, and the result is then divided by that range.
#
# cogmod_lba2() appends its own decision lpdf to this (see .LBA2_STAN_PRELUDE in
# model_lba2.R) rather than repeating the two kernels, so the names here are
# family-neutral.
#' @keywords internal
.LBA_STAN_PRELUDE <- "
// LBA defective density divided by the start-point range A, with
// z2 = z1 + delta and delta = A / (sigma * t). Both differences below vanish
// linearly in delta, so evaluating them directly and dividing by A loses every
// significant digit for a small start-point range. The Taylor expansion is used
// there instead; its truncation error is ~1e-12 at the switch, where the direct
// form still has ~12 good digits.
real cogmod_lba_dens_over_A(real drift, real sigma, real st, real z1, real delta) {
  real phi1 = exp(-0.5 * square(z1)) * 0.3989422804014327;
  if (delta < 1e-4) {
    real series = (drift + sigma * z1)
      - (delta / 2) * (drift * z1 + sigma * (square(z1) - 1))
      + (square(delta) / 6) * (drift * (square(z1) - 1)
                               + sigma * (z1 * square(z1) - 3 * z1));
    return phi1 * series / st;
  }
  real z2 = z1 + delta;
  // upper tail when both are positive, so Phi(z2) - Phi(z1) does not cancel
  real dPhi = z1 > 0 ? (Phi(-z1) - Phi(-z2)) : (Phi(z2) - Phi(z1));
  real dphi = phi1 * -expm1(-delta * (z1 + z2) / 2);
  return (drift * dPhi + sigma * dphi) / (delta * st);
}

// Log probability that an accumulator has NOT finished by t, given that its
// drift is positive: log(S - q) - log(1 - q), with S the untruncated survival
// (g(z2) - g(z1)) / delta, g(z) = z Phi(z) + phi(z), and q = Phi(-drift /
// sigma) the probability of a negative drift. With a negative drift both
// 1 - q and the untruncated CDF F are tiny, so the same number is taken as
// (1 - q) - F from the upper tail, F = (h(z2) - h(z1)) / delta with
// h(z) = z Phibar(z) - phi(z). The Taylor branches are for the cancellation in
// the quotients as delta -> 0; see .lba_lsurv_trunc() for the derivation.
real cogmod_lba_lsurv_trunc(real drift, real sigma, real z1, real delta) {
  real phi1 = exp(-0.5 * square(z1)) * 0.3989422804014327;
  real ratio = drift / sigma;
  real num;
  if (drift >= 0) {
    real s;
    if (delta < 1e-4) {
      s = Phi(z1) + (delta / 2) * phi1 - (square(delta) / 6) * z1 * phi1;
    } else {
      real z2 = z1 + delta;
      real phi2 = exp(-0.5 * square(z2)) * 0.3989422804014327;
      s = ((z2 * Phi(z2) + phi2) - (z1 * Phi(z1) + phi1)) / delta;
    }
    num = s - Phi(-ratio);
  } else {
    real f;
    if (delta < 1e-4) {
      f = Phi(-z1) - (delta / 2) * phi1 + (square(delta) / 6) * z1 * phi1;
    } else {
      real z2 = z1 + delta;
      real phi2 = exp(-0.5 * square(z2)) * 0.3989422804014327;
      f = ((z2 * Phi(-z2) - phi2) - (z1 * Phi(-z1) - phi1)) / delta;
    }
    num = Phi(ratio) - f;
  }
  if (num <= 0) return negative_infinity();
  return fmin(log(num) - std_normal_lcdf(ratio), 0);
}

// Log density of the single-accumulator LBA decision time (no shift).
real cogmod_lba1_decision_lpdf(real t, real drift, real sigma, real sigmabias, real boundary) {
  real st = fmax(sigma * t, 1e-10);
  real z1 = (boundary - drift * t) / st;
  real f = cogmod_lba_dens_over_A(drift, sigma, st, z1, sigmabias / st);
  if (f <= 0) return negative_infinity();
  return log(f) - log1m_exp(std_normal_lcdf(-drift / sigma));
}
"

# LogNormal accumulator with a start-point range --------------------------

# The single accumulator that cogmod_lognormal() is and that cogmod_lnr()
# races two of. It runs from a start point z ~ Uniform(0, A) to the threshold
# 1 + A at a rate v ~ LogNormal(-meanlog, sigma), so its finishing time is
# (1 + A - z) / v, with the distance D = 1 + A - z ~ Uniform(1, 1 + A). The
# threshold *offset* above the highest start point is pinned at 1 - the LBA's
# `boundary` convention with boundary = 1 - because for a LogNormal rate the
# scale cannot be pinned by sigma: rescaling the evidence axis shifts the rate's
# location and scales A and the threshold, and leaves sigma alone. At A = 0 the
# distance is 1, T = 1 / v, and log T ~ Normal(meanlog, sigma): the shifted
# LogNormal as it always was, with `mu` and `sigma` meaning what they meant.
# See ?rcogmod_lognormal and ?rcogmod_lnr.
#
# Write a = (meanlog - log t) / sigma, so that F = Phi(-a) and S = Phi(a) at
# A = 0, and c = log1p(A) / sigma, the width of the start-point range on the
# z-scale. Averaging over the Uniform distance turns the density into a partial
# first moment of the rate distribution, which for a LogNormal is a difference
# of two normal CDFs, and integrating the CDF by parts gives the survival:
#
#   f(t) = exp(-meanlog + sigma^2 / 2) * [Phi(a + c - sigma) - Phi(a - sigma)] / A
#   S(t) = Phi(a + c) + [Phi(a + c) - Phi(a)] / A  -  t * f(t)
#
# Both differences vanish linearly in c, so dividing by A loses every digit as
# A -> 0; below c = 1e-4 the series in c is used instead, at the same switch as
# .lba_dens_over_A(). At A = 0 exactly the plain lognormal functions are
# called, so a model that pins `sigmabias = 0` computes what the family
# computed before the parameter existed, at the same cost.
#
# The survival's two middle terms are each of order |a| / sigma in the late tail
# and cancel to something of order 1, so both tails are assembled in log space
# from ratios of CDFs, on whichever side of a = 0 keeps the bracket well
# behaved: S directly for a < 0, F directly for a >= 0 (where F <= 1/2), and the
# other tail as log(1 - exp(.)) of it. `cens()` needs the CDF written that way
# round rather than as 1 - S; see the registry notes on `stan_lcdf`.

# log(Phi(y + c) - Phi(y)) for c > 0, from whichever tail keeps the two terms
# from cancelling: the upper tail when y > 0, where both CDFs sit near 1, the
# lower tail otherwise. Exact in the far tails, where the plain difference is
# 0 - 0. Vectorised over both arguments.
#' @keywords internal
.lognormal_ldiff_pnorm <- function(y, c) {
  n <- max(length(y), length(c))
  y <- rep_len(y, n)
  c <- rep_len(c, n)
  out <- rep(-Inf, n)
  up <- y > 0
  if (any(up)) {
    l1 <- stats::pnorm(y[up], lower.tail = FALSE, log.p = TRUE)
    l2 <- stats::pnorm(y[up] + c[up], lower.tail = FALSE, log.p = TRUE)
    ok <- l2 < l1
    out[up][ok] <- l1[ok] + .log1m_exp(l2[ok] - l1[ok])
  }
  if (any(!up)) {
    l1 <- stats::pnorm(y[!up] + c[!up], log.p = TRUE)
    l2 <- stats::pnorm(y[!up], log.p = TRUE)
    ok <- l2 < l1
    out[!up][ok] <- l1[ok] + .log1m_exp(l2[ok] - l1[ok])
  }
  out
}


# Log density of the accumulator's finishing time at t > 0. `A` may be 0.
#' @keywords internal
.lognormal_acc_ldens <- function(t, meanlog, sigma, A) {
  n <- max(length(t), length(meanlog), length(sigma), length(A))
  t <- rep_len(t, n)
  meanlog <- rep_len(meanlog, n)
  sigma <- rep_len(sigma, n)
  A <- rep_len(A, n)
  out <- rep(-Inf, n)

  zero <- A == 0
  if (any(zero)) {
    out[zero] <- stats::dlnorm(t[zero], meanlog[zero], sigma[zero], log = TRUE)
  }
  if (all(zero)) return(out)

  i <- which(!zero)
  ti <- t[i]
  mi <- meanlog[i]
  si <- sigma[i]
  Ai <- A[i]
  a <- (mi - log(ti)) / si
  c <- log1p(Ai) / si
  x <- a - si
  val <- rep(-Inf, length(i))

  small <- c < 1e-4
  if (any(small)) {
    # [Phi(x + c) - Phi(x)] / A = (c / A) phi(x) [1 - c x / 2 + c^2 (x^2 - 1) / 6]
    # and exp(-meanlog + sigma^2 / 2) phi(x) = phi(a) / t, so the density is
    # the LogNormal's times log1p(A) / A times the series.
    j <- small
    series <- 1 - c[j] * x[j] / 2 + c[j]^2 * (x[j]^2 - 1) / 6
    ok <- series > 0
    v <- rep(-Inf, sum(j))
    v[ok] <- stats::dlnorm(ti[j][ok], mi[j][ok], si[j][ok], log = TRUE) +
      log(log1p(Ai[j][ok]) / Ai[j][ok]) + log(series[ok])
    val[j] <- v
  }
  if (any(!small)) {
    j <- !small
    val[j] <- -mi[j] + si[j]^2 / 2 + .lognormal_ldiff_pnorm(x[j], c[j]) -
      log(Ai[j])
  }
  out[i] <- val
  out
}


# Both tails of the accumulator's finishing time at t > 0, as a list of the
# log-CDF and the log-survival, each computed directly on the side where it is
# the small one. `A` may be 0.
#' @keywords internal
.lognormal_acc_ltails <- function(t, meanlog, sigma, A) {
  n <- max(length(t), length(meanlog), length(sigma), length(A))
  t <- rep_len(t, n)
  meanlog <- rep_len(meanlog, n)
  sigma <- rep_len(sigma, n)
  A <- rep_len(A, n)
  lF <- rep(-Inf, n)
  lS <- rep(-Inf, n)

  zero <- A == 0
  if (any(zero)) {
    lF[zero] <- stats::plnorm(t[zero], meanlog[zero], sigma[zero], log.p = TRUE)
    lS[zero] <- stats::plnorm(t[zero], meanlog[zero], sigma[zero],
                              lower.tail = FALSE, log.p = TRUE)
  }
  if (all(zero)) return(list(lcdf = lF, lccdf = lS))

  i <- which(!zero)
  ti <- t[i]
  mi <- meanlog[i]
  si <- sigma[i]
  Ai <- A[i]
  a <- (mi - log(ti)) / si
  c <- log1p(Ai) / si
  vF <- rep(-Inf, length(i))
  vS <- rep(-Inf, length(i))

  small <- c < 1e-4
  if (any(small)) {
    # Expanding (sigma / A) * integral_0^c Phi(a + u) exp(sigma u) du in u:
    #   S = Phi(a) + c r phi(a) [1/2 + c (2 sigma - a) / 6],  r = log1p(A) / A,
    # and F = Phi(-a) minus the same correction. phi / Phi is the inverse
    # Mills ratio, bounded by about |a| + 1 on either side.
    j <- small
    r <- log1p(Ai[j]) / Ai[j]
    corr <- c[j] * r * (0.5 + c[j] * (2 * si[j] - a[j]) / 6)
    lPa <- stats::pnorm(a[j], log.p = TRUE)
    lQa <- stats::pnorm(-a[j], log.p = TRUE)
    lphi <- stats::dnorm(a[j], log = TRUE)
    vS[j] <- lPa + log1p(corr * exp(lphi - lPa))
    dF <- 1 - corr * exp(lphi - lQa)
    w <- rep(-Inf, sum(j))
    w[dF > 0] <- lQa[dF > 0] + log(dF[dF > 0])
    vF[j] <- w
  }
  if (any(!small)) {
    j <- which(!small)
    aj <- a[j]
    cj <- c[j]
    sj <- si[j]
    Aj <- Ai[j]
    lA <- log(Aj)
    # D1 = [Phi(a + c) - Phi(a)] / A and D2 = t f(t), both in log space
    lD1 <- .lognormal_ldiff_pnorm(aj, cj) - lA
    lD2 <- -mi[j] + sj^2 / 2 + log(ti[j]) +
      .lognormal_ldiff_pnorm(aj - sj, cj) - lA
    wF <- rep(-Inf, length(j))
    wS <- rep(-Inf, length(j))
    late <- aj < 0
    if (any(late)) {
      # S = Phi(a + c) [1 + D1 / Phi(a + c) - D2 / Phi(a + c)]
      lP <- stats::pnorm(aj[late] + cj[late], log.p = TRUE)
      br <- 1 + exp(lD1[late] - lP) - exp(lD2[late] - lP)
      ok <- br > 0
      w <- rep(-Inf, sum(late))
      w[ok] <- pmin(lP[ok] + log(br[ok]), 0)
      wS[late] <- w
      wF[late] <- ifelse(w < 0, .log1m_exp(w), -Inf)
    }
    if (any(!late)) {
      # F = Phi(-a) [R + (R - 1) / A + D2 / Phi(-a)], R = Phi(-a-c) / Phi(-a)
      e <- !late
      lQ <- stats::pnorm(-aj[e], log.p = TRUE)
      R <- exp(stats::pnorm(-aj[e] - cj[e], log.p = TRUE) - lQ)
      br <- R - (1 - R) / Aj[e] + exp(lD2[e] - lQ)
      ok <- br > 0
      w <- rep(-Inf, sum(e))
      w[ok] <- pmin(lQ[ok] + log(br[ok]), 0)
      wF[e] <- w
      wS[e] <- ifelse(w < 0, .log1m_exp(w), -Inf)
    }
    vF[j] <- wF
    vS[j] <- wS
  }
  lF[i] <- pmin(vF, 0)
  lS[i] <- pmin(vS, 0)
  list(lcdf = lF, lccdf = lS)
}


#' @keywords internal
.lognormal_acc_lcdf <- function(t, meanlog, sigma, A) {
  .lognormal_acc_ltails(t, meanlog, sigma, A)$lcdf
}


#' @keywords internal
.lognormal_acc_lccdf <- function(t, meanlog, sigma, A) {
  .lognormal_acc_ltails(t, meanlog, sigma, A)$lccdf
}


# The Stan side of the helpers above, line for line. cogmod_lnr() appends its
# race to this (see .LNR_STAN_PRELUDE in model_lnr.R).
#' @keywords internal
.LOGNORMAL_STAN_PRELUDE <- paste0(.LOG_PHI_STAN_PRELUDE, "
// log(Phi(y + c) - Phi(y)) for c > 0, from whichever tail keeps the two terms
// from cancelling: the upper one when y > 0, where both CDFs sit near 1, the
// lower one otherwise. Taken as the larger tail plus log(1 - ratio) in log
// space, the way .lognormal_ldiff_pnorm() does it - the quotient u2 / u1 this
// replaces divides two minute numbers and lost its own accuracy long before
// either underflowed (8e-4 against central differences, against 1e-6 here).
real cogmod_lognormal_ldiff_Phi(real y, real c) {
  real hi;
  real lo;
  if (y > 0) {
    hi = cogmod_log_Phi(-y);
    lo = cogmod_log_Phi(-(y + c));
  } else {
    hi = cogmod_log_Phi(y + c);
    lo = cogmod_log_Phi(y);
  }
  return lo < hi ? hi + log1m_exp(lo - hi) : negative_infinity();
}

// Log density of the accumulator's finishing time with start-point range A.
// At A = 0 this is the LogNormal itself, at the LogNormal's cost.
real cogmod_lognormal_acc_ldens(real t, real meanlog, real sigma, real A) {
  if (A == 0) return lognormal_lpdf(t | meanlog, sigma);
  real a = (meanlog - log(t)) / sigma;
  real c = log1p(A) / sigma;
  real x = a - sigma;
  if (c < 1e-4) {
    real series = 1 - c * x / 2 + square(c) * (square(x) - 1) / 6;
    if (series <= 0) return negative_infinity();
    return lognormal_lpdf(t | meanlog, sigma) + log(log1p(A) / A) + log(series);
  }
  return -meanlog + square(sigma) / 2 + cogmod_lognormal_ldiff_Phi(x, c) - log(A);
}

// [log F, log S] of the accumulator's finishing time, each computed directly
// on the side where it is the small one. See .lognormal_acc_ltails().
vector cogmod_lognormal_acc_ltails(real t, real meanlog, real sigma, real A) {
  real a = (meanlog - log(t)) / sigma;
  // At A = 0 this is the plain LogNormal, whose two tails are Phi(-a) and
  // Phi(a). Written that way rather than as lognormal_lcdf()/lognormal_lccdf(),
  // which are erfc alone and so reach log(0) with non-finite partials around
  // |a| = 38 - see cogmod_log_Phi() above for what that costs.
  if (A == 0) return [cogmod_log_Phi(-a), cogmod_log_Phi(a)]';
  real c = log1p(A) / sigma;
  if (c < 1e-4) {
    real r = log1p(A) / A;
    real corr = c * r * (0.5 + c * (2 * sigma - a) / 6);
    real lPa = cogmod_log_Phi(a);
    real lQa = cogmod_log_Phi(-a);
    real lphi = std_normal_lpdf(a);
    real lS = fmin(lPa + log1p(corr * exp(lphi - lPa)), 0);
    real dF = 1 - corr * exp(lphi - lQa);
    real lF = dF > 0 ? fmin(lQa + log(dF), 0) : negative_infinity();
    return [lF, lS]';
  }
  real lA = log(A);
  real lD1 = cogmod_lognormal_ldiff_Phi(a, c) - lA;
  real lD2 = -meanlog + square(sigma) / 2 + log(t)
             + cogmod_lognormal_ldiff_Phi(a - sigma, c) - lA;
  if (a < 0) {
    real lP = cogmod_log_Phi(a + c);
    real br = 1 + exp(lD1 - lP) - exp(lD2 - lP);
    if (br <= 0) return [0, negative_infinity()]';
    real lS = fmin(lP + log(br), 0);
    return [lS < 0 ? log1m_exp(lS) : negative_infinity(), lS]';
  }
  real lQ = cogmod_log_Phi(-a);
  real R = exp(cogmod_log_Phi(-a - c) - lQ);
  real br = R - (1 - R) / A + exp(lD2 - lQ);
  if (br <= 0) return [negative_infinity(), 0]';
  real lF = fmin(lQ + log(br), 0);
  return [lF, lF < 0 ? log1m_exp(lF) : negative_infinity()]';
}

// Above A = 0 the two tails share lD1 and lD2, so building the pair and taking
// one of them is the cheap way round. At A = 0 they share nothing - each is a
// single cogmod_log_Phi() of the same standardized time - and the pair would
// put a whole discarded tail on the autodiff tape for every observation.
// cogmod_lnr() reads the survival alone, once per trial per loser, so that is
// the hot path of the family.
real cogmod_lognormal_acc_logcdf(real t, real meanlog, real sigma, real A) {
  if (A == 0) return cogmod_log_Phi((log(t) - meanlog) / sigma);
  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[1];
}

real cogmod_lognormal_acc_logsurv(real t, real meanlog, real sigma, real A) {
  if (A == 0) return cogmod_log_Phi((meanlog - log(t)) / sigma);
  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[2];
}
")


# Log-Gamma helpers -------------------------------------------------------

# The two functions the Log-Gamma decision density needs. They are pulled out
# of the lpdf as a prelude for the same reason as the LBA's: the density itself
# is one expression in the registry, and the maths that makes it stable near
# shape = 0 does not fit in one.
#' @keywords internal
.LOGGAMMA_STAN_PRELUDE <- "
// Normalizing constant of the standardized log-gamma with shape k = 1 / shape^2,
// with its leading k already cancelled against the equal and opposite k hidden
// in the kernel below. Doing that cancellation on paper is what keeps the
// density accurate near shape = 0, where both terms reach 1e8 by |shape| = 1e-4.
// Below |shape| = 1e-3 the Stirling series is used, where 1 / (12k) = shape^2 / 12; the
// two branches agree to about 1e-9 where they meet, so the density is smooth
// through the switch and finite at shape = 0 exactly.
real cogmod_loggamma_lconst(real shape) {
  if (abs(shape) < 1e-3) {
    return -0.5 * log(2 * pi()) - square(shape) / 12 + pow(shape, 6) / 360;
  }
  real k = inv_square(shape);
  return (k - 0.5) * log(k) - lgamma(k) - k;
}

// (shape * w - expm1(shape * w)) / shape^2, the log-gamma kernel with the same k removed.
// At shape = 0 this is exactly -w^2 / 2, the Normal kernel, which is why shape = 0
// gives back the LogNormal rather than approaching it.
real cogmod_loggamma_lkernel(real shape, real w) {
  real u = shape * w;
  if (abs(u) < 1e-6) {
    return -0.5 * square(w) - shape * pow(w, 3) / 6 - square(shape) * pow(w, 4) / 24;
  }
  return (u - expm1(u)) / square(shape);
}
"


# The outlier component ---------------------------------------------------

# Shared by every shifted family, so it lives here rather than with any one of
# them.

# The outlier component is a **half Normal** on [0, Inf) with scale
# `.POUTLIER_SCALE`, in seconds.
#
# Flat at the origin (zero derivative), which is the property that matters
# most: the fastest responses are the ones least plausibly decisions, so the
# component must not thin out there. A LogNormal or Gamma vanishes at 0 and an
# Exponential peaks there with maximal slope; all three encode a claim nobody
# would defend, that anticipations are far likelier at 150 ms than at 3 ms.
#
# Two choices about it, made for the same reason.
#
# It is a half Normal and not a half Student-t with 3 degrees of freedom. That
# tail is heavier than every decision density in the package bar the
# drift-variability Wald, so far-out slow responses end up better explained by
# the outlier component than by the model: against a shifted LogNormal at
# poutlier = 2%, a 5 s response is attributed to a half-t outlier component
# with probability 0.86, the crossover sits at 3.86 s, and `ndt` is pulled up
# behind it. A Gaussian never gets there - the same responsibility is 0.000 out
# to 30 s - and it costs nothing in the region the component actually exists
# for, because exp(-x^2 / 2s^2) kills the far tail at *any* scale, so flatness
# near zero and tail weight are not traded against each other. The slow tail
# is the decision family's own business, which is what `shape` in
# cogmod_loggamma() and `sigmadrift` in cogmod_invgaussian() are for.
#
# The scale is a constant rather than a user-supplied one in the unit of the
# data. Such an argument would make the likelihood equivariant to that unit,
# but the equivariance would be fictional end to end: cogmod_priors() could
# shift only the `ndt` prior with it, while the `sigmandt` priors, the
# `sigmadrift` prior and the `mu` priors are all in seconds outright - and
# cogmod_priors() is not optional. The package works in seconds; the constant
# says so.
#
# 0.2 s is where it sits because the component has to stay flat across the
# whole range `ndt` plausibly occupies, which the `ndt` prior puts at 0.20-0.45
# s. It holds 76% of its peak density at 0.15 s and 46% at 0.25 s, against 85%
# and 66% for a half-t of scale 0.3; 68% of its mass falls below 0.2 s and 92%
# below 0.35 s. A contaminant landing just under a late `ndt` of 0.42 s still
# has a log-density of -4.6 where that half-t gives -4.0, so nothing a heavier
# tail would cover is starved.
.POUTLIER_SCALE <- 0.2

# Density of the half-Normal component. Shape is preserved, so this works on
# draws x observations matrices.
#' @keywords internal
.dcontam <- function(x, log = FALSE) {
  ld <- log(2) + stats::dnorm(x, 0, .POUTLIER_SCALE, log = TRUE)
  ld[x <= 0] <- -Inf
  if (log) ld else exp(ld)
}

# Mean of the half-Normal, s * sqrt(2 / pi).
#' @keywords internal
.mcontam <- function() {
  .POUTLIER_SCALE * sqrt(2 / pi)
}


# Resolve whether predictions should include the outlier component. An explicit
# argument wins; otherwise the flag carried on the family is used (set by
# cogmod_lognormal() or with_outliers()), defaulting to excluding it.
#' @keywords internal
.predict_outliers <- function(predict_outliers, prep) {
  if (!is.null(predict_outliers)) {
    return(isTRUE(predict_outliers))
  }
  isTRUE(prep$family$predict_outliers)
}


# Switching the outlier component on and off ------------------------------

# User-facing, and shared by every family built on the mixture, so they belong
# with the mixture rather than with the LogNormal they were first written for.

#' Include or exclude the outlier component in predictions
#'
#' @description
#' Switches the `predict_outliers` flag on a model fitted with [cogmod_lognormal()]
#' or [cogmod_loggamma()], controlling whether `posterior_predict()` and `posterior_epred()` describe the
#' fitted mixture or the decision process alone.
#'
#' Predictions **exclude** the outlier component by default, because for almost
#' every downstream use it is a nuisance: it pulls expected values toward its own
#' mean (0.16 s) and adds a
#' spike of implausibly fast draws to posterior predictive samples. It is also a deliberately fixed regularizer rather than a
#' claim about how guesses are distributed, so simulating from it means
#' simulating from something the model does not assert.
#'
#' `with_outliers()` restores the mixture. The main reason to want it is
#' `brms::pp_check()`: on untrimmed data the decision-only predictive has no fast
#' spike to match the one in the data, which reads as misfit. Use
#' `pp_check(with_outliers(m))` for a like-for-like check.
#'
#' The flag is stored **on the model** rather than passed as an argument, because
#' `brms` and the packages built on it (`insight`, `modelbased`,
#' `marginaleffects`, `emmeans`) do not forward extra arguments down to a custom
#' family's prediction methods - `posterior_epred()` reaches the family method
#' with `prep` and nothing else. Carrying it on the object is what makes it work
#' through all of them. The same flag can be set up front with
#' `cogmod_lognormal(predict_outliers = TRUE)`.
#'
#' `log_lik()` is unaffected and has no equivalent switch: the likelihood *is*
#' the mixture, and dropping a component from it would not be a different summary
#' of the same model but a different model. One consequence worth knowing is that
#' `posterior_predict()` and `log_lik()` do not describe the same distribution by
#' default. This also desyncs `loo_pit()`, `loo_predict()` and `bayes_R2()` from
#' `loo()`, not just hand-rolled checks - anything that compares a simulated
#' replicate against the likelihood should be run on `with_outliers()`.
#'
#' @param object A `brmsfit` fitted with [cogmod_lognormal()], [cogmod_loggamma()]
#'   or any other family built on the outlier mixture - see the *Supported
#'   families* section of [cogmod_priors()] for the full list, which includes
#'   the choice-and-RT families such as [cogmod_lnr()].
#'
#' @return The model, with the flag set. The fit itself is untouched - only how
#'   predictions are summarised changes.
#'
#' @examples
#' \donttest{
#' # Fitting needs cmdstanr, which lives outside CRAN - see the package website.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   df <- data.frame(
#'     RT = rcogmod_lognormal(200, ndt = 0.3, poutlier = 0.05),
#'     Condition = rep(c("A", "B"), each = 100)
#'   )
#'   f <- brms::bf(RT ~ Condition, ndt ~ 1, poutlier ~ 1,
#'     family = cogmod_lognormal()
#'   )
#'   m <- brms::brm(f,
#'     data = df, stanvars = cogmod_stanvars(f),
#'     prior = cogmod_priors(f, df), init = cogmod_inits(f, df),
#'     backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
#'   )
#'
#'   # the decision process alone - the default, everywhere downstream
#'   head(brms::posterior_epred(m)[, 1])
#'
#'   # the fitted mixture, e.g. for a like-for-like predictive check
#'   m2 <- with_outliers(m)
#'   head(brms::posterior_epred(m2)[, 1])
#'
#'   without_outliers(m2) # back to the default
#' }
#' }
#'
#' @export
with_outliers <- function(object) {
  .set_predict_outliers(object, TRUE)
}


#' @rdname with_outliers
#' @export
without_outliers <- function(object) {
  .set_predict_outliers(object, FALSE)
}


#' @keywords internal
.set_predict_outliers <- function(object, value) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.", call. = FALSE)
  }
  if (!isTRUE(object$family$name %in% .OUTLIER_FAMILIES)) {
    stop("`object` must have been fitted with one of the families carrying a ",
         "`poutlier` parameter: ",
         paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
         call. = FALSE)
  }
  # brms rebuilds prep$family from object$formula$family, so the flag has to
  # live there to survive prepare_predictions(). object$family is set too, so
  # that inspecting the fit tells the same story.
  object$formula$family$predict_outliers <- value
  object$family$predict_outliers <- value
  object
}


#' Per-trial outlier probabilities
#'
#' @description
#' Posterior probability that each response was generated by the outlier
#' component rather than by the decision process, for a model fitted with any
#' family built on the outlier mixture - [cogmod_lognormal()], [cogmod_loggamma()],
#' [cogmod_lnr()] and the rest listed in the *Supported families* section of
#' [cogmod_priors()]. This is the mixture *responsibility*
#' `poutlier * g(rt) / (poutlier * g(rt) + (1 - poutlier) * f(rt - ndt))`,
#' averaged over posterior draws.
#'
#' Responses faster than `ndt` come out at 1, responses in the heart of the
#' distribution near 0, and responses in either tail somewhere in between - the
#' model discriminates by evidence rather than by a cutoff. A response in the
#' middle can still be an outlier; a low probability means the data cannot tell,
#' not that the trial is clean.
#'
#' Averaging the responsibility over draws gives `P(trial i came from the
#' outlier component | data)` directly, so the posterior mean *is* the quantity
#' of interest and there is no interval to report alongside it. Pass
#' `summary = FALSE` for the raw draws if you need the spread.
#'
#' @param object A `brmsfit` fitted with [cogmod_lognormal()], [cogmod_loggamma()]
#'   or any other family built on the outlier mixture - see the *Supported
#'   families* section of [cogmod_priors()].
#' @param summary Logical; if `TRUE` (default) returns a data frame with one row
#'   per observation. If `FALSE`, returns the full draws x observations matrix.
#'
#' @return A data frame with columns `rt` and `p_outlier`, in the order the
#'   observations appear in the model frame, or a draws x observations matrix if
#'   `summary = FALSE`.
#'
#' @examples
#' \donttest{
#' # Fitting needs cmdstanr, which lives outside CRAN - see the package website.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   df <- data.frame(RT = rcogmod_lognormal(200, ndt = 0.3, poutlier = 0.05))
#'   f <- brms::bf(RT ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_lognormal())
#'   m <- brms::brm(f,
#'     data = df, stanvars = cogmod_stanvars(f),
#'     prior = cogmod_priors(f, df), init = cogmod_inits(f, df),
#'     backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
#'   )
#'   head(p_outlier(m))
#' }
#' }
#'
#' @export
p_outlier <- function(object, summary = TRUE) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.")
  }
  prep <- brms::prepare_predictions(object)
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in the fitted model.")
  }
  y <- prep$data$Y
  n <- length(y)

  # Models fitted before the family carried a name fall back to the LogNormal,
  # which is the only family this ever supported.
  fam <- prep$family$name
  if (is.null(fam)) fam <- "cogmod_lognormal"
  spec <- .mixture_spec(fam)

  mu <- brms::get_dpar(prep, "mu")
  n_draws <- if (is.matrix(mu)) nrow(mu) else brms::ndraws(object)
  as_mat <- function(x) {
    if (is.matrix(x)) x else matrix(x, nrow = n_draws, ncol = n)
  }

  # Whatever decision dpars this family happens to have - the registry knows,
  # so nothing here needs to name them.
  dpars <- lapply(spec$dpars, function(d) as_mat(brms::get_dpar(prep, d)))
  names(dpars) <- spec$dpars
  ndt <- as_mat(brms::get_dpar(prep, "ndt"))
  poutlier <- as_mat(brms::get_dpar(prep, "poutlier"))
  Y <- matrix(y, nrow = n_draws, ncol = n, byrow = TRUE)

  # Done in logs throughout: both components underflow to exactly 0 in the tails
  # on the natural scale, which turns the ratio into 0/0 well before either
  # component is genuinely negligible.
  #
  # A choice family's contaminant has to produce a choice as well as a time, so
  # the outlier density carries the 1 / K that spreads the guess over the
  # response options, and the decision density needs the observed choice.
  if (is.null(spec$K)) {
    lp_out <- .dcontam(Y, log = TRUE)
    lp_dec <- .ldec(fam, Y - ndt, dpars)
  } else {
    resp <- matrix(.dec_from_prep(prep), nrow = n_draws, ncol = n,
                   byrow = TRUE)
    lp_out <- .lout_choice(Y, spec$K)
    lp_dec <- .ldec_choice(fam, Y - ndt, resp, dpars)
  }

  log_num <- log(poutlier) + lp_out
  log_den <- .log_mix(poutlier, lp_out, lp_dec)
  r <- exp(log_num - log_den)
  # Both components vanish (log_den = -Inf): nothing but the outlier could have
  # produced the response.
  r[!is.finite(r)] <- 1

  if (!summary) return(r)

  data.frame(rt = y, p_outlier = colMeans(r))
}
