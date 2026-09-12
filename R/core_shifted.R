# Shared machinery for the shifted RT families
# =============================================
#
# Every family in this group has the same structure: a decision-time
# distribution shifted by `ndt`, mixed with a half Student-t outlier component
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
    dpars = c("mu", "sigma"),
    links = c("identity", "softplus"),
    lb = c(NA, 0), ub = c(NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = "lognormal_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) stats::dlnorm(t, p$mu, p$sigma, log = TRUE),
    rng = function(n, p) stats::rlnorm(n, p$mu, p$sigma),
    mean = function(p) exp(p$mu + p$sigma^2 / 2),
    init = list(mu = -0.7, sigma = 0.5),
    dpar_doc = c(
      "mu: mean of the decision time on the log scale (meanlog).",
      "sigma: SD of the decision time on the log scale (> 0)."
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
    dpars = c("mu", "boundary", "sigmadrift"),
    links = c("softplus", "softplus", "softplus"),
    lb = c(0, 0, 0), ub = c(NA, NA, NA),
    # A zero drift SD is the plain Wald, so that bound is closed - unlike the
    # drift and the threshold, which a Wald needs strictly positive.
    lb_open = c(TRUE, TRUE, FALSE),
    stan_check = "mu <= 0 || boundary <= 0 || sigmadrift < 0",
    stan_dens = "cogmod_invgaussian_decision_lpdf(t_adj | mu, boundary, sigmadrift)",
    prelude = ".WALD_STAN_PRELUDE",
    ldens = function(t, p) .dwald_raw(t, p$mu, p$boundary, p$sigmadrift),
    rng = function(n, p) .rwald_raw(n, p$mu, p$boundary, p$sigmadrift),
    # E[T] = boundary / mu holds only while the drift is fixed. Once it varies
    # the density decays as t^-2 - drifts arbitrarily close to zero take
    # arbitrarily long - so the mean diverges, exactly as it does for
    # cogmod_lba1(). Inf is what posterior_epred() then returns.
    mean = function(p) ifelse(p$sigmadrift > 0, Inf, p$boundary / p$mu),
    init = list(mu = 3, boundary = 0.5, sigmadrift = 0.5),
    # Same flat direction as cogmod_lba1()'s `sigmabias` and cogmod_ddm()'s own
    # `sigmadrift`: the floor at zero is only reached by the softplus link at
    # minus infinity, and the likelihood stops changing well before then, which
    # is an improper posterior wherever brms would leave it flat. There is a
    # second reason here - scaling (mu, boundary, sigmadrift) up by a common
    # factor sends the Wald to the reciprocal-normal (LATER) limit, so large
    # values of all three describe very nearly the same distribution. The prior
    # is what keeps the sampler out of both.
    prior = list(
      sigmadrift = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                     slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "mu: drift rate, the average speed of evidence accumulation (> 0).",
      "boundary: decision threshold, the evidence needed to respond (> 0).",
      paste("sigmadrift: between-trial SD of the drift rate (>= 0), which is",
            "drawn from a"),
      "   Normal(mu, sigmadrift) truncated at zero. 0 is the plain Wald."
    ),
    note = c(
      paste("sigmadrift and poutlier both fatten the right tail, and they are",
            "only weakly"),
      paste("distinguishable: on 2000 simulated trials at mu = 3, boundary =",
            "0.5, sigmadrift ="),
      paste("0.8, estimating sigmadrift buys about 2 log-likelihood units over",
            "fixing it at"),
      "zero. Fix it (`sigmadrift = 0` in bf()) unless the design can identify it."
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
    # density itself rather than an approximation to it. .lba_surv_raw()
    # likewise collapses to Phi(z1), the survival of a deterministic start
    # point, which is what cogmod_lba2() needs.
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

# Generates the `<name>_lpdf` Stan function: the same mixture skeleton for every
# family, with only the decision density swapped in. The outlier component's
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
  # posterior unchanged. The same reasoning applied to the half Student-t this
  # replaced; a half Normal simply leaves less to fold.
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
  sprintf(
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
    // Student-t onto [0, Inf).
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
  # It replaces a bare `pmax(t, 1e-300)`, which covered the non-positive case
  # but let NA through to a branch that cannot take it. See `.dens_mask()`.
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

  ll <- do.call(.dshifted, c(
    list(name = name, x = rep(y, length.out = n_draws), ndt = ndt,
         poutlier = poutlier, log = TRUE),
    dec
  ))
  ll[is.na(ll)] <- -Inf
  ll
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
#' @keywords internal
.dwald_raw <- function(t, drift, boundary, sigmadrift = 0) {
  s2 <- sigmadrift^2
  D <- 1 + s2 * t
  ld <- log(boundary) - 0.5 * (log(2 * pi) + 3 * log(t) + log(D)) -
    (boundary - drift * t)^2 / (2 * t * D)

  # The division below is by `s`, which stands in as 1 wherever sigmadrift is 0
  # so that the truncation term stays finite; multiplying by the indicator then
  # drops it. Written this way rather than by subsetting because `t` and the
  # parameters may be draws x observations matrices, whose shape has to survive.
  s <- sigmadrift + (sigmadrift <= 0)
  ld + (sigmadrift > 0) *
    (.lpnorm_upper((boundary * s2 + drift) / (s * sqrt(D))) -
      .lpnorm_upper(drift / s))
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
# random stream of the plain Wald exactly as it was.
#' @keywords internal
.rwald_raw <- function(n, drift, boundary, sigmadrift = 0) {
  if (any(sigmadrift > 0, na.rm = TRUE)) {
    drift <- .rnorm_truncated(n, mean = drift, sd = sigmadrift, lower = 0)
  }
  ig_mu <- boundary / drift
  lambda <- boundary^2
  y <- stats::rnorm(n)^2
  z <- y * (ig_mu / lambda)
  x1_over_mu <- 1 + z / 2 * (1 - sqrt(1 + 4 / z))
  u <- stats::runif(n)
  ig_mu * ifelse(u < 1 / (1 + x1_over_mu), x1_over_mu, 1 / x1_over_mu)
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
# reason for it: at sigmadrift = 0 the two normal CDFs are 0/0 rather than 1.
# Both are taken through the lower tail (`log1m_exp(std_normal_lcdf(-z))`)
# because with a positive drift and threshold both arguments are positive, which
# is where Phi itself saturates.
#' @keywords internal
.WALD_STAN_PRELUDE <- "
// Log density of the Wald decision time (no shift), with the drift rate drawn
// once per trial from Normal(mu, sigmadrift) truncated at zero. sigmadrift = 0
// is the plain Wald: an inverse Gaussian with mean boundary / mu and shape
// boundary^2.
real cogmod_invgaussian_decision_lpdf(real t, real mu, real boundary, real sigmadrift) {
  real base = log(boundary) - 0.5 * (log(2 * pi()) + 3 * log(t));
  if (sigmadrift <= 0) {
    return base - square(boundary - mu * t) / (2 * t);
  }
  real s2 = square(sigmadrift);
  real D = 1 + s2 * t;
  return base - 0.5 * log(D) - square(boundary - mu * t) / (2 * t * D)
    + log1m_exp(std_normal_lcdf(-(boundary * s2 + mu) / (sigmadrift * sqrt(D))))
    - log1m_exp(std_normal_lcdf(-mu / sigmadrift));
}
"


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
# the two implementations are generated from one table.
#' @keywords internal
.EXWALD_STAN_PRELUDE <- local({
  num <- function(v) formatC(v, format = "g", digits = 17, width = 1)
  sprintf("
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

// log of the Wald CDF with drift v >= 0 and threshold a > 0, at t > 0. The
// exp(2 a v) factor overflows on its own long before the product it belongs to
// stops being representable, so it is folded into the exponent instead.
real cogmod_exwald_lwaldcdf(real t, real v, real a) {
  real st = sqrt(t);
  return log_sum_exp(
    std_normal_lcdf((v * t - a) / st),
    2 * a * v + std_normal_lcdf(-(v * t + a) / st)
  );
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
      + cogmod_exwald_lwaldcdf(t, k, boundary);
  }

  real st = sqrt(t);
  real rw = cogmod_re_faddeeva(sqrt(-k2) * st / sqrt2(), boundary / (st * sqrt2()));
  if (rw <= 0) return negative_infinity();
  return log(g) - square(boundary - mu * t) / (2 * t) + log(rw);
}
",
    .WEIDEMAN$N, num(.WEIDEMAN$L), .WEIDEMAN$N,
    paste(num(.WEIDEMAN$a), collapse = ", "), .WEIDEMAN$N)
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


# CDF of the half Student-t outlier component.
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


# Survival of one accumulator at `t`: the probability it has NOT yet reached the
# threshold. Only cogmod_lba2() needs it - a race scores the winner's density
# against the loser's survival - but it is the same integral over the same
# start-point range, so it belongs beside the density.
#
# Writing b - A - v t = z1 * st and b - v t = z2 * st turns the textbook CDF
#
#   F(t) = 1 + ((b - A - vt)/A) Phi(z1) - ((b - vt)/A) Phi(z2)
#            + (t s / A) (phi(z1) - phi(z2))
#
# into `1 - (g(z2) - g(z1)) / delta` with `g(z) = z Phi(z) + phi(z)`, so the
# survival is that quotient directly rather than `1 - F`, which cancels to
# nothing whenever the accumulator is unlikely to have finished. The quotient
# cancels in its own right as delta -> 0, hence the same Taylor branch as the
# density, the derivatives being Phi, phi and -z phi.
#' @keywords internal
.lba_surv_raw <- function(drift, sigma, st, z1, delta) {
  n <- max(length(drift), length(sigma), length(st), length(z1), length(delta))
  z1 <- rep_len(z1, n)
  delta <- rep_len(delta, n)

  phi1 <- stats::dnorm(z1)
  out <- numeric(n)

  small <- delta < 1e-4
  if (any(small)) {
    d <- delta[small]
    z <- z1[small]
    out[small] <- stats::pnorm(z) + (d / 2) * phi1[small] -
      (d^2 / 6) * z * phi1[small]
  }
  if (any(!small)) {
    i <- !small
    z2 <- z1[i] + delta[i]
    g <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
    out[i] <- (g(z2) - g(z1[i])) / delta[i]
  }
  pmin(pmax(out, 0), 1)
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

# Stan counterpart of .lba_dens_over_A(), .lba_surv_raw() and .dlba1_raw(). Same
# branches, and the same reason for them: both differences cancel as the
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

// Probability that an accumulator has NOT finished by t, as
// (g(z2) - g(z1)) / delta with g(z) = z Phi(z) + phi(z). Taken directly rather
// than as 1 - CDF, which cancels to nothing whenever the accumulator is
// unlikely to have finished; the Taylor branch is for the other cancellation,
// the one in the quotient as delta -> 0.
real cogmod_lba_surv(real z1, real delta) {
  real phi1 = exp(-0.5 * square(z1)) * 0.3989422804014327;
  real out;
  if (delta < 1e-4) {
    out = Phi(z1) + (delta / 2) * phi1 - (square(delta) / 6) * z1 * phi1;
  } else {
    real z2 = z1 + delta;
    real phi2 = exp(-0.5 * square(z2)) * 0.3989422804014327;
    out = ((z2 * Phi(z2) + phi2) - (z1 * Phi(z1) + phi1)) / delta;
  }
  return fmin(fmax(out, 0), 1);
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
# Two things about it were different before 0.2.1, and both changed for the
# same reason.
#
# It was a half Student-t with 3 degrees of freedom. That tail is heavier than
# every decision density in the package bar the drift-variability Wald, so
# far-out slow responses ended up better explained by the outlier component
# than by the model: against a shifted LogNormal at poutlier = 2%, a 5 s
# response was attributed to the outlier component with probability 0.86, the
# crossover sat at 3.86 s, and `ndt` was pulled up behind it. A Gaussian never
# gets there - the same responsibility is 0.000 out to 30 s - and it costs
# nothing in the region the component actually exists for, because
# exp(-x^2 / 2s^2) kills the far tail at *any* scale, so flatness near zero and
# tail weight are no longer traded against each other. The slow tail is now the
# decision family's own business, which is what `shape` in cogmod_loggamma()
# and `sigmadrift` in cogmod_invgaussian() are for.
#
# Its scale was `minrt`, a user-supplied constant in the unit of the data,
# which made the likelihood equivariant to that unit. That equivariance was
# already fictional end to end: cogmod_priors() shifted only the `ndt` prior
# with it, while cogmod_ddm()'s `sigmandt` prior, cogmod_invgaussian()'s
# `sigmadrift` prior and the `mu` priors are all in seconds outright - and
# cogmod_priors() is not optional. The package works in seconds; the scale is
# now a constant that says so, and `minrt` is gone from every signature.
#
# 0.2 s is where it sits because the component has to stay flat across the
# whole range `ndt` plausibly occupies, which the `ndt` prior puts at 0.20-0.45
# s. It holds 76% of its peak density at 0.15 s and 46% at 0.25 s, against 85%
# and 66% for the half-t it replaces; 68% of its mass falls below 0.2 s and 92%
# below 0.35 s. A contaminant landing just under a late `ndt` of 0.42 s still
# has a log-density of -4.6 where the half-t gave -4.0, so nothing that used to
# be covered is starved now.
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
