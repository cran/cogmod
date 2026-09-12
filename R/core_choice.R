# Shared machinery for the shifted choice + RT families
# =====================================================
#
# The counterpart of core_shifted.R for the models that produce a *choice* as
# well as a reaction time. The structure is the same - a decision-time
# distribution shifted by `ndt`, mixed with a half Student-t outlier component
# of weight `poutlier` and a fixed scale - with one genuine difference, which is
# why these cannot simply be added to `.SHIFTED`.
#
# A choice model has a *defective* density per response option `k`, with
# `sum_k integral f_k(t) dt = 1`. A contaminant process therefore has to produce
# a choice as well as an RT. The construction used here is a guess: the choice
# is uniform over the `K` options and the RT comes from the same half Student-t
# as in the RT-only families.
#
#   f(t, k) = poutlier * (1 / K) * g(t) + (1 - poutlier) * f_k(t - ndt)
#
# The `1 / K` is what keeps this normalized:
#
#   integral sum_k [ p*(1/K)*g(t) + (1-p)*f_k(t-ndt) ] dt
#     = p * sum_k (1/K) * integral g  +  (1-p) * sum_k integral f_k
#     = p * 1 + (1-p) * 1 = 1
#
# Without it the total comes to `1 + p` rather than `1`, which is silent: every
# individual density is still positive and still looks reasonable.
#
# See ?rcogmod_lognormal for the account of why `ndt` is expressed directly
# rather than as a fraction of an observed minimum, what the outlier component
# is for, and why the outlier component's scale is a constant rather than a
# dpar. The outlier component itself (`.dcontam()`, `.POUTLIER_SCALE`,
# `.predict_outliers()`, `with_outliers()`, `p_outlier()`) is shared with the
# RT-only families and lives in core_shifted.R.


# The registry ------------------------------------------------------------

# Each entry describes the *decision* component only.
#
#  dpars/links/lb/ub : the distributional parameters, before ndt and poutlier
#  lb_open           : optional, one flag per dpar - is the lower bound an open
#                      one? Defaults to TRUE, which is what every bound in
#                      .SHIFTED is. cogmod_rdm() is the exception: a drift of
#                      exactly zero is a legitimate value there
#  K                 : number of response options
#  vars              : the brms `vars` the family needs for the response index
#  stan_check        : Stan expression that is TRUE for invalid parameters
#  stan_dens         : Stan expression for the decision log-density at `t_adj`,
#                      which may branch on `dec`
#  ldens             : R log-density of the decision component, at t > 0, given
#                      the response index `k`
#  rng               : R sampler; returns list(rt, response), unshifted. A bare
#                      list, not a data frame: `.rchoice()` only reads the two
#                      elements off it, and building a data frame per call is
#                      the single most expensive thing in `posterior_predict()`
#  init              : natural-scale starting values, for cogmod_inits()
#  prior             : optional, per-dpar priors for cogmod_priors() to fill,
#                      where the likelihood has a flat direction brms would
#                      otherwise leave improper. Same shape as in .SHIFTED
#  label             : human-readable name, used in the generated Stan comments
#' @keywords internal
.CHOICE <- list(
  cogmod_lnr = list(
    # `mu` is nuzero: brms requires the first dpar of a custom family to be
    # called `mu`. Both nu are minus the log-scale mean, so larger = faster.
    dpars = c("mu", "nuone", "sigmazero", "sigmaone"),
    links = c("identity", "identity", "softplus", "softplus"),
    lb = c(NA, NA, 0, 0), ub = c(NA, NA, NA, NA),
    K = 2L,
    vars = "dec[n]",
    stan_check = "sigmazero <= 0 || sigmaone <= 0",
    stan_dens = paste0(
      "dec == 0\n",
      "      ? lognormal_lpdf(t_adj | -mu, sigmazero)",
      " + lognormal_lccdf(t_adj | -nuone, sigmaone)\n",
      "      : lognormal_lpdf(t_adj | -nuone, sigmaone)",
      " + lognormal_lccdf(t_adj | -mu, sigmazero)"
    ),
    ldens = function(t, k, p) {
      win_ml <- ifelse(k == 0, -p$mu, -p$nuone)
      win_sd <- ifelse(k == 0, p$sigmazero, p$sigmaone)
      los_ml <- ifelse(k == 0, -p$nuone, -p$mu)
      los_sd <- ifelse(k == 0, p$sigmaone, p$sigmazero)
      # lower.tail = FALSE gives the log-survival of the loser directly, which
      # is stable where log(1 - exp(log F)) is not.
      stats::dlnorm(t, meanlog = win_ml, sdlog = win_sd, log = TRUE) +
        stats::plnorm(t, meanlog = los_ml, sdlog = los_sd,
                      lower.tail = FALSE, log.p = TRUE)
    },
    rng = function(n, p) {
      d0 <- stats::rlnorm(n, meanlog = -p$mu, sdlog = p$sigmazero)
      d1 <- stats::rlnorm(n, meanlog = -p$nuone, sdlog = p$sigmaone)
      list(rt = pmin(d0, d1), response = as.numeric(d0 >= d1))
    },
    init = list(mu = 0.7, nuone = 0.7, sigmazero = 0.5, sigmaone = 0.5),
    # Push an accumulator's rate down far enough and it stops finishing first
    # ever; the density then depends on it only through the loser's survival
    # term, which has already saturated at 1. Beyond about nuone = -6 the
    # log-likelihood is *exactly* constant - and so is that accumulator's sigma,
    # which has nothing left to act on. Both links reach that plateau at a
    # finite distance (identity) or at minus infinity (softplus), so a flat
    # prior over it is an improper posterior, exactly as for `ndt` and
    # `poutlier`. The outlier component makes it worse rather than better: it
    # floors the trials the retreating accumulator can no longer explain, so the
    # plateau is reached even when both responses are observed. See
    # ?rcogmod_lnr.
    #
    # `mu` - nuzero - has the mirror-image plateau, but it is the response's own
    # intercept, so brms already gives it a proper student_t default and it is
    # left alone. Model a rarely-chosen option and it is worth mirroring the
    # `nuone` prior onto it by hand.
    prior = list(
      nuone = c(link = "normal(0.7, 1.5)", nat = "normal(0.7, 1.5)",
                slope = "normal(0, 0.5)"),
      sigmazero = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                    slope = "normal(0, 0.5)"),
      sigmaone = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                   slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "dec: the observed choice, 0 or 1.",
      paste("mu: nuzero, the processing speed of accumulator 0",
            "(meanlog = -mu, so larger is faster)."),
      "nuone: the same for accumulator 1.",
      "sigmazero, sigmaone: the log-scale SD of each accumulator (> 0)."
    ),
    label = "Log-Normal Race"
  ),
  cogmod_rdm = list(
    # Two Wald accumulators race to a common threshold b = boundary + sigmabias,
    # each starting from z ~ Uniform(0, sigmabias). `mu` is the drift rate of
    # accumulator 0: brms requires the first dpar of a custom family to be
    # called `mu`.
    dpars = c("mu", "driftone", "sigmabias", "boundary"),
    links = c("softplus", "softplus", "softplus", "softplus"),
    lb = c(0, 0, 0, 0), ub = c(NA, NA, NA, NA),
    # A drift of exactly zero is legitimate here, unlike anywhere else in either
    # registry: driftless Brownian motion still reaches a positive level with
    # probability one, so a zero-drift accumulator is slow rather than absent and
    # can still win the race. The Stan check below rejects the same set.
    #
    # A start-point range of exactly zero is legitimate too, and for the same
    # kind of reason: it is a model rather than an invalid parameter. Both
    # accumulators then start at 0 on every trial and the race is a plain Wald
    # one - Tillman et al. (2020) equation 2, which is what the density's
    # midpoint branch already collapses to at A = 0. cogmod_lba2() and
    # cogmod_lba1() have always allowed it; the RDM excluding it was an
    # inconsistency, not a decision. Note this does not make the
    # sigmabias -> 0 ridge below any less flat: the softplus link still reaches
    # zero only at minus infinity, so cogmod_priors() still has to fence it.
    lb_open = c(FALSE, FALSE, FALSE, TRUE),
    K = 2L,
    vars = "dec[n]",
    stan_check = "boundary <= 0 || sigmabias < 0 || mu < 0 || driftone < 0",
    stan_dens = paste0(
      "cogmod_rdm_wald_ldens(t_adj, dec == 0 ? mu : driftone,",
      " boundary, sigmabias)\n",
      "      + cogmod_rdm_wald_lsurv(t_adj, dec == 0 ? driftone : mu,",
      " boundary, sigmabias)"
    ),
    prelude = ".RDM_STAN_PRELUDE",
    ldens = function(t, k, p) .rdm_ldens(t, k, p),
    rng = function(n, p) .rdm_rng(n, p),
    init = list(mu = 3, driftone = 3, sigmabias = 0.3, boundary = 0.5),
    # `sigmabias` and `boundary` enter the threshold only through their sum,
    # b = boundary + sigmabias, and trade off along a very flat ridge: on
    # simulated data with 4000 trials the profile log-likelihood moves by only
    # about 5 units as sigmabias ranges from 0 to half the threshold, with
    # boundary sliding from 0.59 to 0.46 to pay for it. Under flat priors the
    # sampler wanders down the sigmabias -> 0 edge - the plain Wald race - and
    # produces divergent transitions. Same failure, and the same fix, as
    # cogmod_lba1(), which shares this parameterization. See ?rcogmod_rdm.
    #
    # `driftone` has the retreating-accumulator plateau described under
    # cogmod_lnr() above, and it is the same failure for the same reason: push
    # the error accumulator's rate down far enough and it never finishes first,
    # so the likelihood depends on it only through the loser's survival, which
    # has already saturated at 1. The softplus link reaches zero only at minus
    # infinity, so a flat prior over that direction is an improper posterior.
    # It is not hypothetical - fitting the two-participant lexical decision data
    # of vignette("decision_making"), where the Accuracy condition carries only
    # 75 errors, ran the intercept to -9.1 (a drift of 0.0001, the floor) with
    # the condition slope going to +8.9 to pay for it, and 49% of transitions
    # divergent. The LNR, which fences the same direction, samples the same data
    # without a single divergence.
    #
    # `mu` has the mirror-image plateau but is the response's own intercept, so
    # brms already gives it a proper student_t default and it is left alone -
    # again as for cogmod_lnr(). Model a rarely-chosen option and it is worth
    # mirroring the `driftone` prior onto it by hand.
    #
    # The slope here is `normal(0, 2)` rather than the `normal(0, 0.5)` the two
    # threshold parameters get: a drift sits around 3 on this link, so a genuine
    # condition effect is worth a unit or two and a tighter prior would fight
    # it. It still fences the runaway, which was an order of magnitude larger.
    prior = list(
      driftone = c(link = "normal(3, 2)", nat = "lognormal(1, 0.75)",
                   slope = "normal(0, 2)"),
      sigmabias = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                    slope = "normal(0, 0.5)"),
      boundary = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                   slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "dec: the observed choice, 0 or 1.",
      "mu: drift rate of accumulator 0 (>= 0).",
      "driftone: drift rate of accumulator 1 (>= 0).",
      "sigmabias: start-point range A; the start point is z ~ Uniform(0, sigmabias).",
      "boundary: threshold offset, so the threshold is b = boundary + sigmabias."
    ),
    note = c(
      "The decision component is carried entirely in log space. That is not",
      "stylistic: for moderate drift the loser's survival underflows to exactly",
      "zero (at drift 6 it does so by t = 4), and log(0) hands the sampler a",
      "zero gradient, which stalls it silently rather than erroring."
    ),
    label = "Racing Diffusion"
  ),
  cogmod_lba2 = list(
    # Two ballistic accumulators race to a common threshold b = boundary +
    # sigmabias from a start point z ~ Uniform(0, sigmabias), each with its own
    # normal drift. `mu` is driftzero: brms requires the first dpar of a custom
    # family to be called `mu`.
    dpars = c("mu", "driftone", "sigmazero", "sigmaone", "sigmabias",
              "boundary"),
    links = c("identity", "identity", "softplus", "softplus", "softplus",
              "softplus"),
    lb = c(NA, NA, 0, 0, 0, 0), ub = c(NA, NA, NA, NA, NA, NA),
    # Zero start-point variability is a model rather than an invalid parameter,
    # exactly as in cogmod_lba1(): both accumulators then start at 0 every
    # trial and only the drifts vary. The same two kernels cover it - the
    # density's Taylor branch at delta = 0 is the recinormal and the survival
    # collapses to Phi(z1) - so nothing here changes but the bound. See the
    # cogmod_lba1() entry in core_shifted.R for the derivation.
    lb_open = c(NA, NA, TRUE, TRUE, FALSE, TRUE),
    K = 2L,
    vars = "dec[n]",
    stan_check = paste("sigmazero <= 0 || sigmaone <= 0 || sigmabias < 0",
                       "|| boundary <= 0"),
    stan_dens = paste0(
      "dec == 0\n",
      "      ? cogmod_lba2_decision_lpdf(t_adj | mu, sigmazero,",
      " driftone, sigmaone, sigmabias, boundary)\n",
      "      : cogmod_lba2_decision_lpdf(t_adj | driftone, sigmaone,",
      " mu, sigmazero, sigmabias, boundary)"
    ),
    prelude = ".LBA2_STAN_PRELUDE",
    ldens = function(t, k, p) .lba2_ldens(t, k, p),
    rng = function(n, p) .lba2_rng(n, p),
    # As for cogmod_lba1(), with two accumulators' worth of parameters on the
    # same unitless evidence scale. Only the RATIO sigmaone / sigmazero is a
    # real quantity, which is why the convention pins sigmazero rather than
    # both. See .warn_scale_ray().
    scale_ray = c("mu", "driftone", "sigmazero", "sigmaone", "sigmabias",
                  "boundary"),
    init = list(mu = 3, driftone = 3, sigmazero = 1, sigmaone = 1,
                sigmabias = 0.5, boundary = 0.5),
    # Two flat directions here, not one.
    #
    # `sigmabias` and `boundary` share the ridge cogmod_lba1() and cogmod_rdm()
    # have: they enter only through the sum b = boundary + sigmabias.
    #
    # `sigmazero` and `sigmaone` are worse. The evidence scale of an LBA is
    # arbitrary - multiply the drifts, their SDs, the start-point range and the
    # threshold by any c > 0 and every finishing time (b - z) / v is unchanged -
    # so the likelihood is *exactly* constant along that ray, which runs to
    # infinity in both directions. Priors make the posterior proper; they do not
    # identify the scale. Fix one SD in the formula (`sigmazero = 1` in bf(),
    # the usual convention) if the individual parameters are to be interpreted
    # rather than the RT distribution they generate. See ?rcogmod_lba2.
    prior = list(
      sigmazero = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                    slope = "normal(0, 0.5)"),
      sigmaone = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                   slope = "normal(0, 0.5)"),
      sigmabias = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                    slope = "normal(0, 0.5)"),
      boundary = c(link = "normal(0, 1)", nat = "lognormal(-0.7, 0.75)",
                   slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "dec: the observed choice, 0 or 1.",
      "mu: driftzero, the mean drift rate of accumulator 0.",
      "driftone: the same for accumulator 1.",
      "sigmazero, sigmaone: between-trial SD of each drift rate (> 0).",
      "sigmabias: start-point range A; the start point is z ~ Uniform(0, sigmabias).",
      "boundary: threshold offset, so the threshold is b = boundary + sigmabias."
    ),
    note = c(
      "A normal drift can come out negative, and such an accumulator never",
      "reaches the threshold. The decision density is therefore conditioned on",
      "at least one of the two drifts being positive - the event the process",
      "itself is conditioned on, since a trial with neither is not a trial.",
      "Without that normalisation the density integrates to the probability of",
      "the event rather than to one, which at low drift rates is a long way",
      "short: 0.83 at drifts of 0.5 and 0.2 with SDs of 1.5."
    ),
    label = "two-accumulator LBA"
  ),
  cogmod_ddm = list(
    # A single diffusion between two absorbing boundaries. `mu` is the drift:
    # brms requires the first dpar of a custom family to be called `mu`. The
    # response coded 1 is the UPPER boundary and 0 the lower, following brms's
    # own wiener() family, which is why the decision density is called with the
    # drift and the starting point flipped for dec == 0.
    dpars = c("mu", "boundary", "bias", "sigmadrift", "sigmabias", "sigmandt"),
    links = c("identity", "softplus", "logit", "softplus", "logit", "log"),
    lb = c(NA, 0, 0, 0, 0, 0), ub = c(NA, NA, 1, NA, 1, NA),
    # The three between-trial variability parameters are legitimately zero -
    # that is the classic 4-parameter DDM - so their lower bounds are closed.
    lb_open = c(NA, TRUE, TRUE, FALSE, FALSE, FALSE),
    K = 2L,
    vars = "dec[n]",
    stan_check = paste(
      "boundary <= 0 || bias <= 0 || bias >= 1 || sigmadrift < 0",
      "|| sigmabias < 0 || sigmabias >= 1 || sigmandt < 0"
    ),
    stan_dens = paste0(
      "dec == 1",
      "\n      ? cogmod_ddm_decision_lpdf(t_adj | mu, boundary, bias,",
      " sigmadrift, sigmabias, sigmandt)",
      "\n      : cogmod_ddm_decision_lpdf(t_adj | -mu, boundary, 1 - bias,",
      " sigmadrift, sigmabias, sigmandt)"
    ),
    prelude = ".DDM_STAN_PRELUDE",
    ldens = function(t, k, p) .ddm_ldens(t, k, p),
    rng = function(n, p) .ddm_rng(n, p),
    # E[RT] has a closed form for the classic 4-parameter DDM, so unlike the
    # race families this one can offer posterior_epred() - see
    # posterior_epred_cogmod_ddm(), which owns the formula rather than the
    # registry, because it is a mean over a *choice* as well as a time.
    init = list(mu = 0, boundary = 1, bias = 0.5, sigmadrift = 0.3,
                sigmabias = 0.1, sigmandt = 0.05),
    # All three variability parameters are flat directions of the same kind as
    # cogmod_lba1()'s `sigmabias`: each has a floor at zero that its link only
    # reaches at minus infinity, and the likelihood stops changing well before
    # then. brms leaves all three flat, which is an improper posterior. They are
    # also the parameters the DDM literature agrees are hardest to recover, so
    # the priors are deliberately tight rather than merely proper.
    #
    # The natural-scale rows matter more here than elsewhere: leaving a
    # variability parameter out of bf() is the common way to write the model,
    # and brms then declares it as a plain auxiliary parameter. Fixing it
    # outright (`sigmadrift = 0` in bf()) removes it instead, and needs no prior.
    prior = list(
      sigmadrift = c(link = "normal(0, 1)", nat = "lognormal(-1, 0.75)",
                     slope = "normal(0, 0.5)"),
      sigmabias = c(link = "normal(-2, 1)", nat = "beta(1, 5)",
                    slope = "normal(0, 0.5)"),
      sigmandt = c(link = "normal(-3, 1)", nat = "lognormal(-3, 1)",
                   slope = "normal(0, 0.5)")
    ),
    dpar_doc = c(
      "dec: the observed choice. 1 is the UPPER boundary, 0 the lower.",
      "mu: drift rate, positive towards the boundary coded 1.",
      "boundary: boundary separation (> 0).",
      "bias: starting point as a proportion of `boundary`, in (0, 1).",
      "sigmadrift: between-trial SD of the drift rate (sv, >= 0).",
      paste("sigmabias: between-trial start-point range as a fraction in",
            "[0, 1) of the widest"),
      "   range that keeps the start point inside the boundaries.",
      paste("sigmandt: between-trial range of the non-decision time (st0),",
            "in the same unit as Y,"),
      "   with `ndt` its lower bound."
    ),
    note = c(
      "sigmandt is st0 itself, in seconds, not a fraction of anything. It used",
      "to be a multiple of a `minrt` dpar that no longer exists, and naming it",
      "after `tau` outlived that parameter too."
    ),
    label = "Drift Diffusion"
  )
)


# The choice families built on this machinery. Folded into .OUTLIER_FAMILIES in
# core_shifted.R, so with_outliers(), p_outlier(), cogmod_priors() and
# cogmod_inits() all reach them.
#' @keywords internal
.CHOICE_FAMILIES <- names(.CHOICE)


#' @keywords internal
.choice_spec <- function(name) {
  spec <- .CHOICE[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the shifted choice families. Supported: ",
      paste0(.CHOICE_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Either registry, for the code that is genuinely family-agnostic: p_outlier()
# and cogmod_inits(). The univariate machinery keeps using .shifted_spec(), so
# it cannot be handed a choice family by accident - the entries are not
# interchangeable, since a choice family's `ldens` needs the response index.
#' @keywords internal
.mixture_spec <- function(name) {
  spec <- .SHIFTED[[name]]
  if (is.null(spec)) spec <- .CHOICE[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the families built on the outlier mixture. ",
      "Supported: ", paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Family constructor ------------------------------------------------------

# Builds the brms custom family for `name`, appending `ndt` and `poutlier` to
# whatever decision dpars the registry lists. Same as .shifted_family() plus the
# `vars` the response index arrives through.
#' @keywords internal
.choice_family <- function(name, links = NULL, predict_outliers = FALSE) {
  spec <- .choice_spec(name)
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
    type = "real",
    vars = spec$vars
  )
  # It rides on the family because that is the only thing brms carries down to
  # a custom family's prediction methods. See ?rcogmod_lognormal.
  fam$predict_outliers <- isTRUE(predict_outliers)
  fam
}


# Stan code ---------------------------------------------------------------

# Generates the `<name>_lpdf` Stan function: the same mixture skeleton for every
# choice family, with only the decision density swapped in. The outlier
# component's scale is a literal because Stan functions cannot see the data
# block, and because a dpar would be estimated whenever the user left it out of
# the formula.
#
# brms appends `vars` after the dpars, so the response index is the *last*
# argument and is an `int`.
#' @keywords internal
.choice_lpdf <- function(name, prelude = "") {
  spec <- .choice_spec(name)
  if (!is.null(spec$prelude)) prelude <- get(spec$prelude)
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(.POUTLIER_SCALE, format = "g", digits = 15, width = 1)

  # As in .shifted_lpdf(): the outlier term has no parameter in it, so its
  # normalising constant is folded into a literal rather than left for Stan to
  # recompute on every leapfrog step. The `-log(K)` that makes the guess uniform
  # over the response options goes into the same constant - at K = 2 it cancels
  # the log(2) that folds the Normal onto [0, Inf) exactly.
  lc <- log(2) - 0.5 * log(2 * pi * .POUTLIER_SCALE^2) - log(spec$K)
  lp_out <- sprintf(
    "%s - %s * square(Y)",
    formatC(lc, format = "g", digits = 17, width = 1),
    formatC(1 / (2 * .POUTLIER_SCALE^2), format = "g", digits = 15, width = 1)
  )
  args <- paste(
    c(sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")), "int dec"),
    collapse = ", "
  )
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
// The outlier component is a guess: the choice is uniform over the %s response
// options and the RT is a half Normal with scale %s s. It keeps the density
// strictly positive below `ndt`, where the shifted decision component has none.
// That is what removes the hard min-RT boundary and lets `ndt` be estimated
// directly rather than as a fraction of an observed minimum. The scale is a
// constant in SECONDS. This family expects reaction times in seconds; give it
// another unit and the component contributes nothing anywhere in the data,
// which silently reinstates the min-RT boundary it exists to remove.
//
// The 1 / %s is what keeps the total density summing to one over the response
// options; without it it comes to 1 + poutlier. The leading constant below
// carries it, along with the normalising constant of the Student-t and the
// log(2) that folds it onto [0, Inf) - all of them constant, so writing them
// out saves Stan recomputing them for every observation on every leapfrog step.
%sreal %s_lpdf(%s) {
    // Parameter checks
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (dec < 0 || dec > %s) return negative_infinity();
    if (Y <= 0) return negative_infinity();

    real lp_out = %s;
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    real lp_dec = %s;
    return log_mix(poutlier, lp_out, lp_dec);
}
",
    prelude, spec$label, dpar_doc, spec$K, scale, spec$K, note, name, args,
    spec$stan_check, spec$K - 1L, lp_out, spec$stan_dens
  )
}


# Shared expose helper: compiles the generated lpdf and hands it back as an R
# function, for checking the Stan code against the R density.
#' @keywords internal
.choice_expose <- function(name) {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {
",
    .choice_lpdf(name),
    "}"
  )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions[[paste0(name, "_lpdf")]]
}


# Density, RNG ------------------------------------------------------------

# Recycle the decision dpars plus ndt/poutlier and the response to a common
# length. `...` holds the decision dpars, by name.
#' @keywords internal
.prepare_choice <- function(name, x = NULL, n = NULL, response = NULL, ndt,
                            poutlier, ...) {
  spec <- .choice_spec(name)
  dec <- list(...)
  missing <- setdiff(spec$dpars, names(dec))
  if (length(missing)) {
    stop("Missing parameter(s) for ", name, "(): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  dec <- dec[spec$dpars]

  # Checked against the registry's own bounds, so the R functions reject exactly
  # what `stan_check` rejects and the two cannot disagree. Lower bounds are open
  # unless the family says otherwise, which is why the default comparison is
  # `<=`; cogmod_rdm() closes both of its drift bounds, a zero drift being a
  # legitimate (if slow) accumulator there.
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
  if (any(ndt < 0, na.rm = TRUE)) stop("`ndt` must be non-negative.")
  if (any(poutlier < 0 | poutlier > 1, na.rm = TRUE)) {
    stop("`poutlier` must be in [0, 1].")
  }
  # The parentheses are load-bearing: `%in%` binds tighter than `-`, so
  # `response %in% seq_len(K) - 1L` compares against 1:K and then subtracts
  # from the logical.
  options <- seq_len(spec$K) - 1L
  if (!is.null(response)) {
    if (any(!(response %in% options), na.rm = TRUE)) {
      stop("`response` must be one of ", paste(options, collapse = ", "), ".",
           call. = FALSE)
    }
  }

  lens <- c(vapply(dec, length, integer(1)), length(ndt), length(poutlier),
            length(response))
  if (!is.null(x)) {
    if (is.null(response)) {
      stop("`response` must be provided alongside `x`.", call. = FALSE)
    }
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
  if (!is.null(response)) params$response <- rep_len(response, m)
  params$ndraws <- m
  params
}


# Log-density of the decision component alone, at t (which may be <= 0, giving
# -Inf) for response `k`. Shape is preserved, so this works on the
# draws x observations matrices p_outlier() builds.
#' @keywords internal
.ldec_choice <- function(name, t, k, p) {
  spec <- .choice_spec(name)
  # As for `.ldec()`: substitute a value the density will accept wherever the
  # time or a parameter is not one, and overwrite those results with -Inf below.
  # See `.dens_mask()` for why handing them the missing value instead throws.
  msk <- .dens_mask(spec, t, p)
  ok <- msk$ok
  # The response selects which accumulator won, so a missing one is no more
  # usable than a missing drift rate. It is not a dpar, so `.dens_mask()` does
  # not see it.
  kf <- !is.finite(k)
  if (any(kf)) {
    ok[] <- ok & rep_len(!kf, length(ok))
    k[kf] <- 0
  }
  ld <- spec$ldens(msk$t, k, msk$p)
  ld[!ok] <- -Inf
  ld[is.na(ld)] <- -Inf
  dim(ld) <- dim(t)
  ld
}


# Log-density of the outlier component: the half-t, thinned by the K options the
# guess is spread over.
#' @keywords internal
.lout_choice <- function(x, K) {
  .dcontam(x, log = TRUE) - log(K)
}


# Mixture density, shared by every choice family's d*() function.
#' @keywords internal
.dchoice <- function(name, x, response, ndt, poutlier,
                     log = FALSE, ...) {
  spec <- .choice_spec(name)
  params <- tryCatch(
    .prepare_choice(name, x = x, response = response, ndt = ndt,
                    poutlier = poutlier, ...),
    error = function(e) {
      warning(conditionMessage(e), ". Returning 0 density / -Inf log-density.")
      list(ndraws = length(x), error = TRUE)
    }
  )
  if (!is.null(params$error)) {
    return(rep(ifelse(log, -Inf, 0), params$ndraws))
  }

  lp_out <- .lout_choice(params$x, spec$K)
  lp_dec <- .ldec_choice(name, params$x - params$ndt, params$response, params)

  ld <- .log_mix(params$poutlier, lp_out, lp_dec)
  ld[is.na(ld)] <- -Inf
  if (log) ld else exp(ld)
}


# Mixture RNG, shared by every choice family's r*() function. Returns a data
# frame with `rt` and `response`, so an outlier draw carries a choice as well as
# a time - drawn uniformly over the K options, matching the density.
#
# `as_matrix = TRUE` returns the same two columns as a plain matrix instead. That
# is for `.posterior_predict_choice()`, which brms calls once per observation
# and which needs a matrix anyway: building a data frame here only to convert it
# back costs more than everything else in the call put together (a data frame is
# ~17x the price of the equivalent `cbind()`, and the round trip was two of them
# plus an `as.matrix()`, or about 63% of the total). The registry's `rng`
# entries return a bare list for the same reason.
#' @keywords internal
.rchoice <- function(name, n, ndt, poutlier, ..., as_matrix = FALSE) {
  spec <- .choice_spec(name)
  params <- .prepare_choice(name, n = n, ndt = ndt, poutlier = poutlier, ...)
  m <- params$ndraws

  out <- function(rt, response) {
    if (as_matrix) cbind(rt = rt, response = response)
    else data.frame(rt = rt, response = response)
  }

  rt <- numeric(m)
  response <- numeric(m)
  if (m == 0) return(out(rt, response))

  is_out <- stats::runif(m) < params$poutlier

  if (any(is_out)) {
    n_out <- sum(is_out)
    rt[is_out] <- abs(stats::rnorm(n_out, 0, .POUTLIER_SCALE))
    response[is_out] <- sample(seq_len(spec$K) - 1L, n_out, replace = TRUE)
  }
  if (any(!is_out)) {
    keep <- !is_out
    p <- lapply(params[spec$dpars], function(v) v[keep])
    sim <- spec$rng(sum(keep), p)
    rt[keep] <- sim$rt + params$ndt[keep]
    response[keep] <- sim$response
  }
  out(rt, response)
}


# brms methods ------------------------------------------------------------

# Pull the decision dpars for observation i (or all of them) off a brmsprep.
#' @keywords internal
.dpars_from_prep_choice <- function(name, prep, i = NULL) {
  spec <- .choice_spec(name)
  out <- lapply(spec$dpars, function(d) {
    if (is.null(i)) brms::get_dpar(prep, d) else brms::get_dpar(prep, d, i = i)
  })
  stats::setNames(out, spec$dpars)
}


# The response index brms carried through `vars`, for observation i.
#' @keywords internal
.dec_from_prep <- function(prep, i = NULL) {
  if (!"dec" %in% names(prep$data)) {
    stop("Decision variable 'dec' not found in prep$data. Add it to the ",
         "formula with `RT | dec(response) ~ ...`.", call. = FALSE)
  }
  if (is.null(i)) prep$data$dec else prep$data$dec[i]
}


#' @keywords internal
.log_lik_choice <- function(name, i, prep) {
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  dec <- .dpars_from_prep_choice(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)
  response <- .dec_from_prep(prep, i)

  n_draws <- max(vapply(c(dec, list(ndt, poutlier)), length, integer(1)))
  if (n_draws == 0) return(numeric(0))

  ll <- do.call(.dchoice, c(
    list(name = name, x = rep(y, length.out = n_draws),
         response = rep(response, length.out = n_draws), ndt = ndt,
         poutlier = poutlier, log = TRUE),
    dec
  ))
  ll[is.na(ll)] <- -Inf
  ll
}


# Returns a draws x 2 matrix: the RT and the simulated choice. That is the shape
# brms expects from a custom family whose response has more than one column.
#' @keywords internal
.posterior_predict_choice <- function(name, i, prep, predict_outliers = NULL) {
  dec <- .dpars_from_prep_choice(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }
  n_draws <- max(vapply(c(dec, list(ndt)), length, integer(1)))

  do.call(.rchoice, c(
    list(name = name, n = n_draws, ndt = ndt, poutlier = poutlier,
         as_matrix = TRUE),
    dec
  ))
}
