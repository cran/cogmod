#' Priors that make a cogmod posterior proper
#'
#' @description
#' Fills in weakly informative priors for the parameters `brms` would otherwise
#' leave flat, for the model you are actually fitting.
#'
#' **For [cogmod_lognormal()] this is not a convenience.** `brms` assigns a flat,
#' improper prior to the intercept of any custom-family parameter it does not
#' recognise, which there means both `ndt` and `poutlier`. The likelihood has
#' two directions in which it is exactly flat: `poutlier` toward 1, where every
#' response is attributed to the outlier component and `mu`, `sigma` and `ndt`
#' drop out of the density altogether; and `ndt` toward 0, where the model
#' reduces to an unshifted LogNormal and the gradient with respect to
#' `log(ndt)` vanishes. Flat prior plus infinite flat region is an improper
#' posterior. The fit does not fail loudly - it returns intercepts around `1e14`
#' with `Rhat` near 2 and an effective sample size of about 5.
#'
#' The second direction has nothing to do with the mixture; it is inherent to
#' putting a positive shift on a log link, which is why a prior on `poutlier`
#' alone is not enough.
#'
#' @details
#' The function starts from `brms::get_prior()` for the model in hand, edits the
#' rows it knows how to improve, and returns the result of
#' `brms::validate_prior()`. Three things follow.
#'
#' Every row comes from the model `brms` is going to build, so a prior matching
#' no parameter is impossible by construction: `0 + Intercept` formulas,
#' interactions, group-level terms and smooths are all handled because none of
#' them are guessed at.
#'
#' The return value is a `brmsprior` object: the same class returned by
#' [brms::get_prior()] and [brms::prior()]. It prints as a table, including the
#' defaults that `brms` will use, but can be combined with `c()` like any other
#' `brms` prior specification.
#'
#' Passing it through `brms::validate_prior()` means a malformed specification
#' errors here, with the offending row in view, rather than deep inside `brm()`.
#'
#' # Setting your own priors
#'
#' Combine it with [brms::prior()] entries. Set `replace = TRUE` to replace one
#' of these defaults, or omit it to add a prior for a different slot:
#'
#' ```r
#' priors <- c(
#'   cogmod_priors(f, df),
#'   brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
#'   replace = TRUE
#' )
#' ```
#'
#' # Centring the priors on a previous fit
#'
#' `warmstart` takes what [cogmod_warmstart()] takes - a `brmsfit`, the table
#' [as.data.frame()] makes of one, or the path to a CSV of it - and re-centres
#' the priors on that fit's posterior: each prior becomes
#' `normal(median, prior_scale * sd)`, the median and SD being the source's
#' for that same parameter. It applies to the population-level intercepts and
#' coefficients, the group-level SDs, and any dpar left out of the formula and
#' so declared as a plain auxiliary parameter. The group-level correlations
#' keep their LKJ, and the standardized effects have no stated prior to
#' change.
#'
#' No transformation is involved, which is the reason this is safe to do
#' automatically: the parameter a prior row is about *is* the parameter the
#' source sampled, matched by its `class`, `dpar`, `coef` and `group`. In
#' particular the `Intercept` prior is stated on the centred intercept in both
#' models, so that is what the median comes from.
#'
#' **This changes the posterior.** It is not in the same category as
#' [cogmod_inv_metric()] and [cogmod_step_size()], which only change how the
#' sampler moves: a prior is part of the model, and a fit with these priors is
#' answering a different question from one with the defaults. Two consequences
#' are worth stating plainly.
#'
#' First, **if the source was fitted to data the new model also contains, this
#' double-counts it.** A pilot on 5 of 10 participants, used to centre the
#' priors for the fit on all 10, uses those 5 participants twice - once as a
#' prior and once as data - and the intervals it produces are too narrow by an
#' amount nothing in the output reveals. That is the standard warm-start
#' setup, and it is exactly the setup in which this argument should not be
#' used for anything you intend to report. It is legitimate when the source is
#' an independent data set - last year's sample, another lab's, a different
#' session of the same task - or when you are deliberately doing a sequential
#' analysis and the priors *are* the previous posterior.
#'
#' Second, `prior_scale` is what stands between the two extremes.
#' `prior_scale = 1` uses the source's posterior as the prior, which is the
#' maximally informative choice and the one that double-counts hardest;
#' letting it grow widens the prior until it is only saying roughly where the
#' parameter lives. The default of 3 gives a prior with nine times the
#' variance of the source's posterior, so it carries something like a ninth of
#' the information - enough to put the sampler in the right region and keep it
#' out of the flat directions these likelihoods have, weak enough that data
#' disagreeing with the pilot will win.
#'
#' A parameter the source never saw - a new predictor, a group-level term it
#' did not have - keeps whatever prior the rest of this function gave it, and
#' [cogmod_warmstart()]'s `print()` says how many of those there are.
#'
#' ```r
#' # last year's sample as the prior for this year's
#' priors <- cogmod_priors(f, df, warmstart = fit_2025)
#' priors <- cogmod_priors(f, df, warmstart = fit_2025, prior_scale = 6)  # weaker
#' ```
#'
#' # Supported families
#'
#' The family is read off `formula`, so build it with
#' `brms::bf(..., family = cogmod_lognormal())`. Every family built on the direct
#' `ndt` + `poutlier` parameterization is edited: [cogmod_lognormal()],
#' [cogmod_logstudent()], [cogmod_loggamma()], [cogmod_invgaussian()], [cogmod_exwald()],
#' [cogmod_bisa()], [cogmod_gamma()],
#' [cogmod_invgamma()],
#' [cogmod_weibull()], [cogmod_invweibull()], [cogmod_logweibull()] and, for the
#' choice-and-RT models, [cogmod_lnr()], [cogmod_rdm()], [cogmod_lba2()] and
#' [cogmod_ddm()]. [cogmod_exgaussian()] and [cogmod_geg()] are edited too,
#' although they are not built on that parameterization - see their own section
#' below - as are the three bounded-scale families for subjective ratings,
#' [cogmod_choco()], [cogmod_betagate()] and [cogmod_betadiscrete()]. Any other
#' family, or a formula carrying none, is passed through: you get a message and
#' `brms`'s own defaults, unchanged, so the call is always safe to leave in a
#' script.
#'
#' What gets set, on the link scale (`log` for `ndt`, `logit` for `poutlier`,
#' `identity` for `shape`):
#'
#' | class | `ndt` | `poutlier` | `shape` |
#' | --- | --- | --- | --- |
#' | `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-1.2, 0.5)` | `normal(-5, 1)` | `normal(0, 0.5)` |
#' | `b` (slopes) | `normal(0, 0.2)` | `normal(0, 0.2)` | `normal(0, 0.2)` |
#' | `sd`, `sds` | `exponential(1)` | `exponential(1)` | `exponential(1)` |
#'
#' `normal(-1.2, 0.5)` centres `ndt` on 0.30 s, with 95% of its mass between
#' about 0.11 and 0.80 s: wide enough for the non-decision times of slower
#' populations and more demanding responses, and still a fence against the
#' `ndt -> 0` direction the likelihood cannot close on its own. Like everything
#' else in these families it is stated in **seconds**, which is the unit the
#' package expects throughout. `poutlier` is a proportion and
#' does not move: `normal(-5, 1)` is centred at about 0.7% and puts roughly 95%
#' of its mass between 0.1% and 5%. That is where the empirical estimates sit -
#' pooling across four lexical-decision megastudies, Miller (2024) puts the
#' outlier proportion below 0.5% and argues that the 5-10% assumed in most
#' simulation work is unrealistically large.
#'
#' `shape` exists only for [cogmod_loggamma()]. `normal(0, 0.5)` is centred on the
#' LogNormal shape and keeps the sampler clear of `sigma * shape >= 1`, where the
#' decision density becomes unbounded at the shift and the likelihood with it.
#'
#' A family may add rows of its own where its likelihood has a flat direction
#' that `brms` would leave improper. Seven do:
#'
#' - [cogmod_logstudent()]: `dof`. The LogNormal is only reached as `dof -> Inf`,
#'   and it is approached slowly - at `dof = 100` the density still differs from
#'   the LogNormal by 66% in the tail - so above about `dof = 30` the likelihood
#'   has almost stopped moving, and a `log` link puts that region at plus
#'   infinity. The prior fences the other end too: a Student-t is symmetric on
#'   the log scale, so a small `dof` piles mass just above `ndt` as well as in
#'   the slow tail, which is what `poutlier` is for. `normal(1.8, 0.7)` centres
#'   `dof` at 6 with 95% between 1.5 and 24.
#' - [cogmod_exwald()]: `tau`, the mean of the exponential residual stage. It is
#'   a length of time in seconds behind a `softplus` link, which `brms` has no
#'   way to know, and it shares a ridge with `ndt` - both delay the response,
#'   and only the shape of the leading edge tells them apart.
#'   `normal(-1.5, 0.7)` is the same statement [cogmod_exgaussian()] gets for the
#'   same quantity: roughly 55 to 630 ms.
#' - [cogmod_lba1()]: `sigmabias` and `boundary`. As the start-point range
#'   approaches zero the LBA converges to the recinormal and the likelihood
#'   stops depending on it, and a `softplus` link reaches zero only at minus
#'   infinity. Without those rows `sigmabias` runs off - `softplus(-10.4)`,
#'   `Rhat` 1.69.
#' - [cogmod_rdm()]: `sigmabias` and `boundary`, for exactly the reason
#'   [cogmod_lba1()] needs them - the two enter the threshold only through the
#'   sum `b = boundary + sigmabias` and trade off along a ridge worth a handful
#'   of log units, and the `softplus` link reaches zero only at minus infinity.
#'   The drift rates are left alone: a zero-drift accumulator still finishes,
#'   and still wins sometimes, so there is no plateau at the bottom of that
#'   direction the way there is for [cogmod_lnr()]'s `nuone`.
#' - [cogmod_lba2()]: `driftone`, `sigmazero`, `sigmaone`, `sigmabias` and
#'   `boundary`. The last two share the threshold ridge above; the two drift SDs
#'   are worse. The evidence scale of an LBA is arbitrary - multiply the drifts,
#'   their SDs, the start-point range and the threshold by any `c > 0` and every
#'   finishing time is unchanged - so the likelihood is *exactly* constant along
#'   that ray. These priors make the posterior proper; only fixing one SD in the
#'   formula (`sigmazero = 1`) identifies the scale. `driftone` has the plateau
#'   [cogmod_lnr()]'s `nuone` has, created here by the truncation of each drift
#'   at zero: once an accumulator rarely wins, its drift and SD are identified
#'   only through `|drift| / sigma^2`, and a flat prior lets the drift run off
#'   to `-12` and beyond. `normal(1, 2)` is centred where LBA drifts sit on the
#'   `sigmazero = 1` scale. See [cogmod_lba2()].
#' - [cogmod_ddm()]: `sigmadrift`, `sigmabias` and `sigmandt`, the three
#'   between-trial variability parameters. Each has a floor at zero that its link
#'   reaches only at minus infinity, and the likelihood stops changing well
#'   before then. They are also the parameters a DDM is least able to recover, so
#'   these priors are deliberately tighter than the rest. Fixing the ones a
#'   design cannot identify (`sigmadrift = 0` in `bf()`) is usually better than
#'   estimating them behind a prior.
#' - [cogmod_invgaussian()]: `sigmadrift` and `sigmandt`, for the same reasons
#'   as the DDM's, with `sigmandt` taking the DDM's prior for the same quantity
#'   verbatim. Both should be fixed at zero unless the design can identify
#'   them; `sigmandt` in particular needs a lot of data, a strong prior, or
#'   both.
#' - [cogmod_lnr()]: `nuone`, `sigmazero` and `sigmaone`. Push an accumulator's
#'   rate far enough down and it never finishes first, so the density depends on
#'   it only through a survival term that has already saturated at 1; past about
#'   `nuone = -6` the log-likelihood is exactly constant, and that accumulator's
#'   `sigma` is unidentified along with it. `mu` - which is `nuzero` - has the
#'   mirror-image plateau, but it is the response's own intercept and `brms`
#'   already gives it a proper `student_t` default, so it is left alone. Model a
#'   rarely-chosen option and it is worth mirroring the `nuone` prior onto it.
#'   `sigmabias`, the start-point range, has the flat direction
#'   [cogmod_lba1()]'s has - the likelihood stops changing as it approaches
#'   zero - and gets the same treatment, rescaled to its unit, the threshold
#'   offset of 1. Fix it at zero unless the design speaks to start-point
#'   variability.
#' - [cogmod_lognormal()]: `sigmabias`, the same start-point range as the
#'   LNR's (the LNR races two of these accumulators), with the same rows.
#'
#' All of them are the same failure as `ndt` and `poutlier`: an infinite flat
#' region under a flat prior. See [cogmod_lba1()], [cogmod_rdm()],
#' [cogmod_lba2()], [cogmod_ddm()] and [cogmod_lnr()].
#'
#' Note that no prior is set on the shape of [cogmod_weibull()] or
#' [cogmod_gamma()], although a shape below 2 makes their `ndt` gradient
#' unbounded. That region is reached because the *likelihood* prefers it, by
#' around 100 log units on the data in the [RT models article](https://dominiquemakowski.github.io/cogmod/articles/rt_models.html), so a prior weak
#' enough to be a sensible default cannot move the posterior out of it - only
#' bias it. `?rcogmod_weibull` sets out what to do instead.
#'
#' # The ex-Gaussian
#'
#' [cogmod_exgaussian()] has neither `ndt` nor `poutlier`, but all three of its
#' parameters are lengths of time in **seconds**, and `brms` has no way to know
#' that. `sigma` and `tau` sit behind a `softplus` link; `mu` is on `identity`.
#' All three intercepts are set.
#'
#' | class | `sigma` | `tau` |
#' | --- | --- | --- |
#' | `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-2.3, 0.7)` | `normal(-1.5, 0.7)` |
#' | `b` (slopes) | `normal(0, 0.5)` | `normal(0, 0.5)` |
#' | `sd`, `sds` | `exponential(1)` | `exponential(1)` |
#'
#' On the `softplus` scale `normal(-2.3, 0.7)` puts `sigma` - the SD of the
#' Gaussian component - between roughly 25 and 330 ms with a median of 96 ms,
#' and `normal(-1.5, 0.7)` puts `tau` - the mean of the exponential tail -
#' between roughly 55 and 630 ms with a median of 201 ms. Both cover the range
#' these parameters occupy across the usual simple- and choice-RT tasks and are
#' wide enough not to fight data that disagrees.
#'
#' `tau` arrives from `brms` flat, so it is filled like any other unrecognised
#' custom dpar. `sigma` arrives with a **non-empty** default,
#' `student_t(3, 0, 2.5)`, because `brms` recognises the name from its own
#' families - centred on `softplus(0) = 0.69 s`, a Gaussian SD wider than most
#' whole RT distributions. That one is overridden rather than filled, the same
#' treatment `shape` and an omitted `ndt` get above.
#'
#' `mu` gets `normal(0.4, 0.25)` on its own intercept - 95% of the mass between
#' -0.09 and 0.89 s. It is not left to `brms`, whose `student_t(3, 0, 2.5)`
#' is a fair statement about a location on a `softplus` link (median 0.69 s) but
#' not on the `identity` link `mu` uses, where it is centred on zero seconds
#' and puts a Gaussian centre of -2 s on a par with one of +2 s. The prior does
#' **not** exclude negative values: `mu` is a location, and for fast
#' heavily-tailed data the Gaussian component genuinely belongs near or below
#' zero with `tau` carrying the mass. Only the intercept is set - the response's
#' slopes are the effects being estimated, and are left to `brms`. Note that `mu`
#' is the centre of the Gaussian component **alone**, so the mean of the
#' distribution it implies is `mu + tau`; see [cogmod_exgaussian()].
#'
#' # The bounded-scale families
#'
#' [cogmod_choco()], [cogmod_betagate()] and [cogmod_betadiscrete()] model
#' subjective ratings on `[0, 1]` or on `1:k` rather than reaction times, and
#' none of them has an infinite flat region of the [cogmod_lognormal()] kind.
#' They are edited anyway, because every distributional parameter they have
#' arrives either flat - improper - or with a `brms` default aimed at a
#' different parameterization.
#'
#' | class | `pmid`, `pzero` | `pex` | `bex`, `confright`, `confleft` | `precright`, `precleft`, `phi` (softplus) | `phi` ([cogmod_betadiscrete()], `log`) |
#' | --- | --- | --- | --- | --- | --- |
#' | `Intercept` | `normal(-2.5, 1)` | `normal(-2, 1)` | `normal(0, 1)` | `normal(2, 1.5)` | `normal(0.7, 0.8)` |
#' | `b` (slopes) | `normal(0, 0.5)` | `normal(0, 0.5)` | `normal(0, 0.5)` | `normal(0, 0.5)` | `normal(0, 0.5)` |
#' | `sd`, `sds` | `exponential(1)` | `exponential(1)` | `exponential(1)` | `exponential(1)` | `exponential(1)` |
#'
#' **The point masses.** `pmid` is the probability of landing exactly on the
#' midpoint of the scale and `pzero` the probability of an extra category
#' outside it. Both are flat on a `logit` link, which is improper in the
#' posterior as well as the prior whenever the event is simply absent from the
#' data: with no exact midpoints anywhere the likelihood in `pmid` increases
#' monotonically all the way to zero, and nothing stops the logit running to
#' minus infinity. That is `poutlier`'s failure exactly, so it gets
#' `poutlier`'s treatment, including the omitted form's **mode at zero** -
#' leaving the dpar out of `bf()` is itself the statement that you do not
#' expect the event. `normal(-2.5, 1)` is centred on 8% with 95% of its mass
#' below 37%, which covers both a slider with a visible midpoint tick and a
#' scale carrying a "not applicable" category.
#'
#' This is not hypothetical, and it has the signature the [cogmod_lognormal()]
#' failure has. Fitting [cogmod_choco()] to 400 slider responses containing no
#' exact midpoints and no extremes, `brms`'s flat defaults put `pex` at
#' `-1.1e14` and `pmid` at `-3.4e13`, both with `Rhat` 2.9 and an effective
#' sample size of 5, alongside 12% divergent transitions and 48% of them
#' hitting the maximum treedepth. With these priors the same data give -4.70
#' and -5.15 - 0.9% and 0.6% on the probability scale, which is what a data
#' set containing none of either should say - at `Rhat` 1.00, 3000 to 4000
#' effective samples, and no divergences.
#'
#' **The extremes.** `pex`, the total probability of a 0 or a 1, deliberately
#' does *not* get the mode-at-zero treatment: extreme responding is a
#' documented response style rather than a contaminant, and a scale nobody ever
#' answers at the endpoints is the unusual case. `normal(-2, 1)` centres it on
#' 12% with 95% between 2% and 49%. `bex`, which end the extremes favour, is
#' centred on symmetric; that also fences the two values at which the model
#' degenerates, since at `bex = 0` or `1` one gate closes entirely and any
#' observation at that endpoint takes the density to `-Inf`.
#'
#' **The Beta shapes.** The precisions are the awkward ones. The underlying
#' Beta has shapes `conf * prec * 2` and `(1 - conf) * prec * 2`, so a
#' precision below about 1 makes it U-shaped and unbounded at both ends, while
#' a large one narrows it sharply: at a precision of 15 the middle 95% of a
#' side spans only a third of it, and a rating a tenth of the way along that
#' side costs 13 log units. `normal(2, 1.5)` on the
#' `softplus` link puts the precision between 0.33 and 4.95 with a median of
#' 2.13, leaving 17% of its mass below 1 - a U-shaped rating distribution is
#' a real thing for a polarising item, so it is fenced rather than excluded.
#'
#' `phi` is the one row of the three that **overrides** a non-empty `brms`
#' default rather than filling an empty one, for the reason
#' [cogmod_exgaussian()]'s `sigma` does: `brms` recognises the name from its
#' own beta family and supplies `student_t(3, 0, 2.5)`. On
#' [cogmod_betagate()]'s `softplus` link that is merely loose. On
#' [cogmod_betadiscrete()]'s `log` link it is close to improper - its 95%
#' interval runs to `phi = 2853`, where the Beta has collapsed onto one rating
#' category and every other category's probability has underflowed.
#' [cogmod_choco()]'s precisions escape it only because they are named
#' `precright` and `precleft`, which `brms` does not recognise, so they arrive
#' flat instead.
#'
#' `mu` is left to `brms` in all three. It is the response's own predictor, and
#' on a `logit` link `student_t(3, 0, 2.5)` is the standard weakly informative
#' choice for exactly that - unlike [cogmod_exgaussian()]'s `mu`, which needed
#' overriding only because an `identity` link made the same prior a statement
#' about seconds.
#'
#' # Parameters left out of the formula
#'
#' Writing `ndt ~ 1` and omitting `ndt` entirely are not the same thing to
#' `brms`, and the difference matters here. A dpar that appears in `bf()` - even
#' as `~ 1` - gets a linear predictor, so it is estimated on the **link** scale
#' and reported under *Regression Coefficients* as `ndt_Intercept`. A dpar left
#' out is declared as a plain auxiliary parameter instead, the same mechanism as
#' `sigma` for `gaussian()`, estimated on the **natural** scale with no link and
#' reported under *Further Distributional Parameters*.
#'
#' Both forms are filled, with priors on whichever scale the parameter actually
#' lives on:
#'
#' | dpar | in `bf()` (link scale) | omitted (natural scale) |
#' | --- | --- | --- |
#' | `ndt` | `normal(-1.2, 0.5)` | `lognormal(-1.2, 0.5)` |
#' | `poutlier` | `normal(-5, 1)` | `exponential(100)` |
#' | `shape` | `normal(0, 0.5)` | `normal(0, 0.5)` |
#' | `sigmabias`, `boundary` ([cogmod_lba1()], [cogmod_rdm()], [cogmod_lba2()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#' | `sigmazero`, `sigmaone` ([cogmod_lba2()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#' | `driftone` ([cogmod_lba2()]) | `normal(1, 2)` | `normal(1, 2)` |
#' | `sigmadrift` ([cogmod_ddm()]) | `normal(0, 1)` | `lognormal(-1, 0.75)` |
#' | `sigmabias` ([cogmod_ddm()]) | `normal(-2, 1)` | `beta(1, 5)` |
#' | `sigmandt` ([cogmod_ddm()], [cogmod_invgaussian()]) | `normal(-3, 1)` | `lognormal(-3, 1)` |
#' | `sigmadrift` ([cogmod_invgaussian()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#' | `nuone` ([cogmod_lnr()]) | `normal(0.7, 1.5)` | `normal(0.7, 1.5)` |
#' | `sigmazero`, `sigmaone` ([cogmod_lnr()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#' | `sigmabias` ([cogmod_lnr()], [cogmod_lognormal()]) | `normal(0, 1)` | `lognormal(-0.35, 0.75)` |
#' | `dof` ([cogmod_logstudent()]) | `normal(1.8, 0.7)` | `lognormal(1.8, 0.7)` |
#' | `tau` ([cogmod_exwald()]) | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` |
#' | `mu` ([cogmod_exgaussian()]) | `normal(0.4, 0.25)` | - (always modelled) |
#' | `sigma` ([cogmod_exgaussian()]) | `normal(-2.3, 0.7)` | `lognormal(-2.3, 0.7)` |
#' | `tau` ([cogmod_exgaussian()]) | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` |
#' | `confright`, `confleft` ([cogmod_choco()]), `bex` | `normal(0, 1)` | `beta(2, 2)` |
#' | `precright`, `precleft` ([cogmod_choco()]), `phi` ([cogmod_betagate()]) | `normal(2, 1.5)` | `lognormal(0.7, 0.7)` |
#' | `phi` ([cogmod_betadiscrete()]) | `normal(0.7, 0.8)` | `lognormal(0.7, 0.8)` |
#' | `pex` ([cogmod_choco()], [cogmod_betagate()]) | `normal(-2, 1)` | `beta(2, 12)` |
#' | `pmid` ([cogmod_choco()]), `pzero` ([cogmod_betadiscrete()]) | `normal(-2.5, 1)` | `exponential(9)` |
#'
#' The `ndt` pair describes the same belief twice: `lognormal` is just `normal`
#' on the log scale, written for the untransformed parameter. The
#' [cogmod_exgaussian()] pairs do the same to within a rounding error, because
#' `softplus(x)` and `exp(x)` agree to three figures for the `x` below `-2` that
#' both parameters live at: `lognormal(-2.3, 0.7)` has median 0.100 against
#' `softplus(-2.3) = 0.096`.
#'
#' The Beta precisions are the one place that reasoning does **not** carry
#' over. They live up at `x = 2`, where `softplus(x)` is approximately `x`
#' rather than `exp(x)`, so the two forms cannot be the same distribution
#' written twice however the numbers are chosen. `lognormal(0.7, 0.7)` is
#' matched to `normal(2, 1.5)` on its median - 2.01 against 2.13 - and on its
#' sense rather than quantile by quantile: its upper tail reaches 7.9 where the
#' link form stops at 4.9. [cogmod_betadiscrete()]'s `phi` is the exception
#' that shows the rule, being on a `log` link, where the pair really is one
#' distribution and carries the same two numbers.
#'
#' If the data were trimmed before fitting, tighten this rather than removing
#' it: `normal(-7, 0.5)` asserts essentially no contamination while keeping the
#' density positive below `ndt`. Fixing `poutlier = 0` outright reinstates the
#' hard min-RT boundary that the component exists to remove. See the *Trimmed
#' data* section of [cogmod_lognormal()].
#'
#' `poutlier` is deliberately *not* the same belief twice. Leaving it out of the
#' formula is itself information - you either trimmed the data already or do not
#' expect outliers - so the omitted form puts its **mode at zero**, which a
#' logit-scale prior cannot do at any location. The centre is unchanged:
#' `exponential(100)` has median 0.0069 against `plogis(-5) = 0.0067`. It is
#' still a prior rather than a constraint, so a genuine spike of fast responses
#' will still pull the rate up; to switch the parameter off entirely, trim the
#' data and it will simply sit near zero.
#'
#' `pmid` and `pzero` are the same argument on a rating scale, and get the same
#' pair: `exponential(9)` has median 0.077 against `plogis(-2.5) = 0.076`. To
#' switch either off outright, fix it in `bf()` - `pmid = 0` or `pzero = 0` -
#' which removes the parameter rather than merely pushing it down. Note that
#' `pmid = 0` makes any response falling exactly on the midpoint impossible,
#' so the fit will fail to initialise if the data contain one.
#'
#' Two rows override a **non-empty** `brms` default rather than filling an empty
#' one. `brms` recognises the name `ndt` from its own shifted families and
#' supplies `uniform(0, min_Y)` - precisely the min-RT bound that
#' [cogmod_lognormal()]'s parameterization exists to remove, reimposed silently
#' and with a warning about an upper bound on an unbounded parameter. It
#' recognises `sigma` too, and gives [cogmod_exgaussian()]'s a half
#' `student_t(3, 0, 2.5)`, whose median of 1.9 s is a Gaussian SD wider than
#' most whole RT distributions. An omitted `poutlier` is left flat
#' over `[0, 1]` by `brms`, which is proper but puts half its mass above 0.5, and
#' an omitted `shape` or `tau` is flat over the whole real line, which is not
#' proper at all.
#'
#' Slope and group-level priors are deliberately narrow. On a log or a logit
#' link a flat slope prior is not as harmless as it looks, and a group-level SD
#' with no prior can wander far enough for individual groups to reach the flat
#' regions above even when the population intercept is well behaved.
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = cogmod_lognormal())`.
#' @param data The data, as passed to `brms::brm()`.
#' @param ... Passed to `brms::get_prior()` and `brms::validate_prior()`, for
#'   arguments such as `data2` or `knots`.
#' @param warmstart Optional. A previous fit to centre the priors on: a
#'   `brmsfit`, a `cogmod_warmstart` object, its data frame, or the path to a
#'   CSV of it, as [cogmod_warmstart()] takes. `NULL` (the default) leaves the
#'   priors where this function put them. **Changes the posterior, and
#'   double-counts the source's data if the new model contains it** - see the
#'   section above.
#' @param prior_scale The prior SD, as a multiple of the source's posterior
#'   SD. Only used with `warmstart`. Larger is weaker; 1 would use the
#'   source's posterior as the prior.
#'
#' @return A `brmsprior` object, to pass to `brms::brm(prior = )`.
#'
#' @references
#' - Miller, J. (2024). Estimating the proportions and latencies of reaction
#'     time outliers: A pooling method and case study of lexical decision
#'     tasks. *Behavior Research Methods*, *56*(7), 7280-7306.
#'     \doi{10.3758/s13428-024-02419-y}
#'
#' @seealso [cogmod_inits()], [cogmod_stanvars()], [cogmod_warmstart()]
#'
#' @examples
#' d <- data.frame(RT = rcogmod_lognormal(50, ndt = 0.3, poutlier = 0.02))
#' f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
#'   family = cogmod_lognormal()
#' )
#' cogmod_priors(f, d)
#'
#' # Replace a default, or append a prior for another parameter.
#' priors <- c(
#'   cogmod_priors(f, d),
#'   brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
#'   replace = TRUE
#' )
#'
#' @export
cogmod_priors <- function(formula, data, ..., warmstart = NULL, prior_scale = 3) {
  family <- .cogmod_family(formula)
  fam <- .family_name(family)

  # Before anything else: is the response the kind of thing this family can be
  # given at all? Milliseconds are the common case and they are silent - see the
  # header of cogmod_checkdata.R for why this runs here rather than in one of
  # the other two generics.
  .cogmod_checkdata(formula, data, family)

  if (isTRUE(fam %in% .OUTLIER_FAMILIES)) {
    out <- .priors_shifted(formula, data, family, ...)
  } else if (isTRUE(fam %in% names(.PRIORS_PLAIN))) {
    spec <- .PRIORS_PLAIN[[fam]]
    out <- .priors_dpars(formula, data, family, spec$prior, spec$override, ...)
  } else {
    message(
      "cogmod_priors() has nothing to add for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'; returning the brms defaults unchanged. ",
      "Currently supported: ", paste(.priors_families(), collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_lognormal())."
    )
    out <- brms::empty_prior()
  }

  # Last, so that it can override any of the above: a previous fit's posterior
  # is a stronger statement about where a parameter is than a default location
  # chosen for the family in general.
  if (!is.null(warmstart)) {
    out <- .priors_warmstart(out, formula, data, family, warmstart, prior_scale, ...)
  }

  args <- list(out, formula, data, ...)
  if (!is.null(family)) args$family <- family
  do.call(brms::validate_prior, args)
}


# Re-centres the priors on a previous fit's posterior: normal(median, scale *
# sd) for every prior row whose parameter the source has a posterior for.
#
# The join is on the address a get_prior() row carries - class, dpar, coef and
# group - which cogmod_warmstart() attaches to each unconstrained parameter as
# it labels the metric. That is what makes the scales line up: the Stan
# parameter a prior row is about is the same parameter the source sampled, so
# no transformation is involved. `Intercept` is the centred intercept in both,
# which is the one brms states the "Intercept" prior on.
#
# `base` is what the rest of cogmod_priors() decided; rows it set for a
# parameter the source knows are replaced, rows for anything else survive.
#' @keywords internal
.priors_warmstart <- function(base, formula, data, family, warmstart, prior_scale, ...) {
  if (!is.numeric(prior_scale) || length(prior_scale) != 1L ||
      !is.finite(prior_scale) || prior_scale <= 0) {
    stop("`prior_scale` must be a single positive number: the multiple of the ",
         "source's posterior SD to use as the prior SD.", call. = FALSE)
  }
  post <- cogmod_warmstart(warmstart, formula = formula, data = data,
                           jitter = 0, ...)$prior
  if (!nrow(post)) {
    message("The warm-start source carries no posterior median and SD, so there ",
            "is nothing to centre the priors on; leaving them as they are. A ",
            "table written before cogmod_priors(warmstart = ) existed has only ",
            "the metric and the means - rewrite it from the fit.")
    return(base)
  }

  args <- list(formula, data = data, ...)
  if (!is.null(family)) args$family <- family
  p <- do.call(brms::get_prior, args)
  # Kept whole: the only place the full coefficient list still exists, which
  # is what says whether a blanket row has been left covering nothing.
  all_rows <- p

  key <- function(d) paste(d$class, d$dpar, d$coef, d$group, sep = "|")
  idx <- match(key(p), key(post))
  hit <- !is.na(idx) & !nzchar(p$resp)
  if (!any(hit)) {
    message("None of the model's prior slots matched a parameter in the warm-start ",
            "source, so the priors are unchanged. Usually this means the source ",
            "was written for a different formula.")
    return(base)
  }

  new <- p[hit, , drop = FALSE]
  new$prior <- sprintf("normal(%s, %s)",
                       .num_str(post$median[idx[hit]]),
                       .num_str(prior_scale * post$sd[idx[hit]]))

  out <- rbind(base[!key(base) %in% key(new), , drop = FALSE], new)
  # Setting every coefficient a blanket row covers leaves the blanket covering
  # nothing, which brms warns about; the same check .priors_dpars() ends with.
  gone <- vapply(seq_len(nrow(out)), .blanket_now_unused, logical(1),
                 p = out, all_rows = all_rows)
  out[!gone, , drop = FALSE]
}


# A number as Stan will read it back. "g" keeps eight significant digits and
# falls to exponent notation only when that is shorter, which Stan accepts.
#' @keywords internal
.num_str <- function(x) trimws(formatC(x, digits = 8, format = "g"))


# Every family cogmod_priors() edits, for the message above.
#' @keywords internal
.priors_families <- function() unique(c(.OUTLIER_FAMILIES, names(.PRIORS_PLAIN)))


# Families that are NOT on the ndt + poutlier mixture but still have dpars brms
# leaves flat, or gives a default that is wrong for a parameter measured in
# seconds.
#
# Same shape as the `prior` slot in the .SHIFTED / .CHOICE registries: one named
# character vector per dpar, with any of `link` (the dpar is modelled in bf(), so
# it lives on its link scale), `nat` (the dpar is omitted from bf(), so brms
# declares it as a plain auxiliary parameter on the natural scale) and `slope`.
#
# `override` names the dpars whose *modelled intercept* row is replaced even
# though brms already supplied a prior for it. Auxiliary rows are always
# replaced - that is what the `nat` entry is for - so they need no listing.
#' @keywords internal
.PRIORS_PLAIN <- list(
  # cogmod_exgaussian() has no ndt and no poutlier, but it does have two
  # parameters that are lengths of time in seconds, behind a softplus link, and
  # brms has no way to know that.
  #
  # `tau` arrives flat: brms does not recognise the name, so a custom dpar gets
  # an improper prior over the whole real line, on a link whose upper reaches
  # are RTs of tens of seconds.
  #
  # `sigma` arrives with student_t(3, 0, 2.5), which brms supplies because it
  # recognises the *name* from its own families - and which is centred on
  # softplus(0) = 0.69 s, a Gaussian SD wider than most whole RT distributions,
  # with the bulk of its mass out past a second. That is why `sigma` is in
  # `override`: the row is not empty, so filling it is not enough.
  #
  # On the softplus scale normal(-2.3, 0.7) puts `sigma` at 25-333 ms with a
  # median of 96 ms, and normal(-1.5, 0.7) puts `tau` at 55-631 ms with a median
  # of 201 ms - the range these two occupy across the usual simple- and
  # choice-RT tasks, and wide enough not to fight data that disagrees.
  #
  # The `nat` pair says the same thing for the untransformed parameter. Because
  # softplus(x) ~ exp(x) for x well below zero - the region both of these live in
  # - a lognormal with the same numbers lands in almost the same place:
  # lognormal(-2.3, 0.7) has median 0.100 against softplus's 0.096.
  #
  # `mu` is not left to brms. Its student_t(3, 0, 2.5) is a reasonable statement
  # about a Gaussian centre on the SOFTPLUS scale, where it has a median of
  # 0.69 s, but `mu` is on an identity link, and the same prior there is centred
  # on zero seconds with a scale of 2.5 s - it says a Gaussian component centred
  # at -2 s is as plausible as one at +2 s. So it needs one of its own, and gets
  # normal(0.4, 0.25): 95% of its mass in -0.09 to 0.89 s.
  #
  # Note it deliberately does NOT exclude negative values. `mu` is a location,
  # and for fast heavily-tailed data the Gaussian component genuinely belongs
  # near or below zero with `tau` carrying the mass - the prior is there to say
  # where RT data usually sit, not to reimpose the constraint the identity link
  # was chosen to remove. Note also that it is the centre of the Gaussian
  # component alone, so the distribution it implies has mean `mu + tau`; see
  # ?rcogmod_exgaussian.
  #
  # Only the intercept is set. The response's slopes are left to brms, because
  # they are the effects the model is being fitted to estimate.
  cogmod_exgaussian = list(
    prior = list(
      mu = c(link = "normal(0.4, 0.25)"),
      sigma = c(link = "normal(-2.3, 0.7)", nat = "lognormal(-2.3, 0.7)",
                slope = "normal(0, 0.5)"),
      tau = c(link = "normal(-1.5, 0.7)", nat = "lognormal(-1.5, 0.7)",
              slope = "normal(0, 0.5)")
    ),
    override = "sigma"
  ),

  # cogmod_geg() is the ex-Gaussian with its CDF raised to `shape`, so mu, sigma
  # and tau get exactly the same treatment - and then `shape` needs a prior of
  # its own for a reason the other three do not have.
  #
  # `shape` is only weakly identified against `mu`: at the maximum likelihood the
  # two correlate about -0.98 on real RT data, because raising the CDF to a power
  # shifts the distribution as well as bending it. Left flat, `shape` will happily
  # wander a long way up that ridge and drag `mu`, `sigma` and `tau` with it,
  # which is how a family that nests the ex-Gaussian ends up reporting a `mu`
  # nowhere near the ex-Gaussian's.
  #
  # On the log link zero is `shape = 1`, the ex-Gaussian, so normal(0, 0.5) is
  # centred on the nested model and puts 95% of its mass in shape 0.37 to 2.7.
  # That is wide enough for the negative-skew and heavy-kurtosis shapes the
  # family exists to reach, and narrow enough that the ridge stays fenced.
  cogmod_geg = list(
    prior = list(
      mu = c(link = "normal(0.4, 0.25)"),
      sigma = c(link = "normal(-2.3, 0.7)", nat = "lognormal(-2.3, 0.7)",
                slope = "normal(0, 0.5)"),
      tau = c(link = "normal(-1.5, 0.7)", nat = "lognormal(-1.5, 0.7)",
              slope = "normal(0, 0.5)"),
      shape = c(link = "normal(0, 0.5)", nat = "lognormal(0, 0.5)",
                slope = "normal(0, 0.5)")
    ),
    override = c("sigma", "shape")
  ),

  # The three bounded-scale families for subjective ratings. None is on the
  # ndt + poutlier mixture and none has an infinite flat region of the
  # cogmod_lognormal() kind, but every dpar they have arrives either flat -
  # improper - or with a brms default aimed at a different parameterization.
  # Four arguments cover all three families, so they are set out once here
  # rather than three times below.
  #
  # THE POINT MASSES. `pmid` (cogmod_choco) and `pzero` (cogmod_betadiscrete)
  # are the probability of landing exactly on the midpoint of the scale, or on
  # an extra category outside it. brms leaves both flat on a logit link, and
  # that is improper in the posterior as well as the prior whenever the event
  # is simply absent from the data: with no exact midpoints anywhere, the
  # likelihood in `pmid` increases monotonically all the way to zero and
  # nothing stops the logit running to minus infinity. That is exactly
  # `poutlier`'s failure, so it gets exactly `poutlier`'s treatment - a normal
  # on the link scale, and on the natural scale a distribution with its MODE at
  # zero, which no logit-scale prior can have at any location. Omitting the
  # dpar from bf() is itself the statement that you do not expect the event.
  # The two forms agree on the centre, as that pair does: exponential(9) has
  # median 0.077 against plogis(-2.5) = 0.076, about 8%, with 95% of the mass
  # below 37%. Sliders with a visible midpoint tick and scales carrying a "not
  # applicable" category both sit in that range. exponential(9) puts 1.2e-4 of
  # its mass above 1, which the declared upper bound truncates away.
  #
  # THE EXTREMES. `pex` is the total probability of a 0 or a 1, and it is NOT
  # given the mode-at-zero treatment, because it is a different kind of
  # quantity: extreme responding is a documented response style rather than a
  # contaminant, and a scale on which nobody ever picks an endpoint is the
  # unusual case. normal(-2, 1) centres it on 12% with 95% between 2% and 49%,
  # which spans the ordinary and the strongly extreme responder. `bex`, which
  # end the extremes favour, is centred on symmetric with normal(0, 1) - 12% to
  # 88%. That also fences the two values at which the model degenerates: at
  # bex = 0 or 1 one of the two gates closes entirely, and any observation at
  # that endpoint takes the density to -Inf.
  #
  # THE BETA SHAPES. `confright` / `confleft` (cogmod_choco) and `mu`
  # (cogmod_betagate) are the Beta mean of each side; normal(0, 1) centres them
  # on the middle of their half of the scale, 12% to 88%. The precisions are
  # the awkward ones. Shapes are `conf * prec * 2` and `(1 - conf) * prec * 2`,
  # so prec below about 1 makes the Beta U-shaped and unbounded at both ends,
  # while a large prec narrows it sharply: at prec 15 the middle 95% of a side
  # spans only a third of it and a rating a tenth of the way along costs 13 log
  # units. normal(2, 1.5) on the softplus link puts prec between 0.33 and 4.95
  # with a median of 2.13, leaving 17% of the mass below 1 - the U shape is a
  # real rating distribution for a polarising item, so it is fenced, not
  # excluded.
  #
  # Note what the `nat` forms here do NOT do. For `ndt` the pair is one belief
  # written twice, because lognormal is normal on the log scale; for
  # cogmod_exgaussian() it is very nearly that, because softplus(x) and exp(x)
  # agree to three figures below x = -2. Neither holds for a precision, which
  # lives up at x = 2, where softplus(x) is approximately x rather than exp(x).
  # So lognormal(0.7, 0.7) is matched to normal(2, 1.5) on its median (2.01
  # against 2.13) and on its sense, not quantile by quantile - its upper tail
  # reaches 7.9 where the link form stops at 4.9. cogmod_betadiscrete()'s `phi`
  # is the exception that shows the rule: it is on a LOG link, so its two forms
  # are the same distribution and carry the same numbers.
  #
  # `phi` is the one row of the three that OVERRIDES rather than fills. brms
  # recognises the name from its own beta family and supplies
  # student_t(3, 0, 2.5), the same trap cogmod_exgaussian()'s `sigma` falls
  # into. On cogmod_betagate()'s softplus link that is merely loose. On
  # cogmod_betadiscrete()'s log link it is close to improper: its 95% interval
  # runs to phi = 2853, where the Beta has collapsed onto a single rating
  # category and every other category's probability has underflowed.
  # cogmod_choco()'s precisions escape only because they are named `precright`
  # and `precleft`, which brms does not recognise, so they arrive flat instead.
  #
  # `mu` is left to brms in all three. It is the response's own predictor, and
  # on a logit link brms' student_t(3, 0, 2.5) is the standard weakly
  # informative choice for exactly that - unlike cogmod_exgaussian()'s `mu`,
  # which needed overriding only because an identity link made the same prior a
  # statement about seconds. The dpar slopes get normal(0, 0.5) rather than the
  # 0.2 default: on a logit scale that is an odds ratio between 0.38 and 2.7
  # per unit of predictor, and these slopes are what the model is fitted to
  # estimate.
  cogmod_choco = list(
    prior = list(
      confright = c(link = "normal(0, 1)", nat = "beta(2, 2)",
                    slope = "normal(0, 0.5)"),
      confleft = c(link = "normal(0, 1)", nat = "beta(2, 2)",
                   slope = "normal(0, 0.5)"),
      precright = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
                    slope = "normal(0, 0.5)"),
      precleft = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
                   slope = "normal(0, 0.5)"),
      pex = c(link = "normal(-2, 1)", nat = "beta(2, 12)",
              slope = "normal(0, 0.5)"),
      bex = c(link = "normal(0, 1)", nat = "beta(2, 2)",
              slope = "normal(0, 0.5)"),
      pmid = c(link = "normal(-2.5, 1)", nat = "exponential(9)",
               slope = "normal(0, 0.5)")
    )
    # Nothing to override: brms recognises none of these names, so every row
    # arrives empty.
  ),

  cogmod_betagate = list(
    prior = list(
      phi = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
              slope = "normal(0, 0.5)"),
      pex = c(link = "normal(-2, 1)", nat = "beta(2, 12)",
              slope = "normal(0, 0.5)"),
      bex = c(link = "normal(0, 1)", nat = "beta(2, 2)",
              slope = "normal(0, 0.5)")
    ),
    override = "phi"
  ),

  cogmod_betadiscrete = list(
    prior = list(
      # On a log link, so the two forms are one distribution: phi between 0.42
      # and 9.7, median 2.0, 19% of the mass below 1. Below 1 with `mu` at 0.5
      # is the U and J shaped rating data this family exists to fit, so it
      # keeps real weight. The upper end is where the scale stops resolving:
      # at phi 10 a seven-point scale already has 95% of its mass inside three
      # of its seven categories, and it only narrows from there.
      phi = c(link = "normal(0.7, 0.8)", nat = "lognormal(0.7, 0.8)",
              slope = "normal(0, 0.5)"),
      pzero = c(link = "normal(-2.5, 1)", nat = "exponential(9)",
                slope = "normal(0, 0.5)")
    ),
    override = "phi"
  )
)


# The two parameters the ndt + poutlier parameterization introduces, plus the
# one decision dpar whose brms default is actively harmful, in the same form as
# a family's own `prior` slot. Merged with that slot in .priors_shifted().
#
# ndt is a location in time, in SECONDS: exp(-1.2) = 0.30 s. Nothing in these
# families is equivariant to the unit of measurement any more - see
# .POUTLIER_SCALE in core_shifted.R - so the location is a constant. The `nat`
# entry is the same distribution written for the untransformed parameter.
#
# The SD is 0.5 on the log scale, i.e. 95% of the mass between 0.11 and
# 0.80 s. A tighter 0.2 (0.20 to 0.44 s) was tried first and was tight enough
# to hurt: on the warm-start benchmark's `shifted` target, whose non-decision
# time is about 0.6 s - 3.5 SDs from the centre under that prior - the RDM
# produced divergent transitions in every run, the reference included, and
# widening the SD to 0.5 removed them (Rhat 1.06 and a minimum ESS of 67
# became 1.01 and 392; benchmarks/rdm_brittleness_report.md). The prior's job
# is to fence the `ndt -> 0` direction, where the likelihood goes flat and a
# flat prior would make the posterior improper; at 0.5 it still does that
# (log(0.01 s) is 6.8 SDs out) without telling the data where in 0.1 to 0.8 s
# the non-decision time is.
#
# poutlier is deliberately NOT the same belief twice. Leaving it out of the
# formula is itself information: the user either trimmed the data already or does
# not expect outliers. So the `nat` form puts its mode at zero, which a
# logit-scale prior cannot do at any location. The median matches regardless -
# exponential(100) has median 0.0069 against plogis(-5) = 0.0067 - so the centre
# is unchanged, only the shape. Truncation at 1 is negligible: exp(-100) of the
# mass is above it. brms's own default here is flat over [0, 1] - proper, but it
# puts half its mass above 0.5.
#
# shape exists only for cogmod_loggamma() and has an identity link, so the two
# scales carry the same prior. It is centred on the LogNormal and tight enough
# that sigma * shape = 1 - where the density blows up at the shift - is several
# SDs away for any realistic sigma. It is overridden rather than filled: brms
# recognises the name from its own gamma/weibull families and supplies
# gamma(0.01, 0.01), a prior with support only on the positives, for a parameter
# that has to be free to go negative (that is the whole inverse Weibull half of
# the family). Stan would reject every negative proposal, silently truncating
# `shape` at 0.
#' @keywords internal
.SHIFTED_BASE_PRIORS <- list(
  ndt = c(link = "normal(-1.2, 0.5)", nat = "lognormal(-1.2, 0.5)"),
  poutlier = c(link = "normal(-5, 1)", nat = "exponential(100)"),
  shape = c(link = "normal(0, 0.5)", nat = "normal(0, 0.5)")
)


# Priors for the shifted families built on `ndt` + `poutlier`: the two the
# parameterization introduces, plus whatever decision dpars the family named for
# itself in the registry's `prior` slot.
#' @keywords internal
.priors_shifted <- function(formula, data, family, ...) {
  own <- .mixture_spec(.family_name(family))$prior
  if (is.null(own)) own <- list()
  # The base entries win, but no family names one of them, so this only decides
  # a collision that cannot currently happen.
  own[names(.SHIFTED_BASE_PRIORS)] <- .SHIFTED_BASE_PRIORS
  .priors_dpars(formula, data, family, own, override = "shape", ...)
}


# Reads get_prior() for the model in hand and fills only the rows brms left
# flat - plus the ones named in `override`, where the default brms supplies is
# worse than nothing - so the result cannot contain a prior matching no
# parameter.
#
# `own` is a named list, one entry per dpar to edit, each a character vector
# with any of `link`, `nat` and `slope`. See .PRIORS_PLAIN above.
#' @keywords internal
.priors_dpars <- function(formula, data, family, own, override = character(0),
                          ..., slope_default = "normal(0, 0.2)") {
  p <- brms::get_prior(formula, data = data, family = family, ...)
  # Kept whole: once `p` has been filtered down to the rows being set, there is
  # no way left to tell "this blanket row covers a coefficient we left alone"
  # from "this blanket row covers nothing at all".
  all_rows <- p

  dpars <- names(own)
  if (is.null(override)) override <- character(0)

  # A dpar reaches get_prior() in one of two forms, depending on whether it
  # appears in the formula at all.
  #
  #  - *Modelled*, even as `~ 1`: brms builds a linear predictor for it and
  #    reports class "Intercept" / "b" / "sd" with dpar = "<name>". The
  #    parameter lives on the LINK scale (log for ndt, logit for poutlier,
  #    softplus for the ex-Gaussian's sigma and tau).
  #
  #  - *Omitted*: brms declares it as a plain auxiliary parameter - the same
  #    mechanism as `sigma` for gaussian() when you do not write `sigma ~ .` -
  #    and reports class = "<name>" with dpar = "". It is declared with the
  #    family's own bounds and NO link applied, so it lives on the NATURAL
  #    scale and needs a different prior, not the same one.
  #
  # Matching on `dpar` alone caught only the first kind, which left an omitted
  # `poutlier` on brms's flat prior over [0, 1] and an omitted `shape` genuinely
  # improper. Auxiliary rows are taken whether or not brms filled them, which is
  # what replaces the uniform(0, min_Y) it puts on an omitted `ndt` - precisely
  # the min-RT bound this parameterization exists to remove - and the
  # student_t(3, 0, 2.5) it puts on an omitted ex-Gaussian `sigma`.
  aux <- p$class %in% dpars & !nzchar(p$dpar)

  link <- p$dpar %in% dpars & !nzchar(p$prior)
  # ...unless a proper blanket row already covers them, as brms gives smooths.
  link <- link & !vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  # The modelled counterpart of the auxiliary override: these rows arrive
  # non-empty, so they have to be replaced rather than filled.
  link <- link | (p$dpar %in% override & p$class == "Intercept")
  # A smooth's wiggliness scale (`sds`) arrives the other way round from a
  # group-level `sd`: brms fills the BLANKET row itself, student_t(3, 0, 2.5),
  # and leaves the per-term rows empty. The empty rows are then covered and the
  # filled one was never a candidate, so until 0.3.3 an `sds` on a dpar kept
  # brms's default while the table in ?cogmod_priors said exponential(1). (The
  # `sd` rows escape because brms's blanket there has an empty group and the
  # rows it sets have the grouping factor's, so the two never match.) The
  # blanket is the row brms will use, so take it and replace it. On a logit- or
  # log-linked dpar a half-t(3, 0, 2.5), median 1.9 on the link scale, lets
  # the smooth alone walk a `sigmabias` across its whole range or move a
  # `sigmandt` by a factor of seven - which undoes the tight intercept the
  # family put there on purpose. The response's own smooth (`dpar == ""`) is
  # not touched, for the same reason its slopes are not.
  link <- link | (p$dpar %in% dpars & p$class == "sds" & !nzchar(p$coef))

  # `mu` is the response's own linear predictor, so brms reports it with an
  # EMPTY dpar - class "Intercept", or class "b" with coef "Intercept" under
  # `0 + Intercept` - and neither branch above can see it. A family that names a
  # `mu` prior in its registry slot is asking for that row specifically, and
  # only cogmod_exgaussian() does: its `mu` is on an identity link, where brms's
  # student_t(3, 0, 2.5) is centred on zero seconds. Every other family leaves
  # `mu` to brms and is unaffected, because none of them names one.
  #
  # Only the INTERCEPT is taken, never the slopes. The response's slopes are the
  # effects the user is there to estimate; shrinking them with a blanket prior
  # is a different decision from giving the intercept a sensible location, and
  # not one this function should make on its own.
  resp <- ("mu" %in% dpars) & !nzchar(p$dpar) &
    (p$class == "Intercept" | (p$class == "b" & p$coef == "Intercept"))

  fill <- aux | link | resp
  if (!any(fill)) {
    return(brms::empty_prior())
  }
  # The dpar to look the prior up under. It is `p$dpar` for everything except
  # the response rows above, whose dpar is empty by construction.
  target <- p$dpar
  target[resp] <- "mu"
  target <- target[fill]
  aux <- aux[fill]
  p <- p[fill, , drop = FALSE]

  p$prior <- vapply(
    seq_len(nrow(p)),
    function(i) {
      cls <- p$class[i]
      # Natural scale: no link is applied to an auxiliary parameter.
      if (aux[i]) return(.own_prior(own, cls, "nat"))
      intercept <- cls == "Intercept" || (cls == "b" && p$coef[i] == "Intercept")
      if (intercept) {
        .own_prior(own, target[i], "link")
      } else if (cls == "b") {
        # a family may widen the blanket slope prior for a dpar of its own
        slope <- .own_prior(own, target[i], "slope")
        if (nzchar(slope)) slope else slope_default
      } else if (cls %in% c("sd", "sds")) {
        "exponential(1)"
      } else {
        "" # cor, and anything else brms already gives a proper default
      }
    },
    character(1)
  )
  p <- p[nzchar(p$prior), , drop = FALSE]

  # get_prior() reports a blanket row (empty coef) alongside one row per
  # coefficient. Setting both makes brms warn that the blanket one is unused, so
  # keep the blanket row and drop what it subsumes. `b` on a coefficient named
  # "Intercept" is the exception: under `0 + Intercept` that IS the intercept and
  # carries a different location from the slopes.
  drop <- vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  drop <- drop & !(p$class == "b" & p$coef == "Intercept")
  kept <- p[!drop, , drop = FALSE]

  # The exception can empty the blanket row out. `ndt ~ 0 + Intercept` has one
  # coefficient, that coefficient is the intercept, and it has just been given a
  # location of its own - leaving the blanket row covering nothing, which brms
  # warns about in turn. So the pair is checked in both directions.
  gone <- vapply(seq_len(nrow(kept)), .blanket_now_unused, logical(1),
                 p = kept, all_rows = all_rows)
  kept[!gone, , drop = FALSE]
}


# One entry from a family's own `prior` slot in the registry, or "" if the
# family named no prior for that dpar, or named one but not on that scale. The
# slot holds named character vectors, and `[[` on one of those errors rather
# than returning NULL when the name is absent, so both lookups are guarded.
#' @keywords internal
.own_prior <- function(own, dpar, what) {
  if (is.null(own) || !dpar %in% names(own)) return("")
  v <- own[[dpar]]
  if (!what %in% names(v)) return("")
  unname(v[[what]])
}


# Is row i of a get_prior() table subsumed by a blanket row (one with an empty
# `coef`) that already carries a prior? brms warns about the blanket row going
# unused if both are set, so only one of the pair should be filled.
#' @keywords internal
.covered_by_blanket <- function(i, p) {
  if (!nzchar(p$coef[i])) {
    return(FALSE)
  }
  any(
    !nzchar(p$coef) & nzchar(p$prior) & p$class == p$class[i] &
      p$dpar == p$dpar[i] & p$group == p$group[i] & p$resp == p$resp[i]
  )
}


# The mirror image: is blanket row i of the *final* table left with nothing to
# cover, because every coefficient it applies to kept an individual prior?
# `all_rows` is the unfiltered get_prior() table, which is the only place the
# full coefficient list still exists.
#' @keywords internal
.blanket_now_unused <- function(i, p, all_rows) {
  if (nzchar(p$coef[i])) {
    return(FALSE)
  }
  same <- function(d) {
    d$class == p$class[i] & d$dpar == p$dpar[i] & d$group == p$group[i] &
      d$resp == p$resp[i]
  }
  coefs <- all_rows$coef[nzchar(all_rows$coef) & same(all_rows)]
  # No coefficients to cover means this is a blanket over something else
  # entirely (a group-level SD, say), which brms is happy to use as it stands.
  length(coefs) > 0 && all(coefs %in% p$coef[nzchar(p$coef) & same(p)])
}
