# cogmod 0.3.3

## New features

* **A distributional parameter pinned at a boundary is now checked against the
  response.** The three bounded-scale families are mixtures of a continuous
  part and one or more point masses, and every weight involved can be fixed in
  `bf()`. Fixing one at a boundary switches its component off, so any
  observation belonging to that component has zero probability and the
  log-likelihood is `-Inf` everywhere - which CmdStan reports only as
  `Initialization failed after 100 attempts`, naming nothing. `cogmod_priors()`
  now errors first, saying how many rows are involved and how to proceed.
  Eleven combinations are covered: for `cogmod_choco()`, `pmid = 0` with a
  response exactly at the midpoint, `pmid = 1` with one anywhere else,
  `pex = 0` with an exact 0 or 1, `bex = 0` with a 1 and `bex = 1` with a 0;
  for `cogmod_betagate()` the same `pex = 0` and `bex` cases plus `pex = 1`
  with any interior response; and for `cogmod_betadiscrete()`, `pzero = 0`
  with a 0 and `pzero = 1` with any rating. The first and the last are the
  ones the documentation invites - `?rcogmod_betadiscrete` offers `pzero = 0`
  as the way to say a scale has no zero category, and `?cogmod_priors`
  suggests fixing `pmid` or `pzero` at 0 to switch the parameter off.

* **`cogmod_priors()` now covers `cogmod_choco()`, `cogmod_betagate()` and
  `cogmod_betadiscrete()`.** Every distributional parameter of the three
  bounded-scale rating families previously arrived either flat - improper - or
  with a `brms` default aimed at a different parameterization. Two of those
  were worse than unhelpful. The point-mass probabilities `pmid` and `pzero`
  were flat on a `logit` link, which is improper in the posterior as well as
  the prior whenever the event is absent from the data: with no exact
  midpoints anywhere, the likelihood in `pmid` increases monotonically all the
  way to zero and nothing stops the logit running to minus infinity. That is
  `poutlier`'s failure, and they now get `poutlier`'s treatment, including an
  omitted form whose mode is at zero. Measured on 400 slider responses with no
  midpoints and no extremes: the flat defaults put `pex` at `-1.1e14` and
  `pmid` at `-3.4e13`, `Rhat` 2.9, ESS 5, with 12% divergent transitions;
  these priors give -4.70 and -5.15, `Rhat` 1.00, ESS 3000-4000, no
  divergences. `phi` arrived with
  `student_t(3, 0, 2.5)`, which `brms` supplies because it recognises the name
  from its own beta family - the same trap `cogmod_exgaussian()`'s `sigma`
  falls into; on `cogmod_betadiscrete()`'s `log` link that prior's 95%
  interval runs to `phi = 2853`, where the Beta has collapsed onto a single
  rating category, so it is overridden rather than filled. The Beta precisions
  are fenced away from both ends: below about 1 the underlying Beta is
  U-shaped and unbounded, and a large one narrows it sharply. `mu` is
  left to `brms` in all three, since on a `logit` link its
  `student_t(3, 0, 2.5)` is the standard weakly informative choice. See
  `?cogmod_priors` for the full table and the reasoning behind each number.

* **`cogmod_inits()` now covers `cogmod_choco()` and `cogmod_betadiscrete()`.**
  Both are bounded-scale families for subjective ratings, and neither has a
  flat region of the kind that makes an init mandatory for the RT families -
  every start is a proper density and a chain begun at the default does move.
  What it has to move away from is a description of a rating scale that no
  data set matches. The logit origin puts `pmid` and `pzero` at 0.5, i.e. half
  of every response exactly on the midpoint or outside the scale altogether;
  they now start at 0.05. `cogmod_choco()`'s two Beta precisions are behind a
  `softplus` link, so `softplus(0) = 0.69` puts both Beta shapes below 1 - a
  U-shaped rating distribution, unbounded at both ends of each half of the
  scale; they now start at 2, in the middle of a plateau where the fit is
  insensitive to the exact value. `pex` starts at 0.1 rather than 0.5. On a
  2000-trial slider data set the old start sat about 1900 log-likelihood units
  from the new one, 983 of them `pmid` alone. For `cogmod_betadiscrete()`,
  `phi = 1` with `mu = 0.5` is the discrete Uniform, which is where the
  origin already is - it is named so that `init = "random"`, which draws `phi`
  anywhere from 0.14 to 7.4, no longer decides it. So this buys warmup rather
  than a fit that would otherwise fail.

* **`cogmod_inits()` lists `cogmod_geg()` among the families it supports.** It
  had targets for it but was not naming it, so the error message for an
  unsupported family said cogmod_geg() was one. The supported list is now
  derived from the targets rather than kept beside them.

* **`cogmod_inits()` jitters the hierarchical blocks a fifth as much as the
  population-level ones, and starts smooths flat.** One `jitter` (default
  0.25) now applies to intercepts and slopes; the standardized group-level
  effects `z_*`, the smooth coefficients `zs_*` and their scales `sd_*` and
  `sds_*` get a fifth of it (0.05), and `sds_*` starts at 0.05 rather than the
  generic 0.25. Two numbers set the two tiers directly. The reason is that a
  unit of noise on a group effect or a spline coefficient is not a unit on the
  linear predictor: it is multiplied by a scale and a design column - tensor
  basis values reach tens - once per participant or basis function. On a
  production model with a tensor smooth on five distributional parameters and
  a participant intercept on six, the old jitter moved linear predictors by
  one to two and a half link units at some rows, which started one
  participant's non-decision time above 98 of their 128 trials and a sigma at
  0.07 s where 0.5 was intended, all of it then explained by the outlier
  component; at 0.05 on those blocks the same model started where the targets
  say. Chains still start apart where it matters for Rhat, on the intercepts
  and slopes. Warm starts (`cogmod_warmstart()`, `warmstart =`) use the same
  rule at their smaller default.

* **The 7-parameter `cogmod_ddm()` is 1.7x cheaper per gradient when both
  `sigmabias` and `sigmandt` are estimated.** Stan's 7-parameter
  `wiener_lpdf()` integrates the start-point and non-decision-time ranges out
  by adaptive cubature - eight integrals per observation, one for the density
  and one per partial derivative - and takes a tolerance for them that the
  Stan code never passed, so it ran at Stan's default of 1e-4. It now passes
  `1e-3` (`.DDM_WIENER_PRECISION`). Measured per observation and gradient on
  CmdStan 2.38.0: with both ranges open, 274 us falls to 163; with one open
  the cost does not move (102 to 100 us), because a one-dimensional integral
  is already settled by cubature's first pass, and nothing much happens past
  1e-3 either (153 us at 1e-2), which is why it stops there. On the
  intercept-only program `benchmarks/gradient_cost.R` times, 5000 trials with
  all seven parameters estimated, a gradient went from 3.25 s to 1.50 s
  (ratio 0.46, median of 21 alternating blocks). The price: the
  Stan density still agrees with the R one to 1e-5 over the test grid, but
  its gradient sits about 2e-5 relative from central differences at a typical
  point where it sat 2e-6 before, and 1.2e-4 with `ndt` pushed 2.5 log units
  above its start - a point no chain visits after warmup, and small next to
  the leapfrog error a sampler already carries, but a change of an order of
  magnitude, so it is recorded with the constant. Any value between 3e-4 and
  1e-3 measures the same gradient error and 1e-3 is the cheaper. A paired fit
  settles it (`benchmarks/ddm_precision/`): the 7-parameter model on 100
  trials, both tolerances from one start, metric and seed, end to end and
  again with adaptation held fixed - the same 24 leapfrog steps per
  iteration, the same step size, acceptance 0.90 against 0.92, fewer
  divergences (7 against 10 and 11), min bulk ESS 248 against 206 and 261
  against 224, at half the cost per gradient (445 against 883 ms): 2.2 to
  2.4 times the effective samples per CPU second. That fit also puts the
  standalone figures above in perspective: at its posterior a gradient cost
  8.3 ms per observation at Stan's default tolerance, thirty times the
  274 us of the benign benchmark point. The same cost measurements are now
  in `?cogmod_ddm`, under the between-trial variability section:
  estimating `sigmadrift` costs about 2.8x the classic model per gradient,
  estimating one of `sigmabias` / `sigmandt` about 18x and both about 30x,
  and the fast path is a test for *exactly* zero, so a tight prior does not
  buy it back - only `sigmandt = 0` in `bf()` does. Nothing user-facing said
  so before, and it is the difference between a day and weeks on a large
  data set.

## Bug fixes

* **`cogmod_priors()` now sets the wiggliness prior (`sds`) of a smooth on a
  distributional parameter.** The table in `?cogmod_priors` has always listed
  `sds` alongside `sd` at `exponential(1)`, but brms fills a smooth's blanket
  `sds` row itself and leaves the per-term rows empty, so a function that
  fills only what arrives empty never reached it: `ndt ~ s(x)` kept
  `student_t(3, 0, 2.5)` on its smooth's scale. On a log or logit link that is
  a half-t with median 1.9 in link units - loose enough for the smooth alone
  to move a `sigmandt` by a factor of seven or walk a `sigmabias` across its
  range, behind an intercept the family had deliberately fenced. The blanket
  row is now replaced, for every family and every dpar `cogmod_priors()`
  handles; the response's own smooth is left to brms, as its slopes are.
  Group-level `sd` rows were never affected.

* **`cogmod_rdm()`'s gradient is now exact.** Every normal tail in the RDM's
  Stan code went through Stan's `std_normal_lcdf()`, whose value is right but
  whose partial derivatives are an approximation, and in a race those partials
  *are* the gradient of the drifts and the boundary. Measured against central
  differences of the log probability, the gradient sat 2e-4 relative from the
  truth at a typical start and as far as 7e-2 where the drift is small, while
  the log probability itself was smooth to 1e-7. HMC stays exact under an
  inexact gradient - the accept step corrects for it - but pays in step size
  and acceptance. All eleven calls now go through `cogmod_log_Phi()`, the
  function introduced below for the LNR, which brings the worst error over the
  same grid to 2e-7 at about 14% more per gradient evaluation. The function
  moved into a prelude of its own (`.LOG_PHI_STAN_PRELUDE`) so that any family
  can take it; the LogNormal's and the RDM's preludes both start with it.

* **`cogmod_invgaussian()`'s gradient is now exact.** The same defect as the
  RDM's, in the family's own Stan code: every normal tail of the Wald - the two
  terms of its CDF and survival, the two integrated ones a non-decision range
  needs, and the truncation factors that make `sigmadrift` a *truncated* normal
  drift - went through `std_normal_lcdf()`. `sigmadrift` is differentiated
  almost entirely through those truncation factors, and its partial sat 1.2e-3
  relative from central differences of the log probability, `mu`'s 2.7e-4,
  against 1e-8 for the families already on `cogmod_log_Phi()`. Six of the
  eighteen points the gradient check walks failed on it. All of them now go
  through `cogmod_log_Phi()`, which brings the worst error over the same grid
  to 5e-8; the truncation factors lose a round trip through `log1m_exp()` in
  the bargain, since `log1m_exp(std_normal_lcdf(-z))` is `cogmod_log_Phi(z)`
  written the long way. `cogmod_exwald()` builds on the same prelude and gets
  the same fix. Values are unchanged.

* **`cogmod_geg()`'s CDF term is no longer Stan's `exp_mod_normal_lcdf()`.**
  The alpha-power construction puts `(shape - 1) * log F_EG` inside the
  *density*, so the ex-Gaussian CDF is differentiated on every evaluation, not
  only when a response is censored. The built-in agrees with the R side to
  1e-8 over the range the tests walk, but outside it - at `sigma / tau` around
  9, where the two terms of `F_EG` cancel hardest - both its value and its
  partials go wrong: Stan's own central differences put `d/d sigma` at 28.5
  where its autodiff said 77.5, and `d/d tau` at 113.6 against 53.0. The Stan
  side now subtracts the two terms in log space itself, through
  `cogmod_log_Phi()`, exactly as `.lcdf_exgaussian()` always has in R; over the
  same grid the worst gradient error is 2e-9, and the log CDF matches
  quadrature of the density to 5e-15. It costs about 10% per gradient
  evaluation (3.1 ms to 3.5 ms over 5000 trials, median of 21 alternating
  blocks). Two things kept it to that. The ex-Gaussian density is the CDF's
  second term over `tau`, so the GEG writes both out from one pair of normal
  tails rather than calling for each; and the pair is written inline rather
  than returned from a helper - a `vector[2]` return measured half again the
  cost of the inline form, through the autodiff-stack allocation it makes on
  every call.

* **`cogmod_exgaussian()`'s Stan density, CDF and survival go through
  `cogmod_log_Phi()`.** The density is the same expression as before -
  `exp_mod_normal_lpdf()` and `dcogmod_exgaussian()` in R evaluate it too, and
  log-probabilities are unchanged to every digit that matters - but it is now
  written out over `cogmod_log_Phi()`, so it stays finite about 38 standardized
  units into the left tail where the built-in's bare `erfc` underflows, and it
  shares its normal tail with the CDF that `cogmod_geg()` is built on. The
  `cens()` path gets the accurate CDF with it: `cogmod_exgaussian_lcdf()` was
  `exp_mod_normal_lcdf()` and had the same wrong partials as the GEG's, and
  `cogmod_exgaussian_lccdf()` was on `std_normal_lcdf()`. About 10% per
  gradient evaluation, the cost of `cogmod_log_Phi()` over the built-in.

* **`cogmod_lnr()` and `cogmod_lognormal()` no longer hand Stan a non-finite
  gradient in the tails.** `cogmod_lognormal_ldiff_Phi()` formed
  `log(Phi(y + c) - Phi(y))` from `erfc` as `log(u1) + log1m(u2 / u1)`, and the
  `A == 0` branch of `cogmod_lognormal_acc_ltails()` called `lognormal_lcdf()`
  and `lognormal_lccdf()`, which are `erfc` alone. `erfc` underflows near
  `x = -38` - about 38 standardized log units from an accumulator's median
  finishing time - and past that the value is `log(0)` and the partials are
  `inf` or `0 / 0`. The outlier mixture then hides it: `log_mix()` stays finite
  with one component at `-inf`, but reverse mode multiplies the zero adjoint
  into the stored partial and `0 * inf` is `NaN`, so one response in a data set
  turned the gradient of the whole model to `NaN`. That reads as
  `Gradient evaluated at the initial value is not finite` at the start of a fit
  and as divergent transitions afterwards. Both now go through one function,
  `cogmod_log_Phi()`: `erfc` in the body of the distribution and, below
  `x = -25`, the asymptotic expansion of the tail, whose leading term is the
  exponent itself, so nothing underflows and the result stays finite and
  differentiable as far as `x = -1e150`; its six terms agree with R's
  `pnorm(log.p = TRUE)` to 4e-16 relative, so the two branches meet with no step
  in the density. `cogmod_lognormal_ldiff_Phi()` now takes the difference in
  logs with `log1m_exp()` rather than as a quotient of two minute numbers.
  Measured over a grid of decision times from 1 ms to 300 s and sigmas from
  0.02 to 1.2, with and without a start-point range: 6 of 72 gradients were
  non-finite before and none are now, and the densities agree with the R
  kernels everywhere, with no value newly truncated to `-inf`.

  `std_normal_lcdf()` is not used for this, although it has the range - its
  value is exact against `pnorm(log.p = TRUE)` as far as `x = -1e7`. Its
  analytic partials are not: on 20000 responses they sat 1.7e-3 from central
  differences of the log probability where the `erfc` route sat 4e-6, and in a
  race those partials *are* the gradient of `nu` and `sigma`.

  Where the old code's gradient was finite it was not always right. On 20000
  responses with a start-point range it sat 8.3e-4 from central differences and
  the new one sits 9.6e-7, which is the finite differences' own noise; without
  a start-point range both sit at 4.1e-6, so nothing there was given up for it.

  Sampling is not slower in the case most models are in.
  `cogmod_lognormal_acc_logcdf()` and `cogmod_lognormal_acc_logsurv()` now take
  the single tail they were asked for when `sigmabias = 0`, instead of building
  the pair and discarding one: the two share no work there, and `cogmod_lnr()`
  reads the survival alone. On 20000 responses a gradient of the plain LNR came
  out about 20% cheaper than before and one with a start-point range about 15%
  dearer, the latter buying the corrected gradient above.

## Documentation

* Help pages that sent the reader to `vignette("rt_models")` or
  `vignette("performance")` now link to those articles on the package website.
  Both are website articles, not installed vignettes, so the call failed.

* The `cogmod_warmstart()` example is now self-contained and runs (with
  CmdStan): it simulates its data, fits the pilot and the warm-started model,
  and writes the warm-start table to a temporary file rather than to the
  working directory.

# cogmod 0.3.2

## New features

* **`cogmod_lognormal()` gains `sigmabias`**, the same between-trial
  start-point range as the LNR's below: the decision time is the LogNormal
  multiplied by a `Uniform(1, 1 + sigmabias)` distance, which is the
  single-accumulator LBA with a LogNormal drift rate and its threshold offset
  pinned at 1. At `sigmabias = 0` - the default of `rcogmod_lognormal()`,
  `dcogmod_lognormal()` and `pcogmod_lognormal()`, placed after `ndt` and
  before `poutlier` as the Wald's `sigmadrift` is, and the value to fix in the
  formula unless
  the design speaks to start-point variability - the family is the shifted
  LogNormal exactly as before, bit for bit and at the same cost. The CDF and
  survival `cens()` needs, the mean `posterior_epred()` reports (which gains a
  factor `1 + sigmabias / 2`) and the Stan functions all carry the range, and
  agree with quadrature over the start point to `1e-7` and with each other to
  `1e-10`. The kernels are shared with `cogmod_lnr()`, which is now a race of
  two of these accumulators in code as well as in theory. As for the LNR, a
  formula that omits `sigmabias` estimates it, `cogmod_priors()` fences the
  flat direction at zero, and fits made with earlier versions have to be
  refit; the two vignette LogNormal models were. `cogmod_logstudent()` and
  `cogmod_loggamma()` do not get the parameter: the density needs a partial
  first moment of the rate distribution, which does not exist for a Student-t
  on the log scale and needs incomplete gamma functions for the log-Gamma.

* **`cogmod_lnr()` gains `sigmabias`**, a between-trial start-point range:
  each accumulator now starts at `Uniform(0, sigmabias)` and runs to a
  threshold `1 + sigmabias` at its LogNormal rate, so its finishing time is
  the distance divided by the rate rather than the reciprocal of the rate
  alone. At `sigmabias = 0` - the default of `rcogmod_lnr()` and
  `dcogmod_lnr()`, and the value to fix in the formula unless the design can
  identify a start-point range - the family is the LNR exactly as before, bit
  for bit and at the same cost, since both the R and the Stan kernels take the
  plain lognormal branch there. Above zero it is the LBA with LogNormal drift
  rates ([Heathcote & Love, 2012](https://doi.org/10.3389/fpsyg.2012.00292)),
  the model the LNR was introduced as a limit of: with a LogNormal *distance*
  as well as a LogNormal rate the two fold into one `sigma`, which is why the
  LNR never had a start-point parameter, whereas a Uniform distance leaves a
  shape the rate alone cannot produce. The threshold offset is pinned at 1
  rather than `sigma` at 1 because rescaling the evidence axis shifts `nu` and
  scales the range and the threshold but leaves a LogNormal rate's `sigma`
  untouched, so `sigma` cannot pin the scale; `sigmabias` is therefore read in
  units of the threshold offset. The density is a difference of two normal
  CDFs and the survival one more, evaluated from a series below a start-point
  range of about `1e-4 * sigma` and from ratios of log-CDFs in the tails, and
  both agree with one-dimensional quadrature over the start point to `1e-7`.
  `cogmod_priors()` fences the flat direction at zero the way it does for
  `cogmod_lba1()`'s `sigmabias`, and `cogmod_inits()` starts it at 0.5.
  A formula that omits `sigmabias` now estimates it, as `brms` does with any
  dpar; fits made with earlier versions cannot be post-processed, because
  their family carries no `sigmabias`, and have to be refit. In
  `rcogmod_lnr()` and `dcogmod_lnr()` the argument sits between `ndt` (or
  `response`) and `poutlier`, where the DDM keeps its between-trial
  variabilities, so a call that passed `poutlier` by position needs it named.
  The vignette LNR model pins it at zero and was refit.

* **`cogmod_priors()` gains `warmstart`**, which re-centres the priors on a
  previous fit: every population-level intercept and coefficient, every
  group-level SD, and every dpar left out of the formula gets
  `normal(median, prior_scale * sd)` from that fit's posterior, matched
  parameter by parameter on the `class`, `dpar`, `coef` and `group` a
  `get_prior()` row carries. No transformation is involved - the parameter a
  prior row is about is the parameter the source sampled, the `Intercept`
  prior being stated on the centred intercept in both - so the scales line up
  by construction. `cogmod_warmstart()` extracts the posterior median and SD
  alongside the means for this, `as.data.frame()` carries them in two new
  columns, and a table written before they existed still reads (the priors are
  then left alone, with a message). The correlations keep their LKJ and the
  standardized effects have no stated prior to change.
  **Unlike the rest of the warm start, this changes the posterior**, and if
  the source was fitted to data the new model also contains it double-counts
  it - a pilot on half the participants used to centre the priors for the fit
  on all of them uses that half twice. `prior_scale`, 3 by default, is what
  stands between a prior that only says roughly where the parameter lives and
  one that is the source's posterior outright; `?cogmod_priors` says when the
  argument is and is not legitimate.

## Bug fixes

* **`cogmod_rdm()` no longer freezes one chain in four on a cold start.**
  `cogmod_inits()` started both drifts at 3; the error accumulator now starts
  at 1. The failure looked like a stuck chain - every transition at the
  maximum treedepth, step size a thousand times smaller than the other
  chains', Rhat 1.5 to 2.9 - and it was traced to the very first warmup
  transition. A Wald density is thin on the fast side and flat on the slow
  side, so on data whose error drift is about 0.2 (the speed-accuracy data of
  `vignette("performance")`) the old start sat 400 log-density units above
  the posterior. The first trajectory converted that into momentum along the
  flat `driftone` direction (the plateau `cogmod_priors()` fences, where the
  likelihood no longer changes) and carried the chain from a link value of +3
  to -20 in eight leapfrog steps; the step size then collapsed to 1e-5 within
  a dozen iterations and the metric windows that followed were estimated from
  a chain that no longer moved. Stan's model methods found nothing numerical
  at the frozen position - the log-density and its gradient are smooth across
  the driftless branch of the survival function, and finite differences agree
  with autodiff - so the fix is where the chain starts, not what it computes.
  On an 800-trial mixed model with a 200-iteration warmup the old start
  froze a chain in 1 run in 4 with Stan's default metric and in 5 of 6 when
  handed a metric adapted to the bulk (4 of 4 when the benchmark cell was
  rerun); the new start has done so in 0 of 6 under the latter, the harder
  case, and each such fit ran in 3 minutes instead of 30. When the error accumulator really is as fast as the correct
  one the start is off by a factor of three on the cheap side, which costs a
  few dozen units and changes nothing.

* **`cogmod_inits()` now starts `ndt` at half the first percentile of the
  observed response times instead of a fixed 0.1 s.** Same mechanism as the
  previous item, other cold start. The fixed value was a third of the prior
  median and safely below any ordinary data, which was the whole argument for
  it; on data whose non-decision time is 0.6 s it sat half a second low, every
  decision time looked far too long, a driftless race then fit better than a
  fast one, and the first trajectory of a cold `cogmod_rdm()` chain threw
  both drifts onto their flat regions - three chains in four on the
  benchmark's `shifted` target when handed a metric adapted to the bulk. Half
  the first percentile is still below essentially every response, so the
  gradient the small start was protecting is intact, and it follows the scale
  of the data. Every `ndt` + `poutlier` family gets it. Starting values do
  not change a posterior, so no fitted model needs revisiting for this.

## Breaking changes

* **The default `ndt` prior is wider: `normal(-1.2, 0.5)` on the log scale,
  `lognormal(-1.2, 0.5)` for an omitted `ndt`, in place of the `0.2` SD.**
  The centre is unchanged at 0.30 s; 95% of the mass now sits between about
  0.11 and 0.80 s instead of 0.20 to 0.44 s. The old SD put a 0.6 s
  non-decision time - not unusual for older participants or more demanding
  responses - 3.5 SDs from the centre, and on the warm-start benchmark's
  `shifted` target (`vignette("performance")`), whose non-decision time is
  about 0.6 s, `cogmod_rdm()` produced divergent transitions in every run
  with that prior, the 600-warmup reference included; the wider prior removed
  them (Rhat 1.06 and a minimum ESS of 67 became 1.01 and 392 on the same
  seed). The prior's job is to fence the `ndt -> 0` direction, where the
  likelihood goes flat and a flat prior would make the posterior improper; at
  0.5 it still does (0.01 s is 6.8 SDs out) without telling the data where in
  0.1 to 0.8 s the non-decision time is. Fits with the default priors will
  move slightly, most where the data put `ndt` far from 0.30 s, which is
  where they should have been free to move. Every `ndt` + `poutlier` family
  is affected.

# cogmod 0.3.1

## New features

* **New `cogmod_warmstart()`** turns what a previous fit's warmup produced -
  the adapted inverse metric, the step size and the posterior means - into the
  `init`, `inv_metric` and `step_size` arguments of a new `brm()` call, so
  that a refit, or the same model on more participants, can run a much shorter
  warmup. Stan adapts one variance per unconstrained parameter, and the
  function labels each with its Stan name (read off the generated program with
  the parser `cogmod_inits()` already has) and joins the two models on those
  names: population-level entries carry over, a pilot participant's
  standardized effects follow it by level name to its position in the bigger
  model, new participants take their effect's average variance and start at
  zero, and anything without a counterpart gets Stan's default variance and a
  generic start, with a count in `print()`. The standardized effects and
  Cholesky factors that `brms` drops from a saved fit are rebuilt from the
  `r_`, `sd_` and `cor_` it keeps. `as.data.frame()` gives a table of a few
  kilobytes that survives `write.csv()` and can be passed back as a file path,
  so a pilot fitted on a laptop can warm-start an array job on a cluster. On
  the other side, `cogmod_inits()` gains a `warmstart` argument and the new
  `cogmod_inv_metric()` and `cogmod_step_size()` share its signature - the
  model's formula and data, then the source, be it a fit, a table or a file -
  so that each argument of `brm()` has one helper and the table is mapped onto
  the model on the way. Whatever of `formula` and `data` is not given is taken
  from the source fit. Works for any `brms` model fitted
  with the `cmdstanr` backend and the diagonal metric. On a mixed LNR and a mixed DDM, a pilot on 4 of 8
  participants warm-started the full fit to about twice the effective draws
  per second of a cold start with the full warmup, and four to six times those
  of a cold start with the same short warmup; the starting values alone bought
  nothing, so the metric and step size are the product
  (`vignette("performance")`).

* **`cogmod_invgaussian()` gains `sigmandt`**, the between-trial range of the
  non-decision time (`st0`): each trial's non-decision time is drawn from
  `Uniform(ndt, ndt + sigmandt)`, so `ndt` becomes its lower bound, exactly as
  in `cogmod_ddm()`. Spreading the shift turns the Wald density into a
  difference of two CDFs and the CDF into a difference of two integrated CDFs,
  both closed form at a fixed drift, so the parameter costs a few normal CDFs
  per observation, works with `cens()` unchanged, and rides the existing drift
  quadrature when `sigmadrift > 0` too. It is on a `log` link with
  `cogmod_ddm()`'s prior for the same quantity. **It is hard to estimate and
  should be fixed at zero for most applications** (`sigmandt = 0` in `bf()`);
  it shares the leading edge of the distribution with `ndt` and `poutlier`, and
  should only be freed with a lot of data, a strong prior, or both. As with
  `sigmadrift`, leaving it out of `bf()` *estimates* it, so existing Wald
  formulas that do not mention it now fit one more parameter unless they add
  `sigmandt = 0`, and fits made before this version cannot be post-processed
  with it; the vignette models were refit. `rcogmod_invgaussian()`,
  `dcogmod_invgaussian()` and `pcogmod_invgaussian()` take `sigmandt` right
  after `sigmadrift`, so a `poutlier` passed by position moves along one.

* **Censored reaction times: `brms`'s `cens()` works on the RT-only families.** `bf(rt | cens(error) ~ ...)` scores an error trial as a
  *right-censored correct response*: its RT is a lower bound on when the
  correct process would have finished, so it contributes that process's
  survival rather than its density. On `cogmod_invgaussian()` this is the
  *simple* censored shifted Wald of [Miller et al. (2018)](https://doi.org/10.1177/0146621617710465), their Eq. 4, the
  `version = "simple"` of the `cswald` model in [`bmm`](https://github.com/popov-lab/bmm).
  It is not their competing-risks variant (Eq. 5, a race of two Wald
  accumulators with drifts `v` and `-v`, as implemented in `rtdists` and
  `bmm`'s `version = "crisk"`), which is a choice model rather than a censoring
  construction; `cogmod_ddm()` with `bias` fixed at 0.5 covers that ground.
  Here censoring is not a family but a
  construction, so the same formula works on every RT-only family with a
  closed-form CDF: `cogmod_lognormal()`, `cogmod_logstudent()`,
  `cogmod_gamma()`, `cogmod_invgamma()`, `cogmod_weibull()`,
  `cogmod_invweibull()`, `cogmod_logweibull()`, `cogmod_bisa()`,
  `cogmod_exgaussian()` and `cogmod_geg()`. Left- and interval-censoring come
  with it, and `log_lik()` - hence `loo()` - honours all three, which `brms`
  leaves to a custom family's own method. `posterior_predict()` predicts the
  latent, uncensored RT, as `brms` does for its own families.
  `cogmod_priors()` and `cogmod_stanvars()` refuse `cens()` on the families
  that cannot take it, and `cogmod_priors()` warns above 20% censored trials -
  a threshold that is exact for timeouts and omissions at any rate but lenient
  for commission errors, where the construction is biased well before it (see
  `?rcogmod_invgaussian`).
  What the model is for, what it assumes and the one check to run before using
  it are in `?rcogmod_invgaussian` and the *Censored Shifted Wald* section of
  `vignette("rt_models")`.

  Under the hood every censorable family gets a `<family>_lcdf` and a
  `<family>_lccdf` beside its `_lpdf`, generated from two new registry slots so
  a family cannot drift out of step with itself. The survivals are written as
  survivals - never as `log(1 - exp(lcdf))` - and the half Normal outlier's
  through `std_normal_lcdf(-z)` rather than `std_normal_lccdf(z)`, which is
  `-inf` from 1.66 s on: the two places `bmm`'s implementation broke. With
  `sigmadrift > 0` the Wald CDF has no closed form and is taken by 64-point
  Gauss-Legendre quadrature over the drift, in R and Stan alike off one node
  table.

* **`pcogmod_*()` for every censorable family.** `pcogmod_lognormal()`,
  `pcogmod_logstudent()`, `pcogmod_gamma()`, `pcogmod_invgamma()`,
  `pcogmod_weibull()`, `pcogmod_invweibull()`, `pcogmod_logweibull()`,
  `pcogmod_bisa()` and `pcogmod_exgaussian()` join `pcogmod_invgaussian()`,
  with `lower.tail` and `log.p`. The upper tail is computed *as* the upper
  tail rather than as `1 - CDF`; these are the R side of the Stan
  `_lcdf`/`_lccdf` pair, and the tests hold the two to each other.
  
* **`cogmod_priors()` now checks the response before returning.** Everything in
  this package is stated in seconds and none of it is unit-equivariant - the
  `ndt` prior means 170-300 ms, `.POUTLIER_SCALE` is 0.2 s - but `brms` fills
  its own defaults from the data, so a column of milliseconds produces a model
  whose two halves silently describe different quantities. It compiles, it
  samples, it converges, and the estimates are meaningless. The check catches
  that and the handful of other mistakes with the same character.

  It **stops** where the offending rows would make the fit impossible or wrong
  in a way `Stan` cannot report: a non-positive reaction time under a family
  that places no density below `ndt`; a response outside `[0, 1]` for
  `cogmod_choco()` or `cogmod_betagate()`; a non-integer rating for
  `cogmod_betadiscrete()`; a non-numeric response; and a third level in `dec()`,
  which the choice families would otherwise fold silently into option 1, since
  their Stan code tests `dec == 0` and takes the else branch for everything
  else.

  It **warns** about the rest: a median implying milliseconds, `NA`s, and either
  tail running past what `poutlier` can absorb. The tails are judged as
  proportions rather than counts, because the outlier component is *supposed* to
  produce the occasional fast response - `rcogmod_lognormal(200, ndt = 0.2,
  poutlier = 0.02)` puts one at 81 ms - and a count-based test fires on the
  package's own generator. Over 20000 draws the component sends 0.8% of
  responses below 0.1 s at `poutlier = 0.02` and 1.9% at 0.05, the top of the
  default prior, so the warning sits at 5%.

  Families with neither `ndt` nor `poutlier` - `cogmod_exgaussian()` and
  `cogmod_geg()` - are exempted from the tail checks, and a non-positive
  response is a warning rather than an error there, their support being the
  whole real line. A formula or family the check cannot read is passed through
  untouched, so `brms`'s own error is what the user sees.

## Bug fixes

* **`dcogmod_ddm()` and everything built on it (`log_lik()`, `loo()`,
  `p_outlier()`, the other R-side post-processing of `cogmod_ddm()` fits) are
  now accurate when the `sigmandt` range reaches down to fast decision
  times.** The R density integrates the non-decision time out with a fixed
  25-node Gauss-Legendre rule, and when that range covers decision times from
  about zero up to the response - a fast response, or a wide `sigmandt` - the
  integrand holds the whole early peak of the first-passage density inside a
  sliver of it, which 25 nodes on the plain time scale cannot resolve. It is
  the defect reported against `rtdists::ddiffusion()` in
  [rtdists issue 28](https://github.com/rtdists/rtdists/issues/28), and it was
  here too: on the issue's own example (`boundary = 0.5`, `drift = 0.5`,
  `bias = 0.3`, `sigmandt = 0.16`, decision time 0.16) the density was out by
  5e-4, by 3% with `bias = 0.1` and `sigmandt = 0.2`, and by 50% with the start
  point almost on the responding boundary. The rule now runs over *log*
  decision time, from the point where the density is dead rather than from
  zero, so the peak is about one log unit wide wherever it sits and the same
  25 nodes resolve it at any time scale, and the rule takes more nodes only
  when the log range is wide enough to need them: the worst error over the
  issue's sweep is now below 1e-9, and below 4e-12 over a much broader grid,
  against a converged 1600-node rule. Densities with `sigmandt = 0` are
  unchanged to the last bit, and the common case - `sigmandt` well inside the
  response time, where the old rule was already accurate - gives the same
  values at the same cost. The Stan likelihood uses Stan's own adaptive
  `wiener_lpdf()` and was never affected, so fitted models are unchanged; only
  their R-side post-processing moves, and only in that regime. The
  R-versus-Stan density test is tightened from a relative 1e-4 to 1e-5
  accordingly, which is Stan's own tolerance.

* **`cogmod_rdm()` no longer produces divergent transitions by the hundred on
  healthy posteriors.** The Stan log-survival of the losing accumulator formed
  its reflection term as `log_diff_exp(log R(b), log R(k))`, and for a response
  less than about half a millisecond above the non-decision time both normal
  CDFs in `R` round to exactly 1, so it evaluated `log_diff_exp(0, 0)`. The
  *value* is fine (`-Inf` for a term that really is negligible there, which
  is why the R-versus-Stan density tests never caught it), but its
  reverse-mode adjoint is `0 / 0`, and that `NaN` propagated into the gradient
  of every parameter. Stan reports a `NaN` gradient as a divergent transition,
  and because `ndt` is estimated a few milliseconds below the fastest responses,
  most trajectories crossed one of those windows: on the lexical decision data
  of the decision-making article, 900 trials from 6 participants, between a
  third and two thirds of the transitions were divergent - with
  population-level effects only or with participant intercepts, under `diag_e`
  or `dense_e` - while `Rhat` and the effective sample sizes said the posterior
  was fine, because it was. The difference is now assembled from three pieces
  that each stay away from the saturated end of the normal CDF, at the cost of
  two extra normal CDFs on early responses only. Values are unchanged to
  `1e-12` on the log scale; the same fits now run without a divergence
  (population-level) or with the handful the other families also show under
  `dense_e` with random effects. The other race families are unaffected. A
  gradient regression test guards it, gated behind `COGMOD_TEST_SLOW` like the
  other tests that compile a model of their own.

* **`rcogmod_ddm()` no longer returns the odd 10-40 s response in place of a
  fast one.** The sampler inverts the CDF with Newton's method, and took a step
  under `1e-10` in log time as convergence. A pass that lands far in the tail
  finds the survival and the density both denormal, and their ratio makes the
  step look tiny while the residual is still hundreds of log units off, so the
  draw was accepted where it stood. It hit about one draw in 4,000 to 20,000
  at short boundaries or strong drifts (`boundary = 0.3`, `bias = 0.3`,
  `drift = -5` is one such cell; `boundary = 2`, `bias = 0.7`, `drift = 6`
  another), which is rare in `rcogmod_ddm()` and a handful of absurd draws per
  observation in `posterior_predict()`, where every observation gets thousands.
  Convergence now also requires a small residual, and a test pushes every draw
  back through the CDF to check it lands on its own quantile.

## Performance

* **`rcogmod_ddm()` is 3x faster, and `posterior_predict()` on a
  `cogmod_ddm()` model up to 16x.** The sampler inverts a series whose length
  is set by the fastest response it could be asked for, and used that length
  for every draw: 41 terms at the default start point where the median draw
  needs 6, 205 at a start point of 0.1. It now runs in stages - 16 terms
  settle the bulk of the draws, and only the responses too fast for that many
  go round again with four times as many, until the full series is reached.
  Nothing is approximated: a draw is only accepted from a stage whose series
  is exact at its root, and the draws agree with the CDF to `1e-12`. The gain
  is largest where the parameters vary across draws, as they do in
  `posterior_predict()`, because one extreme posterior draw used to set the
  series length for all of them. Converged draws now also drop out of the
  Newton iteration, and a fast response starts from the single-barrier
  small-time approximation rather than from the floor of the bracket.

* **The choice families' `posterior_predict_*()` methods take a vector of
  observations**, returning the draws stacked with those for `i[1]` first.
  `brms::posterior_predict()` calls the method once per observation, and with
  a few dozen draws per call about half of each call is fixed cost and the loop
  adds as much again; predicting in chunks of ~50 observations from a prepared
  `brmsprep` instead runs a posterior predictive check on 2,500 DDM trials in
  about a third of the time. The recipe is in `?posterior_predict_cogmod_ddm`.
  The DDM sampler's own fixed cost per call is also down by about 15%, from
  indexed assignment in place of `ifelse()` and no column copies while every
  draw is still active; the draws are bit-identical.

* **The R-side DDM density no longer goes through `brms::dwiener()`.**
  `dcogmod_ddm()`, and with it `log_lik()`, `loo()`, `p_outlier()` and every
  other post-processing method that evaluates the likelihood in R, now use a
  vectorised Navarro and Fuss (2009) series written in log space. The
  4-parameter density is about eight times cheaper per element and agrees with
  `brms::dwiener()` to `1e-12` on the log scale. The 7-parameter density, which
  evaluates that series 625 times per observation under Gauss-Legendre
  quadrature, goes from about 5 ms to about 0.6 ms per draw-observation - a
  LOO over 4000 draws of 500 trials drops from close to three hours to about
  twenty minutes. Both
  now return a finite log-density in the far tails where `dwiener()` returns
  `log(0)`. The Stan likelihood is unchanged. `RWiener` is still needed by the
  test suite, which uses `dwiener()` as the reference.

## Documentation

* The performance article is reorganised from the suggestions with no downside
  to the ones that need judgement, and gains four sections. **Compiler
  optimizations**: stanc's `O1` and CmdStan's `STAN_CPP_OPTIMS` and
  `STAN_NO_RANGE_CHECKS`, passed through `stan_model_args`, and what each one
  does. **Mass matrix adaptation** (`metric = "dense_e"`): why the
  boundary/ndt and drift/boundary trade-offs of evidence accumulation models
  make the default diagonal metric a poor fit, what the dense metric costs as
  the number of parameters grows, and how to pass it through either backend.
  **Warm starts**: reusing the adapted metric and step size that `brms` keeps
  in a fit's metadata to shorten the warmup of a refit; carrying a pilot fit's
  metric, step size and posterior means over to the same model on more
  participants, which has more parameters, by mapping the metric across by
  parameter name (about twice the effective draws per second of a cold start
  on a mixed LNR and a mixed DDM, where the pilot's initial values alone
  bought nothing); and a `cmdstanr`-level pipeline that initializes MCMC from
  Pathfinder draws and their unconstrained covariance, then wraps the result
  back into a `brmsfit`, with the reasons never to fix the metric to a
  variational approximation. The approximation
  section now also covers the **Laplace approximation**
  (`algorithm = "laplace"`) and how it compares with Pathfinder. Each section
  reports what the option bought on the DDM, LBA, LNR and RDM in a local
  benchmark; the scripts behind those numbers live in `benchmarks/` (not part
  of the installed package) and can be rerun on any model.

## Breaking changes

* **`cogmod_lba2()` now truncates each drift rate at zero**, the convention of
  `rtdists` (`posdrift = TRUE`), `DMC`, `EMC2` and `ggdmc`. Previously the pair
  of drifts was conditioned on at least one being positive and a losing
  accumulator was allowed a negative rate, which it kept forever. The two are
  different models of the same race wherever a drift is small relative to its
  SD: densities up to about 40-50% apart in the tails, choice probabilities a
  few percentage points apart. The change makes `cogmod` LBA estimates directly
  comparable with those packages and with the literature built on them, and
  `dcogmod_lba2()` now reproduces `rtdists::dLBA()` at the same parameter
  values. It also simplifies the sampler, which draws each drift from its
  truncated Normal rather than splitting the conditional law into cases.
  **Fits made with earlier versions cannot be post-processed with this one**,
  and their estimates are not comparable with new ones at low drift rates;
  the vignette model was refit. The loser's survival is computed as
  `P(v > 0, unfinished) / P(v > 0)` from whichever tail keeps its digits, so
  the density stays accurate for a loser with a strongly negative drift, where
  both quantities are tiny, and in the far tail, where `1 - CDF` would cancel.
  `cogmod_lba1()` is unaffected: with one accumulator the two conventions
  coincide.

  The truncation has one cost, and `cogmod_priors()` now covers it: once an
  accumulator rarely wins, its `drift` and `sigma` are identified only through
  `|drift| / sigma^2` (the truncated Normal converges to an Exponential along
  that ray), so a flat prior lets the drift run off - the vignette's error
  accumulator sat at `-12` with an interval of `-23` to `-6.5`. `cogmod_priors()`
  therefore puts `normal(1, 2)` on `driftone` and `normal(0, 1.5)` on its
  slopes, the treatment `cogmod_lnr()`'s `nuone` already had. Existing formulas
  that left `driftone` to `brms` get this prior on their next
  `cogmod_priors()` call.

# cogmod 0.3.0

- CRAN Publication.