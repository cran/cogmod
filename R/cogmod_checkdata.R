# Checking the data before anything compiles
# ==========================================
#
# Everything here is stated in SECONDS and none of it is unit-equivariant:
# cogmod_priors() centres `ndt` on 0.30 s, `.POUTLIER_SCALE` is 0.2 s, and
# every registry `init` is a number of seconds. Hand these families a column of
# milliseconds and nothing complains - Stan compiles, the sampler runs, and what
# comes back is a fit to a distribution a thousand times too wide. That is a
# twenty-minute way to learn about a units mistake, and it is the single most
# likely thing to be wrong with a data frame arriving here.
#
# The checks live in cogmod_priors() because it is the one of the three generics
# that sees the data AND is documented as near-mandatory for these families - a
# flat prior on `ndt` and `poutlier` is an improper posterior, so a call that
# skips it is already broken for a different reason.
#
# Fatal versus advisory is decided by one question: does the offending row send
# the total log-likelihood to -Inf, or silently select the wrong branch of the
# density? Stan cannot recover from either, so those stop(). Everything else
# warns, because a script that has run before should not start failing over a
# judgement call about what a plausible reaction time is.

# What the response column has to look like, per family. The two RT groups come
# off the registries so a new entry there is covered without touching this; the
# other two are listed literally because there is no registry behind them.
#' @keywords internal
.checkdata_class <- function(fam) {
  if (is.null(fam)) {
    return(NULL)
  }
  if (fam %in% .CHOICE_FAMILIES) {
    return("choice")
  }
  if (fam %in% c(names(.SHIFTED), "cogmod_exgaussian", "cogmod_geg")) {
    return("rt")
  }
  if (fam %in% c("cogmod_choco", "cogmod_betagate")) {
    return("unit")
  }
  if (identical(fam, "cogmod_betadiscrete")) {
    return("rating")
  }
  NULL
}

# Reduce a dec() column to the 0/1 integers the Stan code branches on, or NULL
# if it cannot be done. Mirrors what brms accepts for `dec`: 0/1, a logical, or
# a two-level factor/character whose SECOND level is the upper boundary.
#' @keywords internal
.checkdata_dec <- function(v) {
  if (is.logical(v)) {
    return(as.integer(v))
  }
  if (is.factor(v) || is.character(v)) {
    lv <- if (is.factor(v)) levels(v) else sort(unique(stats::na.omit(v)))
    if (length(lv) != 2L) {
      return(NULL)
    }
    return(match(as.character(v), lv) - 1L)
  }
  if (is.numeric(v)) {
    return(v)
  }
  NULL
}

#' @keywords internal
.cogmod_checkdata <- function(formula, data, family = NULL) {
  if (is.null(family)) {
    family <- .cogmod_family(formula)
  }
  cls <- .checkdata_class(.family_name(family))
  if (is.null(cls) || !is.data.frame(data) || !nrow(data)) {
    return(invisible(NULL))
  }

  # A formula this cannot read is not an error here: cogmod_priors() is about to
  # hand the same object to brms, which will say something far more useful about
  # it than anything guessed at from this end.
  bterms <- tryCatch(brms::brmsterms(formula), error = function(e) NULL)
  if (is.null(bterms)) {
    return(invisible(NULL))
  }
  resp <- all.vars(bterms$respform)
  if (length(resp) != 1L || !resp %in% names(data)) {
    return(invisible(NULL))
  }

  y <- data[[resp]]
  if (!is.numeric(y) && !is.logical(y)) {
    stop(
      "The response `",
      resp,
      "` is ",
      class(y)[1],
      ", and every cogmod ",
      "family needs a numeric one. ",
      if (cls %in% c("rt", "choice")) {
        "Reaction times go in as seconds, not as a factor or a string."
      } else {
        "Convert it before fitting."
      },
      call. = FALSE
    )
  }
  y <- as.numeric(y)

  n_na <- sum(is.na(y))
  if (n_na) {
    warning(
      n_na,
      " of ",
      length(y),
      " values of `",
      resp,
      "` are NA - brms ",
      "will drop those rows.",
      call. = FALSE
    )
  }
  y <- y[!is.na(y)]
  if (!length(y)) {
    stop("Every value of the response `", resp, "` is NA.", call. = FALSE)
  }

  switch(
    cls,
    rt = ,
    choice = .checkdata_rt(y, resp, family),
    unit = .checkdata_unit(y, resp),
    rating = .checkdata_rating(y, resp)
  )
  if (identical(cls, "choice")) {
    .checkdata_choice(bterms, data, family)
  }
  # Runs after the per-class response checks, so that a response off the
  # support is reported as that rather than as a pinned-weight problem.
  .checkdata_pinned(formula, y, resp, family)
  .checkdata_cens(bterms, data, family)
  invisible(NULL)
}

# Is there a `cens()` addition term on the formula? Tolerant of anything that is
# not a formula, because cogmod_stanvars() also accepts a family object or a
# fitted model, and neither carries one.
#' @keywords internal
.has_cens <- function(formula) {
  if (!inherits(formula, "brmsformula")) {
    return(FALSE)
  }
  bt <- tryCatch(brms::brmsterms(formula), error = function(e) NULL)
  !is.null(bt) && !is.null(bt$adforms$cens)
}

# Refuse `cens()` on a family that cannot take it, from wherever the formula is
# first seen. Left alone, brms generates the `<family>_lccdf()` call regardless
# and the user meets a Stan compilation error about an undefined function some
# minutes later.
#' @keywords internal
.check_cens_family <- function(fam) {
  if (isTRUE(fam %in% .CENS_FAMILIES)) {
    return(invisible(NULL))
  }
  if (isTRUE(fam %in% .CHOICE_FAMILIES)) {
    # A choice model's likelihood is a set of defective densities summing to one
    # over the options, and the contaminant's 1 / K exists to keep it so. A
    # censored likelihood - a density for one outcome, a bare survival for the
    # other - breaks that identity by construction, so it does not belong here
    # even where the Stan code could be written.
    stop(
      fam,
      "() already models the errors, through dec(): each response ",
      "option has its own accumulator and the two are scored jointly. ",
      "`cens()` belongs on the RT-only families - bf(rt | cens(error) ~ ..., ",
      "family = cogmod_invgaussian()) - where an error is treated as a ",
      "right-censored correct response instead. See ?rcogmod_invgaussian.",
      call. = FALSE
    )
  }
  stop(
    fam,
    "() has no closed-form CDF, so it cannot take `cens()`. Families ",
    "that can: ",
    paste0(.CENS_FAMILIES, "()", collapse = ", "),
    ".",
    call. = FALSE
  )
}

# The cens() indicator. brms accepts 0/1, a logical, -1/0/1/2, or the strings
# "left"/"none"/"right"/"interval"; this mirrors that reading so the proportion
# can be judged before anything compiles.
#' @keywords internal
.checkdata_cens <- function(bterms, data, family) {
  censform <- bterms$adforms$cens
  if (is.null(censform)) {
    return(invisible(NULL))
  }
  fam <- .family_name(family)
  .check_cens_family(fam)

  cvar <- all.vars(censform)[1]
  if (is.na(cvar) || !cvar %in% names(data)) {
    return(invisible(NULL))
  }
  v <- data[[cvar]]
  code <- if (is.logical(v)) {
    as.integer(v)
  } else if (is.factor(v) || is.character(v)) {
    key <- c(left = -1L, none = 0L, right = 1L, interval = 2L)
    unname(key[tolower(as.character(v))])
  } else if (is.numeric(v)) {
    v
  } else {
    NULL
  }
  code <- code[!is.na(code)]
  bad <- if (is.null(code)) TRUE else !all(code %in% c(-1, 0, 1, 2))
  if (bad) {
    stop(
      "`",
      cvar,
      "` cannot be read as a censoring indicator. brms takes ",
      "0/1 or a logical (1 = right-censored), -1/0/1/2, or the strings ",
      "\"left\", \"none\", \"right\" and \"interval\".",
      call. = FALSE
    )
  }
  if (!length(code)) {
    return(invisible(NULL))
  }

  # The model these families fit under cens() says a censored trial tells you
  # ONE thing about the correct process: that it had not finished yet. For a
  # timeout or an omission that is exactly true at any rate. For a commission
  # error it is true only if the error came from a process independent of the
  # correct one; under a diffusion process the two responses share one evidence
  # path and the error RTs are as fast as the correct ones, and the survival
  # scores them wrong. In simulation that bias is already ~15% on the drift at
  # 5% errors, which is why Miller et al. (2018) argue for this "simple"
  # censored Wald only above ~95% accuracy. The threshold stays at 20% - bmm's
  # for the same model - because the indicator cannot tell an omission from an
  # error, and warning at 5% would fire on every deadline design; the message
  # carries the stricter figure instead.
  share <- mean(code == 1)
  if (share > 0.2) {
    warning(
      round(100 * share, 1),
      "% of the trials are right-censored. Under ",
      "`cens()` a censored trial says only that the decision process ",
      "had not finished by then. In practice, Miller et al. (2018) argue ",
      "for the potential suitability of this model for less than 5% of errors, ",
      "as with more it can severely bias the drift. Consider a race that models ",
      "errors as such, and compare the ",
      "censored trials' RTs to the others': censoring can only ever ",
      "produce them slower.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Reaction times, in seconds.
#' @keywords internal
.checkdata_rt <- function(y, resp, family) {
  fam <- .family_name(family)
  # The families on the ndt + poutlier mixture put ZERO density below `ndt`, and
  # `ndt` is itself bounded below by zero, so a non-positive reaction time sends
  # the whole log-likelihood to -Inf and no chain can initialise. The
  # ex-Gaussian and the gamma-ex-Gaussian have support on the entire real line,
  # so there the same row is implausible rather than fatal.
  shifted <- isTRUE(fam %in% .OUTLIER_FAMILIES)
  bad <- sum(y <= 0)
  if (bad) {
    msg <- paste0(
      bad,
      " of ",
      length(y),
      " values of `",
      resp,
      "` are zero or negative."
    )
    if (shifted) {
      stop(
        msg,
        " ",
        fam,
        "() places no density below `ndt`, which is itself ",
        "positive, so these rows make the log-likelihood -Inf and the ",
        "sampler cannot start. Drop or correct them.",
        call. = FALSE
      )
    }
    warning(
      msg,
      " These are reaction times in seconds, so that is almost ",
      "certainly a coding or a merge error.",
      call. = FALSE
    )
  }

  # Milliseconds are diagnosed from the MEDIAN, not from the maximum. Testing
  # `any(rt > 10)` - the obvious thing to write - fires on any data set holding
  # a single slow trial, and would cry wolf on most real ones. A units mistake
  # moves the entire column by three orders of magnitude, so the median is what
  # carries the signal: no task in scope here has a median response anywhere
  # near a minute.
  pos <- y[y > 0]
  if (!length(pos)) {
    return(invisible(NULL))
  }
  med <- stats::median(pos)
  if (med > 50) {
    warning(
      "The median of `",
      resp,
      "` is ",
      signif(med, 4),
      ". cogmod states everything in SECONDS - the `ndt` prior, the ",
      "outlier scale, every starting value - so this looks like ",
      "milliseconds. Divide by 1000 and refit; left as it is the model ",
      "will fit, and the estimates will be meaningless.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  # The two tails are only worth a word for the families carrying `poutlier`,
  # because it is `poutlier` that decides whether an implausible response is a
  # problem or the thing the model was built to absorb. The ex-Gaussian and the
  # gamma-ex-Gaussian have neither `ndt` to be pushed down nor `poutlier` to do
  # the absorbing, so neither message would be true of them.
  if (!shifted) {
    return(invisible(NULL))
  }

  # A COUNT is the wrong test here, and the package's own generator is the
  # proof: rcogmod_lognormal(200, ndt = 0.2, poutlier = 0.02) puts a response at
  # 81 ms, because that is exactly what the outlier component is for. What
  # matters is the PROPORTION, measured against how much of it `poutlier` can
  # plausibly claim. Over 20000 draws the outlier component sends 0.8% of
  # responses below 0.1 s at poutlier = 0.02, and 1.9% at poutlier = 0.05 - the
  # top of the default normal(-5, 1) prior. Past 5% the mixture cannot account
  # for them and `ndt` has to come down under every one instead, so that is
  # where the warning sits.
  # How far past the plausible range the outlier component can be asked to go.
  reach <- 0.05
  fast <- mean(pos < 0.1)
  if (fast > reach) {
    warning(
      round(100 * fast, 1),
      "% of `",
      resp,
      "` is under 0.1 s. The ",
      "outlier component reaches about 2% there at the top of its ",
      "default prior, so it cannot account for this many and `ndt` will ",
      "be pushed down under all of them instead. Filter them, or widen ",
      "the `poutlier` prior if they are genuine fast guesses.",
      call. = FALSE
    )
  }
  # The same test on the other side, where the outlier component is much weaker
  # still: its scale is .POUTLIER_SCALE = 0.2 s, so it sends essentially nothing
  # past 10 s and every one of these has to be carried by the decision
  # distribution's own tail.
  slow <- mean(pos > 10)
  if (slow > reach) {
    warning(
      round(100 * slow, 1),
      "% of `",
      resp,
      "` exceeds 10 s. The outlier ",
      "component has a scale of ",
      .POUTLIER_SCALE,
      " s and reaches ",
      "almost none of that, so the decision distribution's own tail has ",
      "to stretch to cover it. Check whether these trials belong in the ",
      "data.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# The dec() index the choice families branch on.
#' @keywords internal
.checkdata_choice <- function(bterms, data, family) {
  decform <- bterms$adforms$dec
  if (is.null(decform)) {
    stop(
      .family_name(family),
      "() needs the chosen option, which goes in the ",
      "formula as an addition term: bf(rt | dec(response) ~ ...).",
      call. = FALSE
    )
  }
  dvar <- all.vars(decform)
  if (length(dvar) != 1L || !dvar %in% names(data)) {
    return(invisible(NULL))
  }

  d <- .checkdata_dec(data[[dvar]])
  if (is.null(d)) {
    stop(
      "`",
      dvar,
      "` cannot be read as a two-option choice. Code it as 0/1, ",
      "as a logical, or as a factor with exactly two levels (the second ",
      "being the upper boundary).",
      call. = FALSE
    )
  }
  d <- d[!is.na(d)]
  bad <- setdiff(unique(d), c(0, 1))
  if (length(bad)) {
    # An error rather than a warning, because the Stan code tests `dec == 0` and
    # takes the else branch for everything else. A third response level would
    # not fail - it would be silently folded into option 1 - and every one of
    # these families has K = 2, so there is no reading of the data under which
    # that is what was meant.
    stop(
      "`",
      dvar,
      "` takes the value",
      if (length(bad) > 1) "s " else " ",
      paste(sort(bad)[seq_len(min(3L, length(bad)))], collapse = ", "),
      " as well as 0 and 1. The choice families are two-option: Stan tests ",
      "`dec == 0` and treats everything else as option 1, so a third level ",
      "is silently absorbed rather than rejected.",
      call. = FALSE
    )
  }
  if (length(unique(d)) < 2L) {
    warning(
      "`",
      dvar,
      "` only ever takes the value ",
      unique(d),
      ". With one option never chosen, the parameters governing the ",
      "other accumulator are informed by nothing but the prior.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Bounded ratings on [0, 1].
#' @keywords internal
.checkdata_unit <- function(y, resp) {
  bad <- sum(y < 0 | y > 1)
  if (bad) {
    rng <- range(y)
    stop(
      bad,
      " of ",
      length(y),
      " values of `",
      resp,
      "` fall outside [0, 1] ",
      "(the observed range is ",
      signif(rng[1], 4),
      " to ",
      signif(rng[2], 4),
      "). These families have no density there, so the log-likelihood is ",
      "-Inf. Rescale to the unit interval first - exact 0s and 1s are fine, ",
      "and are what the zero-one components are for.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Discrete ratings on 0, 1, ..., k. The upper bound `k` arrives in the formula
# as vint() and is already checked by a reject() in the Stan code; what is worth
# catching from this end is the shape of the column.
#' @keywords internal
.checkdata_rating <- function(y, resp) {
  if (any(y != trunc(y))) {
    stop(
      "`",
      resp,
      "` has non-integer values, and cogmod_betadiscrete() is a ",
      "probability mass function over 0, 1, ..., k.",
      call. = FALSE
    )
  }
  if (any(y < 0)) {
    stop(
      "`",
      resp,
      "` has negative values. cogmod_betadiscrete() runs over ",
      "0, 1, ..., k, where 0 is the separate zero-response category and the ",
      "ratings themselves start at 1.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# Mixture weights switched off by a fixed value
# =============================================
#
# The three bounded-scale families are mixtures of a continuous part and one or
# more point masses, and every weight involved can be pinned in bf(). Pin one at
# a boundary and the component it controls has probability zero - so every
# observation belonging to that component has density zero, the log-likelihood
# is -Inf at every point in the parameter space, and CmdStan gives up with
# "Initialization failed after 100 attempts" and nothing about which column is
# to blame. It is a data-and-formula problem rather than a data one, which is
# why it is checked here rather than in .checkdata_unit() / .checkdata_rating().
#
# Both halves of the mistake are easy to make. `?rcogmod_betadiscrete` offers
# `pzero = 0` as the way to say "my scale has no zero category", and
# ?cogmod_priors recommends fixing `pmid` or `pzero` at 0 to switch the
# parameter off outright - neither of which is safe if the column still has one
# of those responses in it.
#
# One entry per (dpar, pinned value) that closes a component, with the responses
# that then become impossible. The gate arithmetic behind the `pex` and `bex`
# rows is `cutzero = pex * (1 - bex)` and `cutone = 1 - pex * bex`: `pex = 0`
# opens both gates fully so neither 0 nor 1 can be produced, and `bex` at either
# end closes one of them. cogmod_choco()'s `mid` is hard-coded at 0.5 in the
# Stan code, so the midpoint is a constant here too. Note that `pex = 1` is
# fatal for cogmod_betagate() but not for cogmod_choco(): it collapses
# betagate's two gates onto each other, leaving no room for the continuous
# part, whereas each of choco's halves uses only one gate.
#' @keywords internal
.CHECKDATA_PINNED <- list(
  cogmod_choco = list(
    list(dpar = "pmid", at = 0, hits = function(y) y == 0.5,
         what = "exactly at the midpoint 0.5"),
    list(dpar = "pmid", at = 1, hits = function(y) y != 0.5,
         what = "anywhere other than the midpoint 0.5"),
    list(dpar = "pex", at = 0, hits = function(y) y == 0 | y == 1,
         what = "exactly 0 or exactly 1"),
    list(dpar = "bex", at = 0, hits = function(y) y == 1,
         what = "exactly 1"),
    list(dpar = "bex", at = 1, hits = function(y) y == 0,
         what = "exactly 0")
  ),
  cogmod_betagate = list(
    list(dpar = "pex", at = 0, hits = function(y) y == 0 | y == 1,
         what = "exactly 0 or exactly 1"),
    list(dpar = "pex", at = 1, hits = function(y) y > 0 & y < 1,
         what = "strictly between 0 and 1"),
    list(dpar = "bex", at = 0, hits = function(y) y == 1,
         what = "exactly 1"),
    list(dpar = "bex", at = 1, hits = function(y) y == 0,
         what = "exactly 0")
  ),
  cogmod_betadiscrete = list(
    list(dpar = "pzero", at = 0, hits = function(y) y == 0,
         what = "exactly 0"),
    list(dpar = "pzero", at = 1, hits = function(y) y >= 1,
         what = "a rating of 1 or more")
  )
)


# The dpars a formula pins to a single number, as a named numeric vector.
# Anything that is not one plain number - an expression, a name, a vector - is
# dropped rather than guessed at, the same tolerance .warn_scale_ray() applies.
#' @keywords internal
.checkdata_pfix <- function(formula) {
  if (!inherits(formula, "brmsformula")) {
    return(stats::setNames(numeric(0), character(0)))
  }
  pfix <- formula$pfix
  if (!length(pfix)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  num <- vapply(
    pfix,
    function(v) {
      out <- tryCatch(suppressWarnings(as.numeric(v)), error = function(e) NA_real_)
      if (length(out) != 1L || !is.finite(out)) NA_real_ else out
    },
    numeric(1)
  )
  num[!is.na(num)]
}


#' @keywords internal
.checkdata_pinned <- function(formula, y, resp, family) {
  fam <- .family_name(family)
  if (is.null(fam) || !fam %in% names(.CHECKDATA_PINNED)) {
    return(invisible(NULL))
  }
  pfix <- .checkdata_pfix(formula)
  if (!length(pfix)) {
    return(invisible(NULL))
  }

  for (e in .CHECKDATA_PINNED[[fam]]) {
    if (!e$dpar %in% names(pfix) || pfix[[e$dpar]] != e$at) {
      next
    }
    bad <- sum(e$hits(y))
    if (!bad) {
      next
    }
    stop(
      bad, " of ", length(y), " values of `", resp, "` are ", e$what,
      ", and `", e$dpar, " = ", e$at, "` in bf() gives those zero ",
      "probability. The log-likelihood is -Inf everywhere, so the fit cannot ",
      "start - CmdStan reports it as `Initialization failed after 100 ",
      "attempts` without saying why. Either drop those rows, or let `",
      e$dpar, "` be estimated: leave it out of bf() to get it as a plain ",
      "auxiliary parameter, or write `", e$dpar, " ~ 1` to model it.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
