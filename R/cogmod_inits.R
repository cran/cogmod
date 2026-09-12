#' Starting values that keep the sampler out of the flat regions
#'
#' @description
#' Builds an `init` argument for `brms::brm()`, for the custom families whose
#' default starting point is a bad one.
#'
#' **This is not a tuning knob to reach for when sampling looks bad.** For
#' [cogmod_gamma()] and [cogmod_weibull()] the usual default is actively harmful, and
#' the failure is silent: the chain does not error, it simply never moves.
#'
#' @details
#' Two regions of the parameter space have no gradient to escape on, and both
#' are easy to start inside.
#'
#' `ndt` **too large.** `brms` initialises on the *unconstrained* scale, so
#' `init = 0` puts `ndt` at `exp(0) = 1` second. For sub-second reaction times
#' that is above nearly every observation, so every response is attributed to
#' the outlier component, the decision parameters drop out of the density
#' entirely, and their gradient is exactly zero. The chain is stuck where it
#' started.
#'
#' **Shape below 1.** For [cogmod_gamma()] and [cogmod_weibull()] the density is
#' unbounded at `ndt` whenever the shape is below 1, and `init = 0` on a
#' `softplus` link starts the shape at `log(2) = 0.69` - inside that region.
#'
#' A single scalar `init` cannot avoid both, because they pull in opposite
#' directions: `ndt = exp(c)` wants `c` around `-1.6`, while
#' `shape = softplus(c)` wants `c` above `1.9`. That is why `init = 0` fails and
#' why no other single number fixes it - the parameters have to be set
#' separately, which is what this function does.
#'
#' # How it works
#'
#' The Stan parameter declarations are read off `brms::make_stancode()` for the
#' model you are actually fitting, and their dimensions off
#' `brms::make_standata()`, rather than reconstructed from the formula. That is
#' what makes it robust to `0 + Intercept`, interactions, group-level terms and
#' smooths: whatever `brms` decided to call things, and however large it decided
#' to make them, that is what is matched.
#'
#' **Every** declared parameter is given a value, not only the ones this
#' function has an opinion about. That is deliberate: CmdStan prints
#' `Init values were only set for a subset of parameters` and lists the rest
#' whenever the list is partial, which is noise for a list that is partial on
#' purpose. The parameters with no family-specific target - regression slopes,
#' standardized group-level effects, group-level SDs, spline coefficients - get
#' generic values that are at least as good as Stan's own `U(-2, 2)`: slopes and
#' standardized effects start at zero, positive parameters just above their
#' lower bound, bounded ones at the midpoint, and Cholesky factors at the
#' identity.
#'
#' The family-specific values go to the intercept of each distributional
#' parameter, and to any dpar left out of the formula (which `brms` declares as
#' a plain auxiliary parameter). Values are set on whichever scale the parameter
#' lives on: the link scale for an intercept - using the links on the family, so
#' `cogmod_gamma(link_mu = "log")` is honoured - and the natural scale for an
#' auxiliary parameter. Under `0 + Intercept` there is no separate intercept
#' parameter, so the coefficient named `Intercept` inside the `b` vector gets it
#' instead.
#'
#' Each chain gets a different draw. The noise is added on the *unconstrained*
#' scale - additive for a free parameter, multiplicative for a positive one,
#' on the logit scale for a doubly bounded one - so a jittered value can never
#' land outside its own bounds, and the chains still start dispersed enough for
#' `Rhat` to mean something.
#'
#' `ndt` starts deliberately **small** (0.1 s, a third of its prior median).
#' The two errors are
#' not symmetric: too small merely means the shift has to grow, which the
#' gradient will do, whereas too large removes the gradient altogether.
#'
#' # Supported families
#'
#' The family is read off `formula`, so build it with
#' `brms::bf(..., family = cogmod_gamma())`. Every family built on the direct
#' `ndt` + `poutlier` parameterization is covered - [cogmod_lognormal()],
#' [cogmod_logstudent()], [cogmod_loggamma()], [cogmod_invgaussian()], [cogmod_exwald()],
#' [cogmod_bisa()], [cogmod_gamma()],
#' [cogmod_invgamma()],
#' [cogmod_weibull()], [cogmod_invweibull()], [cogmod_logweibull()], [cogmod_lba1()] and,
#' for the choice-and-RT models, [cogmod_lnr()], [cogmod_rdm()],
#' [cogmod_lba2()] and [cogmod_ddm()] - plus
#' [cogmod_exgaussian()], whose three parameters are all on the RT scale behind a
#' `softplus` link and so are equally badly served by starting at `log(2)`.
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = cogmod_gamma())`.
#' @param data The data, as passed to `brms::brm()`.
#' @param jitter SD of the noise added on the unconstrained scale, so that
#'   chains start at different points. Set to `0` for identical starts.
#' @param ... Passed to `brms::make_stancode()` and `brms::make_standata()`, for
#'   arguments such as `data2`.
#'
#' @return A function of one argument, suitable for `brms::brm(init = )`. Each
#'   call returns a named list of starting values, one for every parameter the
#'   generated Stan program declares.
#'
#' @seealso [cogmod_priors()], [cogmod_stanvars()]
#'
#' @examples
#' d <- data.frame(RT = rcogmod_gamma(50, ndt = 0.3, poutlier = 0.02))
#' f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_gamma())
#' inits <- cogmod_inits(f, d)
#' inits(1)
#'
#' \donttest{
#' # Fitting needs cmdstanr, which lives outside CRAN - see the package website.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   m <- brms::brm(f,
#'     data = d, prior = cogmod_priors(f, d),
#'     stanvars = cogmod_stanvars(f), init = cogmod_inits(f, d),
#'     backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
#'   )
#' }
#' }
#'
#' @export
cogmod_inits <- function(formula, data, jitter = 0.25, ...) {
  family <- .cogmod_family(formula)
  fam <- .family_name(family)
  targets <- if (is.null(fam)) NULL else .init_targets(family)
  if (is.null(targets)) {
    stop(
      "cogmod_inits() has nothing to offer for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'. Supported: ", paste(.init_families(), collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_gamma()).",
      call. = FALSE
    )
  }
  links <- .family_links(family)

  # Whatever brms decided to call the parameters, and however large it decided
  # to make them, that is what we work from. Warnings are suppressed because
  # only the structure is wanted: the code is generated with brms' own default
  # priors, which is exactly the situation cogmod_priors() exists to fix, and
  # complaining about it here would be noise.
  code <- suppressWarnings(
    brms::make_stancode(formula, data = data, family = family, ...)
  )
  sdata <- suppressWarnings(
    brms::make_standata(formula, data = data, family = family, ...)
  )
  plan <- .init_plan(.stan_param_decls(code), as.list(sdata), targets, links)

  function(chain_id = 1) {
    out <- lapply(plan, function(e) {
      v <- e$value
      if (jitter > 0) {
        v <- switch(
          e$kind,
          bounds = .jitter_bounded(v, e$lower, e$upper, jitter),
          sorted = sort(v + stats::rnorm(length(v), 0, jitter)),
          v # "fixed": a structured value that jitter would invalidate
        )
      }
      if (length(e$dim) > 1) dim(v) <- e$dim
      v
    })
    stats::setNames(out, vapply(plan, `[[`, character(1), "name"))
  }
}


# Natural-scale starting values ------------------------------------------

# Which families cogmod_inits() has targets for.
#' @keywords internal
.init_families <- function() c(.OUTLIER_FAMILIES, "cogmod_exgaussian")

# The natural-scale value each dpar should start at, or NULL for a family with
# no opinion attached. The shifted families take theirs from the registry in
# shifted.R and add the two the parameterization introduces; anything else
# is listed here.
#' @keywords internal
.init_targets <- function(family) {
  fam <- .family_name(family)
  if (is.null(fam)) return(NULL)
  if (fam %in% .OUTLIER_FAMILIES) {
    if (is.character(family)) {
      stop(
        "cogmod_inits() needs the family object rather than its name, ",
        "because the links ride on it. Build the formula with ",
        "bf(..., family = ", fam, "()).",
        call. = FALSE
      )
    }
    # In seconds, like the `ndt` prior: a third of its 0.30 s median.
    return(c(.mixture_spec(fam)$init, list(ndt = 0.1, poutlier = 0.02)))
  }
  switch(
    fam,
    # mu and tau are the Gaussian centre and the exponential mean, both in
    # seconds; sigma is the Gaussian SD. Behind the softplus link that `sigma`
    # and `tau` use, the default start puts both at log(2) = 0.69 s, which makes
    # sigma alone wider than most whole RT distributions. `mu` is on identity,
    # where the generic start of 0 would put the Gaussian centre at zero
    # seconds - harmless but pointless when 0.4 s is known to be the right
    # neighbourhood.
    cogmod_exgaussian = list(mu = 0.4, sigma = 0.1, tau = 0.2),
    # cogmod_geg() starts from the ex-Gaussian it nests: shape = 1 is that
    # model exactly, and it is the one point on the shape axis known not to be
    # somewhere up the shape/mu ridge.
    cogmod_geg = list(mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1),
    NULL
  )
}


# Turning declarations into values ----------------------------------------

# One entry per declared parameter: the value on its declared scale, its shape,
# its bounds, and how it may be jittered.
#
# A declaration whose dimensions cannot be resolved is dropped rather than
# guessed at. That costs only the CmdStan message about a partial init list,
# whereas a wrong dimension is a hard error at startup.
#' @keywords internal
.init_plan <- function(decls, sdata, targets, links) {
  plan <- lapply(decls, function(d) {
    dims <- vapply(d$dims, .eval_stan_expr, numeric(1), sdata = sdata)
    if (anyNA(dims)) return(NULL)
    dims <- as.integer(dims)
    n <- if (length(dims)) prod(dims) else 1L
    bounds <- .stan_bounds(d$bounds, sdata)
    lower <- bounds[["lower"]]
    upper <- bounds[["upper"]]

    entry <- list(name = d$name, dim = dims, lower = lower, upper = upper,
                  kind = "bounds")

    # Structured types: a valid value is not just a number in a range, so the
    # bounds-based default and the jitter are both skipped.
    structured <- switch(
      d$type,
      cholesky_factor_corr = , cholesky_factor_cov = , corr_matrix = ,
      cov_matrix = as.vector(diag(dims[1])),
      simplex = rep(1 / dims[1], dims[1]),
      unit_vector = c(1, rep(0, dims[1] - 1)),
      NULL
    )
    if (!is.null(structured)) {
      entry$value <- structured
      entry$kind <- "fixed"
      if (d$type %in% c("cholesky_factor_corr", "cholesky_factor_cov",
                        "corr_matrix", "cov_matrix")) {
        entry$dim <- c(dims[1], dims[1])
      }
      return(entry)
    }
    # An ordered vector has to stay ordered, so it is jittered and re-sorted.
    if (d$type %in% c("ordered", "positive_ordered")) {
      base <- if (identical(d$type, "ordered")) 0 else 0.5
      entry$value <- base + seq_len(n) * 0.5
      entry$kind <- "sorted"
      return(entry)
    }

    entry$value <- .init_value(d, n, lower, upper, targets, links, sdata)
    entry
  })
  plan[!vapply(plan, is.null, logical(1))]
}


# The value for a plain real / vector / matrix parameter, preferring a
# family-specific target where one applies.
#' @keywords internal
.init_value <- function(d, n, lower, upper, targets, links, sdata) {
  scalar <- identical(d$type, "real") && !length(d$dims)

  # brms calls the main parameter's intercept `Intercept`, and every other
  # dpar's `Intercept_<dpar>`. A dpar omitted from the formula becomes a bare
  # auxiliary parameter named after itself, on the natural scale. Only a scalar
  # `real` is ever a dpar: brms also declares vectors named after its own
  # machinery (`vector[Ks] bs` for a smooth, for one), and matching those by
  # name would hand a spline coefficient a dpar's starting value.
  if (scalar) {
    dpar <- if (identical(d$name, "Intercept")) {
      "mu"
    } else if (startsWith(d$name, "Intercept_")) {
      sub("^Intercept_", "", d$name)
    } else {
      NA_character_
    }
    if (!is.na(dpar) && dpar %in% names(targets)) {
      v <- .apply_link(targets[[dpar]], links[dpar])
      if (!is.na(v)) return(v)
    }
    if (d$name %in% names(targets)) {
      return(targets[[d$name]]) # natural scale, no link applied
    }
  }

  # Under `0 + Intercept` brms declares no separate intercept: the column named
  # `Intercept` inside the (uncentered) design matrix carries it, so the target
  # goes to that element of `b`. The centered parameterization drops the first
  # column instead, which is how the two are told apart - there `b` is one
  # shorter than the design matrix and every element really is a slope.
  if (identical(d$type, "vector") && length(d$dims) == 1 &&
      (identical(d$name, "b") || startsWith(d$name, "b_"))) {
    dpar <- if (identical(d$name, "b")) "mu" else sub("^b_", "", d$name)
    out <- rep(0, n)
    cn <- colnames(sdata[[if (identical(dpar, "mu")) "X" else paste0("X_", dpar)]])
    if (length(cn) == n && dpar %in% names(targets) && "Intercept" %in% cn) {
      v <- .apply_link(targets[[dpar]], links[dpar])
      if (!is.na(v)) out[match("Intercept", cn)] <- v
    }
    return(out)
  }

  rep(.default_value(lower, upper), n)
}


# Where no target applies: zero for a free parameter (a slope of zero and a
# standardized group-level effect of zero are both sensible), just inside the
# bound for a constrained one, the midpoint when bounded on both sides.
#' @keywords internal
.default_value <- function(lower, upper) {
  if (is.na(lower) && is.na(upper)) return(0)
  if (is.na(upper)) return(lower + 0.25)
  if (is.na(lower)) return(upper - 0.25)
  lower + (upper - lower) / 2
}


# Jitter on the scale Stan samples on, so the result cannot leave the bounds:
# additive when free, multiplicative on the distance to a single bound, and on
# the logit scale between two.
#' @keywords internal
.jitter_bounded <- function(v, lower, upper, jitter) {
  e <- stats::rnorm(length(v), 0, jitter)
  if (is.na(lower) && is.na(upper)) return(v + e)
  if (is.na(upper)) return(lower + pmax(v - lower, 1e-8) * exp(e))
  if (is.na(lower)) return(upper - pmax(upper - v, 1e-8) * exp(e))
  w <- pmin(pmax((v - lower) / (upper - lower), 1e-8), 1 - 1e-8)
  lower + (upper - lower) * stats::plogis(stats::qlogis(w) + e)
}


# Reading the Stan program ------------------------------------------------

# Every declaration in the `parameters` block, as name / type / dimension
# expressions / bounds expression. Dimensions and bounds are left as text
# because they are Stan expressions over the data block (`Kc`, `knots_1[1]`),
# and are resolved against make_standata() output later.
#' @keywords internal
.stan_param_decls <- function(code) {
  lines <- strsplit(code, "\n")[[1]]
  start <- grep("^\\s*parameters\\s*\\{", lines)
  if (!length(start)) return(list())
  close <- grep("^\\s*\\}\\s*$", lines)
  close <- close[close > start[1]]
  if (!length(close)) return(list())
  block <- sub("//.*$", "", lines[(start[1] + 1):(close[1] - 1)])

  stmts <- trimws(strsplit(paste(block, collapse = " "), ";", fixed = TRUE)[[1]])
  out <- lapply(stmts[nzchar(stmts)], .parse_stan_decl)
  out[!vapply(out, is.null, logical(1))]
}


#' @keywords internal
.STAN_PARAM_TYPES <- c(
  "real", "vector", "row_vector", "matrix", "simplex", "ordered",
  "positive_ordered", "unit_vector", "cholesky_factor_corr",
  "cholesky_factor_cov", "corr_matrix", "cov_matrix"
)


# One declaration. Hand-rolled rather than a single regex because dimensions
# nest (`vector[knots_1[1]] zs_1_1`), which a regex cannot bracket-match.
# Anything unrecognised returns NULL and is skipped.
#' @keywords internal
.parse_stan_decl <- function(s) {
  s <- trimws(s)
  dims <- character(0)

  if (grepl("^array\\s*\\[", s)) {
    g <- .take_bracket(sub("^array\\s*", "", s))
    if (is.null(g)) return(NULL)
    dims <- .split_top_commas(g$inside)
    s <- trimws(g$rest)
  }

  type <- regmatches(s, regexpr("^[A-Za-z_][A-Za-z0-9_]*", s))
  if (!length(type) || !type %in% .STAN_PARAM_TYPES) return(NULL)
  s <- trimws(substring(s, nchar(type) + 1L))

  bounds <- ""
  if (startsWith(s, "<")) {
    i <- regexpr(">", s, fixed = TRUE)
    if (i < 0) return(NULL)
    bounds <- substring(s, 2L, i - 1L)
    s <- trimws(substring(s, i + 1L))
  }

  if (startsWith(s, "[")) {
    g <- .take_bracket(s)
    if (is.null(g)) return(NULL)
    dims <- c(dims, .split_top_commas(g$inside))
    s <- trimws(g$rest)
  }

  if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", s)) return(NULL)
  list(name = s, type = type, dims = dims, bounds = bounds)
}


# Contents of the bracketed group `s` opens with, plus whatever follows it.
#' @keywords internal
.take_bracket <- function(s) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(chars) || chars[1] != "[") return(NULL)
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[i] == "[") depth <- depth + 1L
    if (chars[i] == "]") {
      depth <- depth - 1L
      if (depth == 0L) {
        return(list(
          inside = paste(chars[seq_len(i - 1L)[-1]], collapse = ""),
          rest = paste(chars[-seq_len(i)], collapse = "")
        ))
      }
    }
  }
  NULL
}


# Split on the commas that separate arguments, ignoring those inside brackets
# or parentheses.
#' @keywords internal
.split_top_commas <- function(s) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(chars)) return(character(0))
  depth <- 0L
  grp <- 1L
  which_grp <- integer(length(chars))
  for (i in seq_along(chars)) {
    if (chars[i] %in% c("[", "(")) depth <- depth + 1L
    if (chars[i] %in% c("]", ")")) depth <- depth - 1L
    if (chars[i] == "," && depth == 0L) {
      grp <- grp + 1L
      next # which_grp stays 0, so the separator itself is dropped
    }
    which_grp[i] <- grp
  }
  keep <- which_grp > 0L
  parts <- split(chars[keep], factor(which_grp[keep], levels = seq_len(grp)))
  trimws(unname(vapply(parts, paste, character(1), collapse = "")))
}


# `lower` and `upper` out of a `<...>` bounds expression, as numbers.
#' @keywords internal
.stan_bounds <- function(bounds, sdata) {
  out <- c(lower = NA_real_, upper = NA_real_)
  if (!nzchar(bounds)) return(out)
  for (part in .split_top_commas(bounds)) {
    kv <- regmatches(part, regexec("^\\s*(lower|upper)\\s*=\\s*(.+)$", part))[[1]]
    if (length(kv) == 3) out[[kv[2]]] <- .eval_stan_expr(kv[3], sdata)
  }
  out
}


# A Stan expression over the data block, as a number. `enclos = baseenv()` so
# that arithmetic resolves but nothing from the caller's workspace does.
#' @keywords internal
.eval_stan_expr <- function(expr, sdata) {
  v <- tryCatch(
    eval(parse(text = expr), envir = sdata, enclos = baseenv()),
    error = function(e) NULL
  )
  if (!is.numeric(v) || length(v) != 1 || is.na(v)) return(NA_real_)
  as.numeric(v)
}


# Natural scale -> link scale. Returns NA for a link with no inverse here,
# which leaves the parameter on its generic default rather than erroring: an
# unhelpful starting value is a smaller problem than a function that refuses to
# run.
#' @keywords internal
.apply_link <- function(x, link) {
  if (is.null(link) || is.na(link)) return(NA_real_)
  out <- switch(
    link,
    identity = x,
    log = log(x),
    logm1 = log(x - 1),
    log1p = log1p(x),
    inverse = 1 / x,
    sqrt = sqrt(x),
    logit = stats::qlogis(x),
    probit = , probit_approx = stats::qnorm(x),
    cloglog = log(-log1p(-x)),
    cauchit = stats::qcauchy(x),
    softplus = log(expm1(x)),
    # squareplus(z) = (z + sqrt(z^2 + 4)) / 2, whose inverse is x - 1 / x
    squareplus = x - 1 / x,
    tan_half = tan(x / 2),
    NA_real_
  )
  if (!is.numeric(out) || length(out) != 1 || !is.finite(out)) NA_real_ else out
}
