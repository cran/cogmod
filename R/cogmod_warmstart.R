#' Warm-start a fit from a previous one: metric, step size and starting values
#'
#' @description
#' Takes what warmup produced in a previous fit - the adapted inverse metric,
#' the step size and the posterior means - and turns it into the `inv_metric`,
#' `step_size` and `init` arguments of a new [brms::brm()] call, so that the
#' new run can get by with a much shorter warmup. The previous fit can be the
#' **same model** (a refit with more draws, another seed, a slightly different
#' prior) or a **pilot on a subset of the participants**: the full model then
#' has more parameters, one standardized random effect per new participant per
#' group-level term, and the metric is carried over **by parameter name**, with
#' a sensible filler for what the pilot never saw.
#'
#' ```r
#' ws <- cogmod_warmstart(pilot, data = data)   # the pilot's model, on all the data
#' m <- brm(formula, data = data, prior = ..., stanvars = ...,
#'          init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size,
#'          warmup = 100, iter = 600, backend = "cmdstanr")
#' ```
#'
#' The sampler keeps adapting from the supplied values during whatever warmup
#' remains, so a poor warm start costs speed, not correctness. Nothing about it
#' changes the posterior being sampled.
#'
#' @details
#' # What is carried over, and how
#'
#' Stan adapts one variance per **unconstrained** scalar parameter, in the
#' order of the program's `parameters` block, and `brms` keeps those variances
#' (one vector per chain, averaged here) and the step size in the fit's
#' metadata. To move them to another model they are first labelled with the
#' Stan parameter names - `Intercept`, `sd_1[1]`, `z_1[1,3]`, and so on - which
#' are read off the generated program, the same way [cogmod_inits()] does it.
#' The target model's labels are built the same way, and the two are joined:
#'
#' - **Population-level entries** (coefficients, intercepts, group-level SDs
#'   and correlations, auxiliary parameters) carry over one to one.
#' - **Standardized group-level effects** (`z_<k>[m, n]`) are indexed by
#'   participant, so a pilot participant's entry moves to that participant's
#'   position in the target model, matched by the level's *name*. Participants
#'   the pilot never saw get the average of the pilot's entries for the same
#'   term (they are standardized effects, so that is a fair guess) and start
#'   at zero.
#' - **Anything without a counterpart** - a predictor or a group-level term
#'   the pilot did not have - gets Stan's default variance of 1 and the generic
#'   starting value [cogmod_inits()] would give it. A count is reported by
#'   `print()`; a large one usually means the two formulas differ more than
#'   intended.
#'
#' Starting values are the pilot's posterior means. `brms` drops the raw `z_`
#' and Cholesky factors from a saved fit, so these are rebuilt from what it
#' keeps: the group-level effects `r_`, their SDs and their correlations. Each
#' chain gets the same values plus a little noise on the unconstrained scale
#' (`jitter`), so that the chains do not start at one point and `Rhat` keeps
#' some meaning. Note that tightly initialised chains are less likely to find
#' a second mode than dispersed ones; if that is a concern, run the cold start
#' once.
#'
#' # Same model, or a different one
#'
#' Whatever is not given is taken from the source fit. With `formula` and
#' `data` both left `NULL`, the target is the fit itself: the result
#' reproduces its own adaptation, and is the way to refit with fewer warmup
#' iterations. With only `data`, the target is the same model on that data -
#' the pilot-to-full-sample case. With only `formula`, it is that model on the
#' source's data - a variant of the model, say with one more predictor, whose
#' shared parameters can start where the first fit left them. A table or file
#' source carries neither and needs both. Any `brms` model fitted with the
#' `cmdstanr` backend and the default diagonal metric can be a source; a
#' `dense_e` fit is refused, and so is the `rstan` backend, which does not
#' store the adaptation.
#'
#' The `inv_metric` returned is a plain vector, as `cmdstanr` wants it; the
#' labels are in `ws$table`. The metric must match the target program exactly,
#' which is why `formula` and `data` are needed rather than just a count of
#' participants: `brms` decides the layout from both.
#'
#' # Storing it, and the four helpers
#'
#' `as.data.frame()` gives a small table (one row per unconstrained parameter:
#' label, group and level, variance, posterior mean, median and SD, step size)
#' that can be written with [utils::write.csv()], so that a pilot fitted on a
#' laptop can warm-start an array job on a cluster with a file of a few
#' kilobytes and no `brmsfit` in sight. On the other side, four functions with
#' the same signature - the model's formula and data first, the source under
#' `warmstart` - each give one argument of the `brm()` call:
#'
#' ```r
#' tab <- read.csv("pilot_warmstart.csv")   # or the path, or the brmsfit itself
#' m <- brm(formula, data = data, stanvars = ...,
#'          prior = cogmod_priors(formula, data, warmstart = tab),
#'          init = cogmod_inits(formula, data, warmstart = tab),
#'          inv_metric = cogmod_inv_metric(formula, data, warmstart = tab),
#'          step_size = cogmod_step_size(formula, data, warmstart = tab),
#'          warmup = 100, iter = 600, backend = "cmdstanr")
#' ```
#'
#' The first of those is the odd one out and is **not** part of a warm start
#' in the sense the rest of this page uses. [cogmod_priors()] re-centres the
#' priors on the source's posterior median and SD, which changes the model
#' rather than the path the sampler takes through it - and double-counts the
#' source's data if the new model contains it. Its own documentation says when
#' that is and is not legitimate; the other three change nothing about the
#' posterior being sampled.
#'
#' Each maps the table onto the model `formula` and `data` describe, so it
#' does not matter which model the table was written for: a table from a pilot
#' on fewer participants is extended, one written for another formula falls
#' back to the defaults with a note. When the table was made for this very
#' model, `tab$inv_metric` and `tab$step_size[1]` are the same numbers (the
#' step size is one number repeated down the column; a whole column there
#' would be read as one step size per chain).
#'
#' # What it is worth
#'
#' On an LNR and a DDM with participant random intercepts, a pilot on 4 of 8
#' participants warm-started the full fit to about twice the effective draws
#' per second of a cold start with a 500-iteration warmup, and four to six
#' times those of a cold start with the same 100-iteration warmup. The
#' starting values alone
#' bought nothing: what a short warmup lacks is an adapted step size and
#' metric, not a good position. Details and the benchmark are in
#' the [performance article](https://dominiquemakowski.github.io/cogmod/articles/performance.html).
#'
#' @param x The source: a `brmsfit` fitted with `backend = "cmdstanr"`, a
#'   `cogmod_warmstart` object, the data frame `as.data.frame()` makes of one,
#'   or the path to a CSV file holding that data frame.
#' @param formula,data The target model, as they will be passed to
#'   [brms::brm()]. Either left `NULL` (the default) is taken from the source
#'   fit: the same formula on new data, the same data under a new formula, or
#'   with both `NULL` the source fit itself. That requires `x` to be a
#'   `brmsfit`; a table or file source needs both.
#' @param jitter SD of the noise added to the starting values on the
#'   unconstrained scale, so that chains start at different points. Smaller than
#'   [cogmod_inits()]'s default because the values come from a converged
#'   posterior; `0` gives identical starts. As there, one number is the SD for
#'   the population-level blocks and the group-level and smooth blocks get a
#'   fifth of it; two numbers set the two tiers directly.
#' @param ... Passed to [brms::make_stancode()], [brms::make_standata()] and
#'   [brms::brm()] (with `empty = TRUE`) when the target model is built, for
#'   arguments such as `data2`.
#'
#' @return An object of class `cogmod_warmstart`: a list with
#' \describe{
#'   \item{`inv_metric`}{Numeric vector, one variance per unconstrained
#'     parameter of the target model, in Stan's order.}
#'   \item{`step_size`}{The step size, averaged over the source's chains.}
#'   \item{`init`}{A function of one argument, for `brms::brm(init = )`,
#'     returning a named list of starting values for every parameter the
#'     target program declares.}
#'   \item{`table`}{A data frame with one row per unconstrained parameter:
#'     `parameter` (the Stan label), `group`, `coef` and `level` (for a
#'     group-level parameter, the grouping factor, the coefficient and, for a
#'     standardized effect, the level it stands for; otherwise `NA`),
#'     `inv_metric`, `mean` (the source posterior mean, `NA` where none
#'     applies) and `step_size` (the same value in every row). This is what
#'     `as.data.frame()` returns.}
#'   \item{`counts`, `missing`}{How many entries came from the source, are new
#'     group levels, or have no counterpart, and the names of the latter; what
#'     `print()` reports.}
#' }
#'
#' @seealso [cogmod_inits()], which supplies the starting values used where
#'   the source has none, and the [performance article](https://dominiquemakowski.github.io/cogmod/articles/performance.html).
#'
#' @examples
#' \donttest{
#' # Fitting needs cmdstanr, which lives outside CRAN - see the package website.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   df <- data.frame(
#'     RT = rcogmod_lognormal(400, ndt = 0.3, poutlier = 0.02),
#'     id = factor(rep(1:8, each = 50))
#'   )
#'   f <- brms::bf(RT ~ 1 + (1 | id), ndt ~ 1, poutlier ~ 1,
#'     family = cogmod_lognormal()
#'   )
#'
#'   # A pilot on some participants, then the full sample
#'   pilot_df <- droplevels(df[df$id %in% 1:3, ])
#'   pilot <- brms::brm(f,
#'     data = pilot_df, prior = cogmod_priors(f, pilot_df),
#'     init = cogmod_inits(f, pilot_df), stanvars = cogmod_stanvars(f),
#'     backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
#'   )
#'   ws <- cogmod_warmstart(pilot, data = df) # same model, all the data
#'   print(ws) # how much of the metric came from the pilot
#'   m <- brms::brm(f,
#'     data = df, prior = cogmod_priors(f, df), stanvars = cogmod_stanvars(f),
#'     init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size,
#'     backend = "cmdstanr", chains = 1, warmup = 100, iter = 300, refresh = 0
#'   )
#'
#'   # Keep it as a table, e.g. for a cluster...
#'   tab <- tempfile(fileext = ".csv")
#'   write.csv(as.data.frame(ws), tab, row.names = FALSE)
#'   # ... and there, one helper per brm() argument, all with the same signature
#'   init <- cogmod_inits(f, df, warmstart = tab)
#'   inv_metric <- cogmod_inv_metric(f, df, warmstart = tab)
#'   step_size <- cogmod_step_size(f, df, warmstart = tab)
#'
#'   # The fourth helper is a different kind of thing: it moves the PRIORS onto
#'   # the pilot's posterior, which changes the model rather than the sampler.
#'   # Only where the pilot is independent of `df` - see ?cogmod_priors.
#'   print(cogmod_priors(f, df, warmstart = tab))
#'   unlink(tab)
#' }
#' }
#'
#' @export
cogmod_warmstart <- function(x, formula = NULL, data = NULL, jitter = 0.05, ...) {
  src <- .warmstart_source(x)

  # Whatever is not given describes the source fit: no formula means the same
  # model (on new data, if any), no data means the same data (under a new
  # formula, if any), neither means a refit of the source itself. Only a
  # brmsfit carries them; a table or file needs both spelled out.
  if (is.null(formula) || is.null(data)) {
    if (!inherits(x, "brmsfit")) {
      missing <- c("`formula`", "`data`")[c(is.null(formula), is.null(data))]
      stop("Give ", paste(missing, collapse = " and "), ": only a brmsfit ",
           "carries the source model's own, and `x` is not one.", call. = FALSE)
    }
    if (is.null(formula)) formula <- x$formula
    if (is.null(data)) data <- x$data
  }
  target <- .warmstart_target(formula, data, ...)
  tab <- .warmstart_join(src, target$table)
  init <- .init_fun(.warmstart_apply_means(target$plan, tab), jitter)

  # Where each entry came from is reported, not stored: the table is what
  # gets written to a file and read back as a source, and it should carry
  # nothing that is not about the target model.
  counts <- table(factor(tab$source, levels = c("pilot", "new level", "default")))
  missing <- unique(sub("\\[.*$", "", tab$parameter[tab$source == "default"]))
  tab$source <- NULL

  prior <- .warmstart_prior_table(tab)
  tab[c("pclass", "pdpar", "pcoef")] <- NULL

  structure(
    list(inv_metric = unname(tab$inv_metric), step_size = tab$step_size[1],
         init = init, prior = prior, table = tab,
         counts = stats::setNames(as.integer(counts), names(counts)), missing = missing),
    class = "cogmod_warmstart"
  )
}


# The posterior behind each prior slot the source can speak to: one row per
# brms prior row that has a parameter with a posterior behind it, addressed
# the way get_prior() addresses it. Everything else drops out - the `z_` and
# `L_` entries, which have no stated prior; a parameter the source never saw,
# which has no posterior; and a parameter whose posterior SD came out at zero
# or missing, where `normal(m, 0)` would be a point mass rather than a prior.
# Only the four classes a Normal is the right shape for are kept: a location
# on its own scale (an intercept, an auxiliary dpar), a coefficient, or a
# group-level SD, whose prior brms does not require to be positive-only.
#' @keywords internal
.warmstart_prior_table <- function(tab) {
  keep <- !is.na(tab$pclass) & !is.na(tab$pcoef) &
    is.finite(tab$median) & is.finite(tab$sd) & tab$sd > 0
  # An auxiliary dpar's class is its own name, so the whitelist is on shape,
  # not on a list of names: what is excluded is `cor`, `L` and `sds`, none of
  # which a Normal describes.
  keep <- keep & !tab$pclass %in% c("cor", "L", "sds", "z")
  data.frame(
    class = tab$pclass[keep], dpar = tab$pdpar[keep], coef = tab$pcoef[keep],
    group = ifelse(is.na(tab$group[keep]), "", tab$group[keep]),
    median = tab$median[keep], sd = tab$sd[keep],
    stringsAsFactors = FALSE
  )
}


#' @rdname cogmod_warmstart
#' @param warmstart The source, as `x` above: a `brmsfit`, a `cogmod_warmstart`
#'   object, its data frame, or the path to a CSV file of it.
#' @export
cogmod_inv_metric <- function(formula = NULL, data = NULL, warmstart, ...) {
  if (missing(warmstart)) {
    stop("`warmstart` is required: the fit, table or file to take the metric from.",
         call. = FALSE)
  }
  cogmod_warmstart(warmstart, formula = formula, data = data, jitter = 0, ...)$inv_metric
}


#' @rdname cogmod_warmstart
#' @export
cogmod_step_size <- function(formula = NULL, data = NULL, warmstart, ...) {
  if (missing(warmstart)) {
    stop("`warmstart` is required: the fit, table or file to take the step size from.",
         call. = FALSE)
  }
  # The step size is a property of the source alone, so the target model is
  # accepted for symmetry with the other helpers and not needed.
  .warmstart_source(warmstart)$step_size[1]
}


#' @rdname cogmod_warmstart
#' @param row.names,optional Ignored; present for compatibility with the
#'   [as.data.frame()] generic.
#' @export
as.data.frame.cogmod_warmstart <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table[, c("parameter", "group", "coef", "level", "inv_metric", "mean",
              "median", "sd", "step_size")]
}


#' @rdname cogmod_warmstart
#' @export
print.cogmod_warmstart <- function(x, ...) {
  n <- nrow(x$table)
  from <- x$counts[["pilot"]]
  newl <- x$counts[["new level"]]
  dflt <- x$counts[["default"]]
  cat(sprintf("<cogmod_warmstart> %d unconstrained parameters, step size %.3g\n", n, x$step_size))
  cat(sprintf("  %d carried over from the source", from))
  if (newl) cat(sprintf(", %d for new group levels (term average, start at 0)", newl))
  if (dflt) cat(sprintf(", %d without a counterpart (variance 1, generic start)", dflt))
  cat("\n")
  if (dflt) cat("  without counterpart:", paste(x$missing, collapse = ", "), "\n")
  if (nrow(x$prior)) {
    cat(sprintf("  %d of them carry a posterior median and SD, for cogmod_priors(warmstart = )\n",
                nrow(x$prior)))
  }
  cat("  Use: brm(..., init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size)\n")
  invisible(x)
}


# The source table -----------------------------------------------------------

# One row per unconstrained parameter of the source, however it was given.
#' @keywords internal
.warmstart_source <- function(x) {
  if (inherits(x, "cogmod_warmstart")) return(x$table)
  if (inherits(x, "brmsfit")) return(.warmstart_extract(x))
  if (is.character(x) && length(x) == 1) {
    if (!file.exists(x)) stop("No such file: ", x, call. = FALSE)
    x <- utils::read.csv(x, stringsAsFactors = FALSE)
  }
  if (is.data.frame(x)) {
    need <- c("parameter", "inv_metric", "step_size")
    miss <- setdiff(need, names(x))
    if (length(miss)) {
      stop("A warm-start table needs the columns ", paste(need, collapse = ", "),
           "; missing: ", paste(miss, collapse = ", "), ". Write one with ",
           "write.csv(as.data.frame(cogmod_warmstart(fit)), file, row.names = FALSE).",
           call. = FALSE)
    }
    for (col in c("group", "coef", "level")) {
      if (is.null(x[[col]])) x[[col]] <- NA_character_
      x[[col]] <- as.character(x[[col]])
      x[[col]][!is.na(x[[col]]) & !nzchar(x[[col]])] <- NA_character_
    }
    # `median` and `sd` came later than the rest, and only cogmod_priors()
    # wants them, so a table written before they existed still reads.
    for (col in c("mean", "median", "sd")) {
      if (is.null(x[[col]])) x[[col]] <- NA_real_
    }
    x$parameter <- as.character(x$parameter)
    return(x)
  }
  stop("`x` must be a brmsfit, a cogmod_warmstart object, its data frame, or ",
       "the path to a CSV of it.", call. = FALSE)
}


# Metric, step size and posterior means out of a brmsfit.
#' @keywords internal
.warmstart_extract <- function(fit) {
  md <- attr(fit$fit, "metadata")
  if (is.null(md) || is.null(md$inv_metric) || is.null(md$step_size)) {
    stop("The fit carries no adaptation. cogmod_warmstart() needs a model ",
         "sampled with backend = \"cmdstanr\" (the rstan backend does not ",
         "keep the metric, and approximations have none).", call. = FALSE)
  }
  if (is.matrix(md$inv_metric[[1]])) {
    stop("The fit used a dense metric (metric = \"dense_e\"). Only the ",
         "diagonal metric can be carried over.", call. = FALSE)
  }
  inv_metric <- Reduce(`+`, md$inv_metric) / length(md$inv_metric)
  step_size <- mean(unlist(md$step_size))

  decls <- .stan_param_decls(fit$model)
  sdata <- as.list(brms::standata(fit))
  levels <- attr(fit$ranef, "levels")
  tab <- .warmstart_labels(decls, sdata, fit$ranef, levels)
  if (nrow(tab) != length(inv_metric)) {
    stop("The fit's metric has ", length(inv_metric), " entries but its ",
         "parameters block was read as ", nrow(tab), " unconstrained scalars. ",
         "The program uses a declaration cogmod_warmstart() cannot size.",
         call. = FALSE)
  }
  tab$inv_metric <- inv_metric
  st <- .warmstart_stats(fit, decls, sdata, tab)
  tab$mean <- st$mean
  tab$median <- st$median
  tab$sd <- st$sd
  tab$step_size <- step_size
  tab
}


# The target: its labels, and the init plan cogmod_inits() would use, so that
# whatever the source cannot supply gets the same generic values.
#' @keywords internal
.warmstart_target <- function(formula, data, ...) {
  family <- tryCatch(.cogmod_family(formula), error = function(e) NULL)
  targets <- tryCatch(.init_targets(family), error = function(e) NULL)
  if (is.null(targets)) targets <- list()
  links <- tryCatch(.family_links(family), error = function(e) NULL)

  # `empty = TRUE` validates the formula against the data, writes the program
  # and works out the group-level structure - levels included - without
  # compiling anything.
  skeleton <- suppressWarnings(
    brms::brm(formula, data = data, empty = TRUE, ...)
  )
  decls <- .stan_param_decls(skeleton$model)
  sdata <- as.list(suppressWarnings(brms::make_standata(formula, data = data, ...)))
  tab <- .warmstart_labels(decls, sdata, skeleton$ranef, attr(skeleton$ranef, "levels"))
  list(table = tab, plan = .init_plan(decls, sdata, targets, links))
}


# Labels ----------------------------------------------------------------------

# One row per unconstrained scalar, in Stan's storage order: declaration
# order, and within a declaration the layout Stan serializes it in. A
# `matrix[M, N]` is column-major; an `array[M] vector[N]` has the array index
# outermost; a `cholesky_factor_corr[K]` has K(K-1)/2 free entries, the strict
# lower triangle by rows. The label of a free entry that is not itself a
# constrained element (the Cholesky factor, a simplex) is only a name to join
# on, and is the same on both sides.
#
# `z_<k>` rows also carry the group level their column stands for, which is
# what lets a participant's entry follow it into a differently sized model.
#' @keywords internal
.warmstart_labels <- function(decls, sdata, ranef, levels) {
  rows <- lapply(decls, function(d) {
    dims <- vapply(d$dims, .eval_stan_expr, numeric(1), sdata = sdata)
    if (anyNA(dims)) {
      stop("Cannot size the declaration of `", d$name, "` from the data.", call. = FALSE)
    }
    dims <- as.integer(dims)
    lab <- .unconstrained_labels(d, dims)
    if (!length(lab)) return(NULL)
    n_lab <- length(lab)
    group <- coef <- level <- rep(NA_character_, n_lab)
    # The brms prior slot the entry belongs to, where it has one: the `class`,
    # `dpar` and `coef` a get_prior() row would carry for it, with `group`
    # above completing the address. This is what lets cogmod_priors(warmstart =)
    # find the posterior of the parameter a prior row is about. Left NA for
    # everything whose prior is fixed by the parameterization rather than
    # stated - the standardized effects `z_` and the Cholesky factors `L_`.
    pclass <- pdpar <- pcoef <- rep(NA_character_, n_lab)
    if (identical(d$type, "real") && !length(d$dims)) {
      # A bare real is a centred intercept or, for a dpar left out of the
      # formula, the auxiliary parameter itself, whose prior class is its name.
      if (identical(d$name, "Intercept")) {
        pclass <- "Intercept"
        pdpar <- ""
      } else if (startsWith(d$name, "Intercept_")) {
        pclass <- "Intercept"
        pdpar <- sub("^Intercept_", "", d$name)
      } else {
        pclass <- d$name
        pdpar <- ""
      }
      pcoef <- ""
    } else if (identical(d$type, "vector") &&
               (identical(d$name, "b") || grepl("^b_", d$name))) {
      dp <- if (identical(d$name, "b")) "" else sub("^b_", "", d$name)
      pclass <- "b"
      pdpar <- dp
      pcoef <- .b_coef_names(sdata, if (nzchar(dp)) dp else "mu", n_lab)
    }
    # Group-level parameters carry the grouping factor, the coefficient and
    # (for a standardized effect) the level they stand for, so that they can
    # be matched on what they mean rather than on where brms happened to put
    # them: sd_<k>[m] is coefficient m of term k, z_<k>[m, n] its value for
    # level n, and L_<k>[i, j] the entry for coefficients j and i.
    if (grepl("^(sd|z|L)_[0-9]+$", d$name) && !is.null(ranef)) {
      k <- as.integer(sub("^[a-zA-Z]+_", "", d$name))
      rows <- ranef[ranef$id == k, , drop = FALSE]
      if (nrow(rows)) {
        cf <- .re_coef(rows)
        idx <- .label_indices(lab)
        m <- vapply(idx, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
        if (!anyNA(m) && max(m) <= length(cf)) {
          group[] <- rows$group[1]
          if (startsWith(d$name, "L_")) {
            j <- vapply(idx, `[`, integer(1), 2)
            coef <- paste0(cf[j], "__", cf[m])
          } else {
            coef <- cf[m]
          }
          if (startsWith(d$name, "z_")) {
            lv <- levels[[rows$group[1]]]
            n <- vapply(idx, `[`, integer(1), 2)
            if (!is.null(lv) && !anyNA(n) && max(n) <= length(lv)) level <- lv[n]
          }
          # A group-level SD is the one of the three that has a prior of its
          # own, addressed by the bare coefficient and the dpar separately
          # rather than by the joined name `coef` carries.
          if (startsWith(d$name, "sd_")) {
            pclass <- "sd"
            pdpar <- ifelse(nzchar(rows$dpar[m]), rows$dpar[m], rows$nlpar[m])
            pcoef <- rows$coef[m]
          }
        }
      }
    }
    data.frame(parameter = lab, group = group, coef = coef, level = level,
               pclass = pclass, pdpar = pdpar, pcoef = pcoef,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(parameter = character(0), group = character(0),
                      coef = character(0), level = character(0),
                      pclass = character(0), pdpar = character(0),
                      pcoef = character(0)))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# The names brms gives the `k` entries of a population-level design matrix,
# which are both the `coef` of their prior rows and the tail of their
# `b_<dpar>_<coef>` variable names. NA if the matrix cannot be found or does
# not have `k` columns to spare - every caller then falls back rather than
# guessing. brms drops the intercept column from `X` when it centres it, but
# keeps its name, hence the second try.
#' @keywords internal
.b_coef_names <- function(sdata, dpar, k) {
  X <- sdata[[if (identical(dpar, "mu")) "X" else paste0("X_", dpar)]]
  cn <- colnames(X)
  if (is.null(cn)) return(rep(NA_character_, k))
  if (length(cn) != k) cn <- cn[cn != "Intercept"]
  if (length(cn) != k) return(rep(NA_character_, k))
  cn
}


# The bracketed indices of a label, integer(0) for a bare name.
#' @keywords internal
.label_indices <- function(lab) {
  lapply(lab, function(p) {
    if (!grepl("[", p, fixed = TRUE)) return(integer(0))
    as.integer(strsplit(sub("^.*\\[(.*)\\]$", "\\1", p), ",", fixed = TRUE)[[1]])
  })
}


# The name brms gives a group-level coefficient: the coefficient itself,
# prefixed by the dpar or nlpar (and the response, in a multivariate model) it
# belongs to, as in `nuone_Intercept`.
#' @keywords internal
.re_coef <- function(rows) {
  pre <- ifelse(nzchar(rows$dpar), rows$dpar, rows$nlpar)
  if (!is.null(rows$resp)) {
    pre <- ifelse(nzchar(rows$resp), paste0(rows$resp, ifelse(nzchar(pre), "_", ""), pre), pre)
  }
  ifelse(nzchar(pre), paste0(pre, "_", rows$coef), rows$coef)
}


# The labels of one declaration's free entries, in Stan's order.
#' @keywords internal
.unconstrained_labels <- function(d, dims) {
  nm <- d$name
  # How many of the sizes belong to the type itself; the rest are array
  # dimensions in front of it. (`cholesky_factor_cov` can take one or two
  # sizes, which makes it ambiguous inside an array; brms never emits it, so
  # it is taken as two only when it is not an array.)
  type_dims <- switch(
    d$type, real = 0L,
    vector = , row_vector = , simplex = , ordered = , positive_ordered = ,
    unit_vector = , cholesky_factor_corr = , corr_matrix = , cov_matrix = 1L,
    matrix = 2L,
    cholesky_factor_cov = min(2L, length(dims))
  )
  n_array <- length(dims) - type_dims
  if (n_array < 0L) {
    stop("Too few sizes in the declaration of `", nm, "`.", call. = FALSE)
  }
  adims <- dims[seq_len(n_array)]
  idims <- if (n_array) dims[-seq_len(n_array)] else dims

  inner <- switch(
    d$type,
    real = list(),
    vector = , row_vector = , ordered = , positive_ordered = , unit_vector =
      lapply(seq_len(idims[1]), function(i) i),
    simplex = lapply(seq_len(max(idims[1] - 1L, 0L)), function(i) i),
    matrix = {
      g <- expand.grid(m = seq_len(idims[1]), n = seq_len(idims[2]))  # column-major
      lapply(seq_len(nrow(g)), function(i) c(g$m[i], g$n[i]))
    },
    cholesky_factor_corr = , corr_matrix = {
      K <- idims[1]
      out <- list()
      for (i in seq_len(K)[-1]) for (j in seq_len(i - 1L)) out[[length(out) + 1L]] <- c(i, j)
      out
    },
    cov_matrix = {
      K <- idims[1]
      lapply(seq_len(K * (K + 1L) / 2L), function(i) i)
    },
    cholesky_factor_cov = {
      M <- idims[1]
      N <- if (length(idims) > 1) idims[2] else M
      lapply(seq_len(N * (N + 1L) / 2L + (M - N) * N), function(i) i)
    },
    stop("Unsupported parameter type `", d$type, "` for `", nm, "`.", call. = FALSE)
  )
  fmt <- function(idx) if (!length(idx)) nm else sprintf("%s[%s]", nm, paste(idx, collapse = ","))
  scalar <- identical(d$type, "real")

  if (!n_array) {
    if (scalar) return(nm)
    # a zero-length declaration (`vector[Kc] b` with Kc = 0) has no entries
    if (!length(inner)) return(character(0))
    return(vapply(inner, fmt, character(1)))
  }
  if (any(adims == 0L) || (!scalar && !length(inner))) return(character(0))
  # Array indices outermost, first index slowest (row-major over the array
  # dimensions), then the inner layout.
  g <- expand.grid(rev(lapply(adims, seq_len)))[, rev(seq_along(adims)), drop = FALSE]
  out <- character(0)
  for (r in seq_len(nrow(g))) {
    a <- unlist(g[r, ])
    if (scalar) out <- c(out, fmt(a))
    else out <- c(out, vapply(inner, function(idx) fmt(c(a, idx)), character(1)))
  }
  out
}


# Posterior summaries ---------------------------------------------------------

# The source posterior mean, median and SD for each labelled entry, NA where
# none applies. brms renames what it saves, and drops `z_` and `L_`, so the
# values are looked up under the saved names and the dropped ones rebuilt: a
# standardized effect from the saved group-level effects r = sd * (L z), a
# Cholesky factor from the saved correlations. Every lookup is by name and
# returns NA when the name is not there, so an unforeseen naming scheme
# degrades to the generic starting value rather than to an error.
#
# Only the mean is rebuilt for `z_` and `L_`. Reconstructing them takes a
# matrix solve and a Cholesky decomposition, which are not operations one can
# apply to a median or an SD and get the median or SD of the result; the
# starting values are all those two entries are wanted for, so their `median`
# and `sd` stay NA rather than being filled with something that is not one.
#' @keywords internal
.warmstart_stats <- function(fit, decls, sdata, tab) {
  vars <- brms::variables(fit)
  draws <- brms::as_draws_matrix(fit)
  pm <- colMeans(draws)
  get <- function(nms) {
    out <- rep(NA_real_, length(nms))
    ok <- nms %in% vars
    out[ok] <- pm[nms[ok]]
    out
  }
  # The three summaries at once, for the entries that are a saved variable and
  # so can be summarized directly.
  summaries <- rbind(
    mean = pm,
    median = apply(draws, 2, stats::median),
    sd = apply(draws, 2, stats::sd)
  )
  get3 <- function(nms) {
    out <- matrix(NA_real_, length(nms), 3,
                  dimnames = list(NULL, c("mean", "median", "sd")))
    ok <- !is.na(nms) & nms %in% vars
    if (any(ok)) out[ok, ] <- t(summaries[, nms[ok], drop = FALSE])
    out
  }
  ranef <- fit$ranef
  out <- matrix(NA_real_, nrow(tab), 3,
                dimnames = list(NULL, c("mean", "median", "sd")))
  base <- sub("\\[.*$", "", tab$parameter)

  coef_names <- function(dpar, k) .b_coef_names(sdata, dpar, k)
  re_coef <- .re_coef
  re_suffix <- function(rows) {
    pre <- ifelse(nzchar(rows$dpar), rows$dpar, rows$nlpar)
    pre <- ifelse(nzchar(rows$resp), paste0(rows$resp, ifelse(nzchar(pre), "_", ""), pre), pre)
    ifelse(nzchar(pre), paste0("__", pre), "")
  }

  for (d in decls) {
    i <- which(base == d$name)
    if (!length(i)) next
    nm <- d$name
    if (identical(d$type, "real") && !length(d$dims)) {
      out[i, ] <- get3(nm)   # Intercept, Intercept_<dpar>, auxiliary dpars
    } else if (identical(d$type, "vector") && (identical(nm, "b") || grepl("^b_", nm))) {
      dpar <- if (identical(nm, "b")) "mu" else sub("^b_", "", nm)
      cn <- coef_names(dpar, length(i))
      out[i, ] <- get3(paste0(if (identical(dpar, "mu")) "b_" else paste0("b_", dpar, "_"), cn))
    } else if (grepl("^sd_[0-9]+$", nm)) {
      rows <- ranef[ranef$id == as.integer(sub("^sd_", "", nm)), , drop = FALSE]
      if (nrow(rows) == length(i)) {
        out[i, ] <- get3(paste0("sd_", rows$group, "__", re_coef(rows)))
      }
    } else if (grepl("^L_[0-9]+$", nm)) {
      L <- .warmstart_L(fit, ranef, as.integer(sub("^L_", "", nm)), get, re_coef)
      if (!is.null(L)) {
        ij <- do.call(rbind, lapply(strsplit(sub("^.*\\[(.*)\\]$", "\\1", tab$parameter[i]), ","), as.integer))
        out[i, "mean"] <- L[ij]
      }
    } else if (grepl("^z_[0-9]+$", nm)) {
      k <- as.integer(sub("^z_", "", nm))
      rows <- ranef[ranef$id == k, , drop = FALSE]
      lv <- attr(ranef, "levels")[[rows$group[1]]]
      M <- nrow(rows)
      if (M && !is.null(lv)) {
        sds <- get(paste0("sd_", rows$group, "__", re_coef(rows)))
        # r_<group>[<level>,<coef>], or r_<group>__<dpar>[<level>,<coef>]
        R <- matrix(NA_real_, M, length(lv))
        for (m in seq_len(M)) {
          R[m, ] <- get(sprintf("r_%s%s[%s,%s]", rows$group[m], re_suffix(rows)[m], lv, rows$coef[m])) / sds[m]
        }
        L <- if (M > 1) .warmstart_L(fit, ranef, k, get, re_coef) else matrix(1, 1, 1)
        if (!is.null(L) && !anyNA(R)) {
          Z <- solve(L, R)
          ij <- do.call(rbind, lapply(strsplit(sub("^.*\\[(.*)\\]$", "\\1", tab$parameter[i]), ","), as.integer))
          out[i, "mean"] <- Z[ij]
        }
      }
    }
    # anything else stays NA and gets the generic starting value
  }
  as.data.frame(out)
}


# The Cholesky factor of a group-level correlation matrix, from the saved
# `cor_<group>__<coef_i>__<coef_j>` means; NULL if any is missing or the
# matrix they make is not positive definite.
#' @keywords internal
.warmstart_L <- function(fit, ranef, k, get, re_coef) {
  rows <- ranef[ranef$id == k, , drop = FALSE]
  M <- nrow(rows)
  if (!M) return(NULL)
  C <- diag(M)
  if (M > 1) {
    cf <- re_coef(rows)
    for (i in seq_len(M)[-1]) for (j in seq_len(i - 1L)) {
      v <- get(paste0("cor_", rows$group[1], "__", cf[j], "__", cf[i]))
      if (is.na(v)) return(NULL)
      C[i, j] <- C[j, i] <- v
    }
  }
  L <- tryCatch(t(chol(C)), error = function(e) NULL)
  L
}


# Joining source and target ------------------------------------------------

# The target table with the source's variances and means moved onto it.
#' @keywords internal
.warmstart_join <- function(src, target) {
  # Population-level entries match on their label. Group-level entries match
  # on what they stand for - kind, grouping factor, coefficient and level -
  # so that a term brms numbers differently in the two programs, or a
  # participant it indexes differently, still finds its counterpart.
  key <- function(tab, with_level = TRUE) {
    k <- tab$parameter
    re <- !is.na(tab$group)
    kind <- sub("_[0-9]+$", "", sub("\\[.*$", "", tab$parameter[re]))
    lv <- if (with_level) ifelse(is.na(tab$level[re]), "", tab$level[re]) else ""
    k[re] <- paste(kind, tab$group[re], tab$coef[re], lv, sep = "|")
    k
  }
  idx <- match(key(target), key(src))
  hit <- !is.na(idx)

  out <- target
  for (col in c("inv_metric", "mean", "median", "sd")) {
    out[[col]] <- rep(NA_real_, nrow(out))
    out[[col]][hit] <- src[[col]][idx[hit]]
  }
  out$source <- ifelse(hit, "pilot", "default")

  # A new level of an effect the source had: the average over the source's
  # levels of that effect, and a start at zero.
  is_z <- grepl("^z_[0-9]+\\[", out$parameter) & !is.na(out$level)
  new_z <- !hit & is_z
  if (any(new_z)) {
    k_t <- key(out, with_level = FALSE)
    k_s <- key(src, with_level = FALSE)
    for (kk in unique(k_t[new_z])) {
      from <- src$inv_metric[k_s == kk]
      if (!length(from)) next
      sel <- new_z & k_t == kk
      out$inv_metric[sel] <- mean(from, na.rm = TRUE)
      out$mean[sel] <- 0
      out$source[sel] <- "new level"
    }
  }
  out$inv_metric[is.na(out$inv_metric)] <- 1
  out$step_size <- src$step_size[1]

  n_default <- sum(out$source == "default")
  if (n_default) {
    message(n_default, " of ", nrow(out), " parameters have no counterpart in the ",
            "source and take Stan's default variance and a generic start; ",
            "see print() for which.")
  }
  out
}


# Starting values -------------------------------------------------------------

# The target's init plan with the source means written in wherever the table
# has them. A Cholesky factor is rebuilt from its strict lower triangle (rows
# of a correlation Cholesky factor have unit norm, so the diagonal follows);
# an invalid result keeps the identity default.
#' @keywords internal
.warmstart_apply_means <- function(plan, tab) {
  base <- sub("\\[.*$", "", tab$parameter)
  lapply(plan, function(e) {
    rows <- which(base == e$name)
    if (!length(rows)) return(e)
    means <- tab$mean[rows]
    if (all(is.na(means))) return(e)
    idx <- .label_indices(tab$parameter[rows])

    if (identical(e$kind, "fixed") && length(e$dim) == 2 && e$dim[1] == e$dim[2] &&
        length(rows) == e$dim[1] * (e$dim[1] - 1) / 2) {
      if (anyNA(means)) return(e)
      L <- diag(e$dim[1])
      for (r in seq_along(rows)) L[idx[[r]][1], idx[[r]][2]] <- means[r]
      dg <- 1 - rowSums(L^2) + diag(L)^2
      if (any(dg <= 0)) return(e)
      diag(L) <- sqrt(dg)
      e$value <- as.vector(L)
      return(e)
    }
    if (!identical(e$kind, "bounds") && !identical(e$kind, "sorted")) return(e)

    v <- e$value
    if (length(e$dim) > 1) dim(v) <- e$dim
    for (r in seq_along(rows)) {
      if (is.na(means[r])) next
      if (!length(e$dim)) v <- means[r]
      else if (length(e$dim) == 1) v[idx[[r]][1]] <- means[r]
      else v[matrix(idx[[r]], 1)] <- means[r]
    }
    v <- as.vector(v)
    # keep inside the bounds Stan will check
    if (!is.na(e$lower)) v <- pmax(v, e$lower + 1e-8)
    if (!is.na(e$upper)) v <- pmin(v, e$upper - 1e-8)
    if (identical(e$kind, "sorted")) v <- sort(v)
    e$value <- v
    e
  })
}
