#' The Stan code a cogmod family needs, read off the model
#'
#' @description
#' Returns the `stanvars` argument for `brms::brm()`, for whichever `cogmod`
#' family the model uses. It is a front end to the per-family
#' `<family>_stanvars()` functions - `cogmod_lognormal_stanvars()`,
#' `cogmod_choco_stanvars()`, `cogmod_ddm_stanvars()` and the rest - which remain available
#' and unchanged.
#'
#' @details
#' `brms` needs the custom likelihood injected into the generated Stan program,
#' and every `cogmod` family ships one. Calling this instead of the family's own
#' function means the family is named once, in `bf()`, rather than twice:
#'
#' ```r
#' f <- brms::bf(RT ~ Condition, ndt ~ Condition,
#'               family = cogmod_lognormal())
#'
#' brms::brm(f, data = df,
#'           prior    = cogmod_priors(f, df),
#'           init     = cogmod_inits(f, df),
#'           stanvars = cogmod_stanvars(f))
#' ```
#'
#' # A warning it may emit
#'
#' [cogmod_lba1()] and [cogmod_lba2()] have a likelihood that is *exactly*
#' constant along the ray that multiplies the drift rates, their SDs, the
#' start-point range and the threshold offset by a common factor. If the
#' formula pins none of them to a constant, this warns: the RT distribution is
#' still identified, but the individual parameters are not, and the fit will
#' converge to whatever the priors say about that direction rather than fail.
#' Fixing any one member in `bf()` - conventionally `sigmazero = 1` - silences
#' it. Leaving a parameter *out* of `bf()` does not count: `brms` estimates it
#' anyway.
#'
#' # What it accepts
#'
#' A `brms::bf()` formula carrying the family, the family object itself, or a
#' fitted `brmsfit` (useful for recompiling or for `update()`).
#'
#' @param formula A `brms::bf()` formula carrying the family, a `cogmod` family
#'   object, or a fitted `brmsfit`.
#' @param ... Passed to the family's own `<family>_stanvars()` function.
#'
#' @return A `stanvars` object, to pass to `brms::brm(stanvars = )`.
#'
#' @seealso [cogmod_priors()], [cogmod_inits()]
#'
#' @examples
#' f <- brms::bf(RT ~ 1, ndt ~ 1, family = cogmod_lognormal())
#' cogmod_stanvars(f)
#'
#' # Equivalent to naming the family a second time:
#' cogmod_lognormal_stanvars()
#'
#' @export
cogmod_stanvars <- function(formula, ...) {
  # The LBA families are identified only up to a common scale factor, and the
  # symptom of leaving it free is a fit that converges to meaningless parameter
  # values rather than one that fails. This is the generic a fit cannot skip,
  # so the warning goes here rather than in cogmod_priors() - a prior makes the
  # posterior proper, but it does not identify a direction the likelihood
  # cannot see.
  .warn_scale_ray(formula)
  family <- .cogmod_family(formula)
  fam <- .family_name(family)
  if (!isTRUE(fam %in% .cogmod_families())) {
    stop(
      "cogmod_stanvars() has no Stan code for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'. Supported: ", paste(.cogmod_families(), collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_lognormal()).",
      call. = FALSE
    )
  }

  fn <- get(paste0(fam, "_stanvars"), envir = asNamespace("cogmod"),
            mode = "function")
  do.call(fn, list(...))
}


# The families that ship Stan code. Built from the shifted-family registry so
# the two cannot drift apart.
#' @keywords internal
.cogmod_families <- function() {
  # unique() because every choice family is now on the outlier mixture too, so
  # they appear in .OUTLIER_FAMILIES as well as being named here.
  unique(c(
    "cogmod_betadiscrete", "cogmod_betagate", "cogmod_choco", "cogmod_ddm",
    "cogmod_lba2", "cogmod_lnr", "cogmod_rdm", "cogmod_exgaussian",
    "cogmod_geg",
    .OUTLIER_FAMILIES
  ))
}
