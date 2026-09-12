# Opt-in gate for the model-fitting tests
# =======================================
#
# Nine blocks in the suite fit a real model with `brms::brm()`. Each compiles its
# own Stan program, which on Windows costs more than every other test in the
# file put together, and they are the reason a full run takes the better part of
# an hour. They are worth having - they are the only end-to-end check that a
# family's Stan code, priors, inits and post-processing agree - but they are not
# worth paying for on every `devtools::test()`.
#
# `skip_if_not_slow()` runs them only when `COGMOD_TEST_SLOW` is set to a true
# value. Set it for one run:
#
#   Sys.setenv(COGMOD_TEST_SLOW = "true"); devtools::test()
#
# or export it in the environment before starting R.
#
# Note what this does *not* change: these tests also need cmdstanr, and the
# R-CMD-check workflow deliberately uninstalls it (see the "Remove cmdstanr"
# step), so they have always been local-only. `COGMOD_TEST_SLOW` is set there
# too, so that a CI job which does provide CmdStan runs them without any further
# change - but on the current workflow they still skip for want of cmdstanr, one
# line earlier than before.
skip_if_not_slow <- function() {
  # as.logical() alone would make COGMOD_TEST_SLOW=1 mean "no", which is not
  # what anyone typing it intends.
  on <- tolower(trimws(Sys.getenv("COGMOD_TEST_SLOW", "")))
  if (!on %in% c("true", "t", "yes", "y", "1", "on")) {
    testthat::skip("model-fitting test; set COGMOD_TEST_SLOW=true to run")
  }
  invisible(TRUE)
}
