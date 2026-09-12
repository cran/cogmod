# A covering subset of a full factorial
# =====================================
#
# Several blocks in the suite sweep a parameter grid. Written as nested loops
# over seven factors they reach 192 combinations, each drawing tens of thousands
# of variates and, in one case, running two numerical integrations - minutes of
# work per block.
#
# Nothing under test couples those parameters: the samplers and densities take
# each one through the same arithmetic whatever the others are, so the only
# thing a full factorial buys over a design that covers every *level* of every
# factor is running time. `covering_grid()` returns such a subset. Every level
# of every column appears at least once, rows matching `always` are kept
# whatever else happens - that is where the cases a comment singles out go, such
# as a drift of exactly zero - and the choice is deterministic, so a failure is
# reproducible from the seed alone.
#
# If a genuine interaction is ever suspected, name it in `always` rather than
# going back to the factorial.
covering_grid <- function(..., always = NULL, seed = 1) {
  g <- expand.grid(..., KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  keep <- if (is.null(always)) rep(FALSE, nrow(g)) else as.logical(always(g))

  levels_of <- lapply(g, unique)
  covered <- function() {
    rows <- which(keep)
    all(vapply(
      seq_along(g),
      function(j) all(levels_of[[j]] %in% g[rows, j]),
      logical(1)
    ))
  }

  # Own RNG stream, restored afterwards, so that thinning the grid cannot move
  # the draws the calling test makes from its own set.seed().
  old_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv)
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, .GlobalEnv),
          add = TRUE)
  set.seed(seed)
  for (i in sample(nrow(g))) {
    if (covered()) break
    keep[i] <- TRUE
  }
  g[keep, , drop = FALSE]
}
