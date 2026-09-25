# EXPERIMENTAL: Rating Process Model (RAPROC)
#
# Generates continuous rating responses from a ballistic accumulation process.
# On each trial, evidence accumulates vertically from a starting bias, with
# trial-specific horizontal drift determining the angle of the trajectory.
# The response is the horizontal position where the trajectory reaches the
# fixed boundary, constrained to the [0, 1] rating scale.
#
# @param n Number of responses to generate.
# @param drift Mean horizontal displacement of the accumulation trajectory
#   per unit of vertical accumulation. A value of 0 corresponds to a
#   vertical trajectory, with positive values indicating drift towards
#   higher ratings and negative values indicating drift towards lower
#   ratings. The corresponding trajectory angle is `atan(drift)`, such
#   that `drift = -1` corresponds to an angle of -45 degrees from the
#   vertical.
# @param sigmadrift Standard deviation of trial-to-trial variability in
#   horizontal drift. Larger values indicate greater variability in the
#   direction of the accumulation trajectory.
# @param bias Starting position on the rating scale, between 0 and 1.
# @inheritParams rcogmod_betagate
#
# @return A numeric vector of simulated rating responses in [0, 1].
# rcogmod_raproc <- function(n, drift = 0, sigmadrift = 0.1, bias = 0.5) {
#   height <- 1 # Threshold height, fixed for identifiability
#
#   slope <- rnorm(n, mean = drift, sd = sigmadrift)
#   x <- bias + height * slope
#
#   pmin(pmax(x, 0), 1)
# }

# @rdname rcogmod_raproc
# dcogmod_raproc <- function(x,
#                            drift = 0,
#                            sigmadrift = 0.1,
#                            bias = 0.5,
#                            log = FALSE) {
#   if (any(sigmadrift <= 0)) {
#     stop("sigmadrift must be > 0.")
#   }
#
#   # Standardised boundaries in drift space
#   z0 <- (0 - bias - drift) / sigmadrift
#   z1 <- (1 - bias - drift) / sigmadrift
#
#   # Allocate output
#   out <- numeric(length(x))
#
#   # Lower boundary: P(X = 0)
#   i0 <- x == 0
#   out[i0] <- pnorm(z0)
#
#   # Upper boundary: P(X = 1)
#   i1 <- x == 1
#   out[i1] <- 1 - pnorm(z1)
#
#   # Continuous interior
#   interior <- x > 0 & x < 1
#   out[interior] <- dnorm(
#     x[interior],
#     mean = bias + drift,
#     sd = sigmadrift
#   )
#
#   # Outside support
#   outside <- x < 0 | x > 1 | is.na(x)
#   out[outside & !is.na(x)] <- 0
#   out[is.na(x)] <- NA_real_
#
#   if (log) {
#     out <- log(out)
#   }
#
#   out
# }
