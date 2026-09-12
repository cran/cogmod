context("Ordered Beta")

# `.rordbeta()`, `.dordbeta()` and `.ordbeta_lpdf()` are internal and
# unexported - there is no `cogmod_ordbeta()` family - so nothing in the package
# exercises them. These tests pin the behaviour they would be promoted with.

test_that(".dordbeta splits its mass into two point masses and a Beta", {
  # The three components are P(y = 0), P(y = 1) and the continuous middle, and
  # by construction they have to sum to one at every parameter value.
  grid <- expand.grid(
    mu = c(0.3, 0.5, 0.7),
    phi = c(4, 10),
    cutzero = c(-2, -1),
    cutone = c(1, 2)
  )

  for (.i in seq_len(nrow(grid))) {
    mu <- grid$mu[.i]
    phi <- grid$phi[.i]
    cutzero <- grid$cutzero[.i]
    cutone <- grid$cutone[.i]

    p0 <- cogmod:::.dordbeta(0, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)
    p1 <- cogmod:::.dordbeta(1, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)
    middle <- stats::integrate(
      function(z) {
        cogmod:::.dordbeta(z, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)
      },
      lower = 0, upper = 1
    )$value

    # mu * phi and (1 - mu) * phi both exceed 1 across the grid, so the Beta
    # density is bounded and integrate() has nothing to trip over.
    expect_equal(p0 + middle + p1, 1, tolerance = 1e-4)
  }
})


test_that(".dordbeta point masses are the cutpoint probabilities", {
  mu <- 0.6
  cutzero <- -1
  cutone <- 1.5
  mu_ql <- stats::qlogis(mu)

  expect_equal(
    cogmod:::.dordbeta(0, mu = mu, phi = 5, cutzero = cutzero, cutone = cutone),
    1 - stats::plogis(mu_ql - cutzero)
  )
  expect_equal(
    cogmod:::.dordbeta(1, mu = mu, phi = 5, cutzero = cutzero, cutone = cutone),
    stats::plogis(mu_ql - cutone)
  )
})


test_that(".dordbeta honours `log` and vectorises over mu", {
  d <- cogmod:::.dordbeta(c(0, 0.5, 1), mu = 0.6, phi = 5)
  ld <- cogmod:::.dordbeta(c(0, 0.5, 1), mu = 0.6, phi = 5, log = TRUE)
  expect_equal(ld, log(d))

  x <- c(0.2, 0.5, 0.8)
  mu <- c(0.3, 0.5, 0.7)
  vectorised <- cogmod:::.dordbeta(x, mu = mu, phi = 5)
  one_at_a_time <- vapply(
    seq_along(x),
    function(i) cogmod:::.dordbeta(x[i], mu = mu[i], phi = 5),
    numeric(1)
  )
  expect_equal(vectorised, one_at_a_time)
})


test_that(".rordbeta draws match the density's point masses", {
  skip_on_cran()
  n <- 20000
  mu <- 0.5
  phi <- 5
  cutzero <- -1
  cutone <- 1

  set.seed(123)
  x <- cogmod:::.rordbeta(n, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)

  expect_length(x, n)
  expect_true(all(x >= 0 & x <= 1))
  # Both boundaries are genuinely attained rather than merely approached - that
  # is the whole point of the ordered Beta over a plain Beta.
  expect_true(any(x == 0))
  expect_true(any(x == 1))

  p0 <- cogmod:::.dordbeta(0, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)
  p1 <- cogmod:::.dordbeta(1, mu = mu, phi = phi, cutzero = cutzero, cutone = cutone)

  expect_equal(mean(x == 0), p0, tolerance = 0.1)
  expect_equal(mean(x == 1), p1, tolerance = 0.1)
  expect_equal(mean(x > 0 & x < 1), 1 - p0 - p1, tolerance = 0.1)
})


test_that(".rordbeta and .dordbeta reject out-of-range parameters", {
  expect_error(cogmod:::.rordbeta(10, mu = 0), "between 0 and 1")
  expect_error(cogmod:::.rordbeta(10, mu = 1), "between 0 and 1")
  expect_error(cogmod:::.rordbeta(10, phi = 0), "greater than 0")
  expect_error(cogmod:::.rordbeta(10, mu = c(0.3, 0.4)), "length 1 or length N")

  expect_error(cogmod:::.dordbeta(c(0.2, 0.5), mu = 0), "between 0 and 1")
  expect_error(cogmod:::.dordbeta(c(0.2, 0.5), phi = -1), "greater than 0")
  expect_error(
    cogmod:::.dordbeta(rep(0.5, 3), mu = c(0.2, 0.3)),
    "length 1 or length N"
  )
})


test_that(".ordbeta_lpdf returns the Stan function it claims to", {
  code <- cogmod:::.ordbeta_lpdf()
  expect_type(code, "character")
  expect_length(code, 1)
  expect_match(code, "real ord_beta_reg_lpdf\\(", fixed = FALSE)
})
