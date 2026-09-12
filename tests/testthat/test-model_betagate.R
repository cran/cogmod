context("Beta-Gate")

test_that("rcogmod_betagate empirical proportions match theoretical probabilities", {
  skip_on_cran()
  n <- 20000 # Increased sample size for better stability
  set.seed(123)

  # Tolerances
  tol_prob <- 0.05 # Allow some sampling variation

  # Parameter grid, thinned to a covering subset - see helper-grid.R. The
  # degenerate corners are pinned in: pex = 0 (no gate at all) and pex = 1 with
  # bex at each end, where one point mass takes everything.
  grid <- covering_grid(
    mu = c(0.2, 0.5, 0.8),
    phi = c(1, 3, 10),
    pex = c(0, 0.1, 0.5, 0.9, 1),
    bex = c(0, 0.3, 0.5, 0.7, 1),
    always = function(g) {
      g$mu == 0.5 & g$phi == 3 & (g$pex == 0 | (g$pex == 1 & g$bex %in% c(0, 1)))
    }
  )

  for (.i in seq_len(nrow(grid))) {
    mu <- grid$mu[.i]
    phi <- grid$phi[.i]
    pex <- grid$pex[.i]
    bex <- grid$bex[.i]

    l <- paste0("rcogmod_betagate(mu=", mu, ", phi=", phi, ", pex=", pex, ", bex=", bex, ")")

    # --- Calculate Theoretical Probabilities (New Logic) ---
    # Cutpoints on probability scale
    cutzero <- pex * (1 - bex)
    cutone <- 1 - pex * bex

    # Handle edge cases for qlogis
    eps <- .Machine$double.eps
    cutzero <- pmax(eps, pmin(1 - eps, cutzero))
    cutone <- pmax(eps, pmin(1 - eps, cutone))
    cutone <- pmax(cutzero, cutone) # Ensure cutzero <= cutone

    # Cutpoints on logit scale
    cutzerolog <- qlogis(cutzero)
    cutonelog <- qlogis(cutone)

    # Location parameter on logit scale
    mu_ql <- qlogis(mu)

    # Theoretical probabilities using plogis
    theo_p0 <- plogis(cutzerolog, location = mu_ql, lower.tail = TRUE)
    theo_p1 <- plogis(cutonelog, location = mu_ql, lower.tail = FALSE)
    theo_pmid <- 1 - theo_p0 - theo_p1
    theo_pmid <- pmax(0, theo_pmid) # Handle floating point inaccuracies

    # Normalize (should be close to 1 already)
    total_p <- theo_p0 + theo_pmid + theo_p1
    if (abs(total_p - 1) > 1e-6) { # Avoid division by zero if total_p is 0
        theo_p0 <- theo_p0 / total_p
        theo_p1 <- theo_p1 / total_p
        theo_pmid <- theo_pmid / total_p
  }

  # Generate sample
  x_sample <- rcogmod_betagate(n, mu = mu, phi = phi, pex = pex, bex = bex)

  # Calculate empirical probabilities
  emp_p0 <- mean(x_sample == 0)
  emp_p1 <- mean(x_sample == 1)
  emp_pmid <- mean(x_sample > 0 & x_sample < 1)

  # --- Check Point Masses ---
  expect_equal(emp_p0, theo_p0, tolerance = tol_prob,
               label = paste(l, "- P(0): Empirical vs Theoretical (plogis)"))
  expect_equal(emp_p1, theo_p1, tolerance = tol_prob,
               label = paste(l, "- P(1): Empirical vs Theoretical (plogis)"))
  expect_equal(emp_pmid, theo_pmid, tolerance = tol_prob,
               label = paste(l, "- P(mid): Empirical vs Theoretical (plogis)"))

  }
})


test_that("dcogmod_betagate matches rcogmod_betagate empirical distribution and integrates correctly", {
  skip_on_cran()
  n <- 50000 # Use a larger sample size
  set.seed(456)

  # Tolerances
  tol_prob <- 0.05
  tol_integration <- 0.01 # Tighter integration tolerance
  tol_quantile <- 0.05 # Tolerance for quantile comparison

  # The same grid as above, thinned to a covering subset - see helper-grid.R. The
  # degenerate corners are pinned in: pex = 0 (no gate at all) and pex = 1 with
  # bex at each end, where one point mass takes everything.
  grid <- covering_grid(
    mu = c(0.2, 0.5, 0.8),
    phi = c(1, 3, 10),
    pex = c(0, 0.1, 0.5, 0.9, 1),
    bex = c(0, 0.3, 0.5, 0.7, 1),
    always = function(g) {
      g$mu == 0.5 & g$phi == 3 & (g$pex == 0 | (g$pex == 1 & g$bex %in% c(0, 1)))
    }
  )

  # Quantiles to check
  quantile_probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  for (.i in seq_len(nrow(grid))) {
    mu <- grid$mu[.i]
    phi <- grid$phi[.i]
    pex <- grid$pex[.i]
    bex <- grid$bex[.i]

    l <- paste0("dcogmod_betagate(mu=", mu, ", phi=", phi, ", pex=", pex, ", bex=", bex, ")")

    # Generate sample
    x_sample <- rcogmod_betagate(n, mu = mu, phi = phi, pex = pex, bex = bex)

    # Calculate empirical probabilities
    emp_p0 <- mean(x_sample == 0)
    emp_p1 <- mean(x_sample == 1)
    emp_pmid_count <- sum(x_sample > 0 & x_sample < 1)

    # Calculate theoretical probabilities/densities at extremes using dcogmod_betagate
    theo_p0_db <- dcogmod_betagate(0, mu = mu, phi = phi, pex = pex, bex = bex)
    theo_p1_db <- dcogmod_betagate(1, mu = mu, phi = phi, pex = pex, bex = bex)

    # --- Check Point Masses (dcogmod_betagate vs empirical) ---
    # This comparison remains valid as dcogmod_betagate uses the new internal logic
    expect_equal(emp_p0, theo_p0_db, tolerance = tol_prob,
                 label = paste(l, "- P(0): Empirical vs dcogmod_betagate"))
    expect_equal(emp_p1, theo_p1_db, tolerance = tol_prob,
                 label = paste(l, "- P(1): Empirical vs dcogmod_betagate"))

    # --- Integration Check ---
    # This check remains valid as it integrates the output of the updated dcogmod_betagate
    integrand <- function(x_int) {
      dcogmod_betagate(x_int, mu = mu, phi = phi, pex = pex, bex = bex)
  }

  integral_val <- 0
  # Only integrate if there's a continuous region (pex < 1 implies prob_mid > 0)
  # Calculate prob_mid theoretically to decide if integration is needed
  cutzero_int <- pex * (1 - bex)
  cutone_int <- 1 - pex * bex
  eps_int <- .Machine$double.eps
  cutzero_int <- pmax(eps_int, pmin(1 - eps_int, cutzero_int))
  cutone_int <- pmax(eps_int, pmin(1 - eps_int, cutone_int))
  cutone_int <- pmax(cutzero_int, cutone_int)
  cutzerolog_int <- qlogis(cutzero_int)
  cutonelog_int <- qlogis(cutone_int)
  mu_ql_int <- qlogis(mu)
  prob_0_int <- plogis(cutzerolog_int, location = mu_ql_int, lower.tail = TRUE)
  prob_1_int <- plogis(cutonelog_int, location = mu_ql_int, lower.tail = FALSE)
  prob_mid_int <- 1 - prob_0_int - prob_1_int
  prob_mid_int <- pmax(0, prob_mid_int)

  if (prob_mid_int > 1e-9) { # Check if middle probability is non-negligible
     integral_result <- tryCatch(
        stats::integrate(integrand, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps,
                         subdivisions = 500, rel.tol = 1e-4)$value, # Tighter rel.tol for integrate
        error = function(e) {
            warning("Integration failed for ", l, ": ", e$message)
            NA
        }
     )
     if (!is.na(integral_result)) {
         integral_val <- integral_result
     }
  }

  # Check if total probability (point masses + integral) is 1
  total_prob <- theo_p0_db + theo_p1_db + integral_val
  expect_equal(total_prob, 1, tolerance = tol_integration,
               label = paste(l, "- Theoretical density integrates to 1"))

  # --- Quantile Comparison for Continuous Part ---
  # Only run if prob_mid > 0 and there are enough middle values empirically
  if (prob_mid_int > 1e-9 && emp_pmid_count > 50) { # Increased required count
    x_middle <- x_sample[x_sample > 0 & x_sample < 1]

    # Calculate empirical quantiles from the middle part
    emp_quantiles <- stats::quantile(x_middle, probs = quantile_probs, type = 8, names = FALSE)

    # Calculate theoretical quantiles from the underlying Beta distribution
    shape1 <- mu * phi * 2
    shape2 <- (1 - mu) * phi * 2

    # Check for valid shapes before calling qbeta
    if (shape1 > 0 && shape2 > 0) {
        theo_quantiles <- tryCatch(
            stats::qbeta(quantile_probs, shape1 = shape1, shape2 = shape2),
            error = function(e) {
                warning("qbeta failed for ", l, ": ", e$message)
                rep(NA, length(quantile_probs))
            }
        )

        # Check if theoretical quantiles could be calculated
        if(any(is.na(theo_quantiles))) {
            warning("Skipping Quantile test for ", l, " due to NA in theoretical quantile calculation (qbeta).")
        } else {
            # Compare empirical and theoretical quantiles
            expect_equal(emp_quantiles, theo_quantiles, tolerance = tol_quantile,
                         label = paste(l, "- Quantiles (Empirical vs Theoretical Beta)"))
        }
    } else {
         warning("Skipping Quantile test for ", l, " due to non-positive Beta shapes.")
    }
  } else if (prob_mid_int <= 1e-9) {
      # If theoretically no middle part, ensure empirically there isn't much either
      expect_lt(emp_pmid_count / n, tol_prob,
                label = paste(l, "- Empirical middle count low when theo_pmid is zero"))
  }
  }
})


context("Beta-Gate - brms")

test_that("Beta-Gate model can recover parameters with brms using variational inference", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan") # Needed by brms for custom families sometimes

  # --- 1. Set up Simulation ---
  true_mu <- 0.7
  true_phi <- 5
  true_pex <- 0.2
  true_bex <- 0.6
  n_obs <- 3000 # Increased sample size for better stability

  # Generate synthetic data with known parameters
  set.seed(456)
  df <- data.frame(
    y = rcogmod_betagate(n = n_obs, mu = true_mu, phi = true_phi,
                  pex = true_pex, bex = true_bex)
  )

  # --- 2. Define and Fit brms Model ---
  f <- brms::bf(
    y ~ 1,
    phi ~ 1,
    pex ~ 1,
    bex ~ 1,
    family = cogmod_betagate()
  )

  m <- brms::brm(
    formula = f,
    data = df,
    init = 0,
    stanvars = cogmod_betagate_stanvars(),
    backend = "cmdstanr",
    seed = 456,
    refresh = 0,
    algorithm = "pathfinder"
  )


  # --- 3. Check Parameter Recovery (point estimate vs true value) ---
  post_summary <- brms::posterior_summary(m)
  # Mu
  mu_est <- brms::inv_logit_scaled(post_summary["b_Intercept", "Estimate"])
  expect_equal(mu_est, true_mu, tolerance = 0.15, label = "Recovered mu")

  # Phi (softplus link)
  phi_est <- log(1 + exp(post_summary["b_phi_Intercept", "Estimate"]))
  expect_equal(phi_est, true_phi, tolerance = 0.5, label = "Recovered phi")

  # Pex
  pex_est <- brms::inv_logit_scaled(post_summary["b_pex_Intercept", "Estimate"])
  expect_equal(pex_est, true_pex, tolerance = 0.15, label = "Recovered pex")

  # Bex
  bex_est <- brms::inv_logit_scaled(post_summary["b_bex_Intercept", "Estimate"])
  expect_equal(bex_est, true_bex, tolerance = 0.15, label = "Recovered bex")

  # --- 4. Test Post-processing Functions ---
  n_pred_draws <- 10
  pred <- brms::posterior_predict(m, ndraws = n_pred_draws)
  expect_equal(nrow(pred), n_pred_draws)
  expect_equal(ncol(pred), n_obs)
  expect_true(all(pred >= 0 & pred <= 1), "Posterior predictions outside [0, 1]")
  expect_false(any(is.na(pred)), "NA values in posterior predictions")

  n_newdata <- 5
  pred_new <- brms::posterior_predict(m, ndraws = n_pred_draws, newdata = df[1:n_newdata, ])
  expect_equal(nrow(pred_new), n_pred_draws)
  expect_equal(ncol(pred_new), n_newdata)
  expect_true(all(pred_new >= 0 & pred_new <= 1), "Posterior predictions (newdata) outside [0, 1]")
  expect_false(any(is.na(pred_new)), "NA values in posterior predictions (newdata)")

  n_ll_draws <- 5
  ll <- brms::log_lik(m, ndraws = n_ll_draws)
  expect_equal(nrow(ll), n_ll_draws)
  expect_equal(ncol(ll), n_obs)
  expect_true(all(is.finite(ll)), "Non-finite values found in log-likelihood")

  n_epred_draws <- 5
  epred <- brms::posterior_epred(m, ndraws = n_epred_draws)
  expect_equal(nrow(epred), n_epred_draws)
  expect_equal(ncol(epred), n_obs)
  expect_true(all(epred >= 0 & epred <= 1, na.rm = TRUE), "Posterior epred outside [0, 1]")
  expect_false(any(is.na(epred)), "NA values in posterior epred")
})


test_that("Stan cogmod_betagate_lpdf matches R dcogmod_betagate function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function if possible
  cogmod_betagate_lpdf <- stan_fun("cogmod_betagate")

  # --- Define parameter grids for testing ---
  y_values <- c(0, 0.01, 0.5, 0.99, 1) # Test boundaries and mid-points
  mu_values <- c(0.2, 0.5, 0.8)
  phi_values <- c(1, 5, 10) # Include phi=1 (uniform beta case)
  pex_values <- c(0, 0.3, 0.7, 1) # Include boundary cases for pex
  bex_values <- c(0, 0.5, 1)   # Include boundary cases for bex

  # --- Loop through parameter combinations ---
  for (mu in mu_values) {
    for (phi in phi_values) {
      for (pex in pex_values) {
        for (bex in bex_values) {
          # Test across different y values for the current parameter set
          for (y in y_values) {
            label <- sprintf("y=%.2f, mu=%.1f, phi=%.1f, pex=%.1f, bex=%.1f",
                             y, mu, phi, pex, bex)

            # Calculate log-density using Stan function
            stan_log_lik <- cogmod_betagate_lpdf(y, mu, phi, pex, bex)

            # Calculate log-density using R function
            r_log_lik <- dcogmod_betagate(y, mu, phi, pex, bex, log = TRUE)

            # Compare log-likelihoods (handles -Inf correctly)
            expect_equal(stan_log_lik, r_log_lik, tolerance = 1e-6,
                         label = paste("Log-likelihood comparison:", label))
          }
        }
      }
    }
  }

  # --- Test invalid parameter handling ---
  # R function should warn and return NaN/-Inf
  expect_error(dcogmod_betagate(0.5, mu = -0.1, phi = 5, pex = 0.1, bex = 0.5), "mu must be strictly between 0 and 1.")
  expect_error(dcogmod_betagate(0.5, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), "phi must be positive.")
  expect_error(dcogmod_betagate(0.5, mu = 0.5, phi = 5, pex = 1.1, bex = 0.5), "pex must be between 0 and 1.")
  expect_error(dcogmod_betagate(0.5, mu = 0.5, phi = 5, pex = 0.1, bex = -0.2), "bex must be between 0 and 1.")

})
