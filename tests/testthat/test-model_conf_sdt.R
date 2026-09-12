# ---------------------------------------------------------------------------
# EXPERIMENTAL / NOT READY: the confidence SDT model is commented out until it
# is finished. Uncomment this file (and the matching entries in NAMESPACE and
# _pkgdown.yml) to re-enable it.
# ---------------------------------------------------------------------------

# test_that("rconf_sdt runs and returns correct structure", {
#   n <- 100
#   dprime <- 1.5
#   c <- 0.1
#   thetazero <- c(-0.5, -1.5)
#   thetaone <- c(0.5, 1.5)
#
#   sim_data <- rconf_sdt(n, dprime, c, thetazero, thetaone)
#
#   expect_s3_class(sim_data, "data.frame")
#   expect_equal(nrow(sim_data), n)
#   expect_equal(names(sim_data), c("truth", "response", "confidence"))
#   expect_true(all(sim_data$truth %in% c(0, 1)))
#   expect_true(all(sim_data$response %in% c(0, 1)))
#   expect_true(is.integer(sim_data$confidence) || all(sim_data$confidence == floor(sim_data$confidence)))
# })
#
# test_that("rconf_sdt respects parameters", {
#   # With very high dprime, response should almost always equal truth
#   n <- 1000
#   dprime <- 10
#   c <- 0
#   thetazero <- c(-0.5, -1.5)
#   thetaone <- c(0.5, 1.5)
#   truth <- sample(c(0, 1), n, replace = TRUE)
#
#   sim_data <- rconf_sdt(n, dprime, c, thetazero, thetaone, truth = truth)
#   expect_gt(mean(sim_data$response == sim_data$truth), 0.99)
#
#   # With high positive c, bias towards response 1
#   dprime <- 1
#   c <- 5
#   # Adjust thetazero to be on the correct side of the new criterion (-5)
#   thetazero_biased <- c(-5.5, -6.5)
#   sim_data <- rconf_sdt(n, dprime, c, thetazero_biased, thetaone, truth = truth)
#   expect_gt(mean(sim_data$response == 1), 0.99) # Positive c -> bias to 1
#
#   # With high negative c, bias towards response 0
#   c <- -5
#   # Adjust thetaone to be on the correct side of the new criterion (5)
#   thetaone_biased <- c(5.5, 6.5)
#   sim_data <- rconf_sdt(n, dprime, c, thetazero, thetaone_biased, truth = truth)
#   expect_gt(mean(sim_data$response == 0), 0.99) # Negative c -> bias to 0
# })
#
#
# test_that("dconf_sdt returns valid probabilities", {
#   dprime <- 1.5
#   c <- 0.1
#   thetazero <- c(-0.5, -1.5)
#   thetaone <- c(0.5, 1.5)
#
#   # Test for one case
#   p <- dconf_sdt(truth = 1, response = 1, confidence = 2, dprime, c, thetazero, thetaone)
#   expect_type(p, "double")
#   expect_gte(p, 0)
#   expect_lte(p, 1)
#
#   # Test log probability
#   log_p <- dconf_sdt(truth = 1, response = 1, confidence = 2, dprime, c, thetazero, thetaone, log = TRUE)
#   expect_equal(log_p, log(p))
# })
#
# test_that("dconf_sdt probabilities sum to 1", {
#   dprime <- 1.5
#   c <- 0.1
#   n_conf_levels <- 3 # 2 thresholds -> 3 levels (0, 1, 2)
#   thetazero <- c(-0.5, -1.5)
#   thetaone <- c(0.5, 1.5)
#
#   # For truth = 1
#   total_p_truth1 <- 0
#   for (response in c(0, 1)) {
#     for (confidence in 0:(n_conf_levels - 1)) {
#       total_p_truth1 <- total_p_truth1 + dconf_sdt(1, response, confidence, dprime, c, thetazero, thetaone)
#     }
#   }
#   expect_equal(total_p_truth1, 1)
#
#   # For truth = 0
#   total_p_truth0 <- 0
#   for (response in c(0, 1)) {
#     for (confidence in 0:(n_conf_levels - 1)) {
#       total_p_truth0 <- total_p_truth0 + dconf_sdt(0, response, confidence, dprime, c, thetazero, thetaone)
#     }
#   }
#   expect_equal(total_p_truth0, 1)
# })
#
# test_that("conf_sdt functions handle parameter validation", {
#   # rconf_sdt
#   expect_error(rconf_sdt(10, 1, c(0, 1), c(-0.5), c(0.5)), "'c' must be a single value")
#   expect_error(rconf_sdt(10, 1, 0, c(-0.5), c(-0.5)), "All thetaone must be > -c")
#   expect_error(rconf_sdt(10, 1, 0, c(0.5), c(0.5)), "All thetazero must be < -c")
#   expect_error(rconf_sdt(10, 1, 0, c(-0.5), c(1.5, 0.5)), "thetaone must be sorted")
#   expect_error(rconf_sdt(10, 1, 0, c(-1.5, -0.5), c(0.5)), "thetazero must be sorted in decreasing order")
#
#   # dconf_sdt
#   expect_error(dconf_sdt(1, 1, 1, 1, c(0, 1), c(-0.5), c(0.5)), "'c' must be a single value")
#   expect_error(dconf_sdt(1, 1, 2, 1.5, 0.1, c(-0.5), c(0.5)), "A confidence value for response=1 is higher than allowed")
#   expect_error(dconf_sdt(1, 0, 2, 1.5, 0.1, c(-0.5), c(0.5)), "A confidence value for response=0 is higher than allowed")
# })
#
#
# context("conf_sdt - brms")
#
# test_that("conf_sdt model can recover parameters with brms", {
#   # Skip on CRAN and when not running full tests
#   skip_on_cran()
#   skip_if_not_installed("brms")
#   skip_if_not_installed("cmdstanr")
#
#   # Set random seed for reproducibility
#   set.seed(12345)
#
#   # --- 1. Generate synthetic data with known parameters ---
#   n <- 2500
#   true_dprime <- 1.5
#   true_c <- 0.2
#   # Absolute criteria for data generation
#   true_thetaone <- c(0.5, 1.5)
#   true_thetazero <- c(-0.6, -1.8)
#
#   # Calculate the true offsets that the model will estimate
#   decision_criterion <- -true_c
#   true_toneone_off <- true_thetaone[1] - decision_criterion
#   true_tonetwo_off <- true_thetaone[2] - true_thetaone[1]
#   true_tzeroone_off <- decision_criterion - true_thetazero[1]
#   true_tzerotwo_off <- true_thetazero[1] - true_thetazero[2]
#
#   # Generate data using the absolute criteria
#   df <- rconf_sdt(n, dprime = true_dprime, c = true_c,
#                   thetazero = true_thetazero, thetaone = true_thetaone)
#
#   # --- 2. Define brms formula and priors for the offset parameters ---
#   # The `vint()` function passes variables to the custom Stan code.
#   # The main formula `confidence ~ 1` models `dprime`.
#   f <- brms::bf(
#     confidence | dec(response) + vint(truth) ~ 1,
#     c ~ 1,
#     tzeroone ~ 1,
#     tzerotwo ~ 1,
#     toneone ~ 1,
#     tonetwo ~ 1
#   )
#
#   # Priors on the true parameter values (offsets for criteria)
#   # The log link in the family definition handles the positivity constraint.
#   # We use lognormal priors as they are appropriate for positive-only parameters.
#   priors <- c(
#     brms::prior(normal(1.5, 0.1), class = "Intercept"), # for dprime
#     brms::prior(normal(0.2, 0.1), class = "Intercept", dpar = "c"),
#     brms::prior(normal(log(0.7), 0.1), class = "Intercept", dpar = "toneone"),
#     brms::prior(normal(log(1.0), 0.1), class = "Intercept", dpar = "tonetwo"),
#     brms::prior(normal(log(0.4), 0.1), class = "Intercept", dpar = "tzeroone"),
#     brms::prior(normal(log(1.2), 0.1), class = "Intercept", dpar = "tzerotwo")
#   ) |>
#     brms::validate_prior(f, data = df, family = conf_sdt_custom_family())
#
#   # --- 3. Fit with brms ---
#   fit <- brms::brm(
#     formula = f,
#     data = df,
#     family = conf_sdt_custom_family(),
#     stanvars = conf_sdt_stanvars(),
#     prior = priors,
#     init = 0,
#     backend = "cmdstanr",
#     algorithm = "sampling",
#     # refresh = 0,
#     cores = 4 # Use multiple cores
#   )
#
#   # --- 4. Check parameter recovery ---
#   post <- brms::posterior_summary(fit)
#   means <- post[, "Estimate"]
#   tolerance <- 0.3
#
#   # dprime (Intercept of the main formula)
#   expect_equal(means[["b_Intercept"]], true_dprime, tolerance = tolerance,
#                label = "dprime recovery")
#
#   # c
#   expect_equal(means[["b_c_Intercept"]], true_c, tolerance = tolerance,
#                label = "c recovery")
#
#   # toneone (offset)
#   expect_equal(means[["b_toneone_Intercept"]], true_toneone_off, tolerance = tolerance,
#                label = "toneone offset recovery")
#
#   # tonetwo (offset)
#   expect_equal(means[["b_tonetwo_Intercept"]], true_tonetwo_off, tolerance = tolerance,
#                label = "tonetwo offset recovery")
#
#   # tzeroone (offset)
#   expect_equal(means[["b_tzeroone_Intercept"]], true_tzeroone_off, tolerance = tolerance,
#                label = "tzeroone offset recovery")
#
#   # tzerotwo (offset)
#   expect_equal(means[["b_tzerotwo_Intercept"]], true_tzerotwo_off, tolerance = tolerance,
#                label = "tzerotwo offset recovery")
#
#   # --- 5. Test Post-processing Functions ---
#   # Test posterior prediction
#   n_pred_draws <- 10
#   # Ensure newdata has the 'dec' column
#   newdata_pred <- df[1:5, ]
#   pred <- brms::posterior_predict(fit, ndraws = n_pred_draws, newdata = newdata_pred)
#   expect_true(is.matrix(pred))
#   expect_equal(nrow(pred), n_pred_draws)
#   expect_equal(ncol(pred), 5)
#   expect_true(all(pred %in% 0:2, na.rm = TRUE), "Predicted confidence should be 0, 1, or 2")
#
#   # Test log-likelihood
#   ll <- brms::log_lik(fit, ndraws = 5)
#   expect_true(is.matrix(ll))
#   expect_equal(nrow(ll), 5)
#   expect_equal(ncol(ll), n)
#   expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")
# })
