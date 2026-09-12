# ---------------------------------------------------------------------------
# EXPERIMENTAL / NOT READY: the confidence SDT model is commented out until it
# is finished. Uncomment this file (and the matching entries in NAMESPACE and
# _pkgdown.yml) to re-enable it.
# ---------------------------------------------------------------------------

# #' @title Confidence Signal Detection Theory Model (EXPERIMENTAL)
# #'
# #' @description
# #' Functions for the Signal Detection Theory (SDT) model of confidence. This model
# #' assumes that a single sample of sensory evidence is used to generate both a
# #' discrimination response and a confidence rating.
# #'
# #' Functions:
# #' - `rconf_sdt()`: Simulates random draws from the SDT model.
# #' - `dconf_sdt()`: Computes the likelihood/density of observed responses and ratings.
# #'
# #' @param n Number of simulated trials.
# #' @param dprime Sensitivity parameter(s) (d'). Can be a single value or a vector.
# #' @param c Response bias parameter. Must be a single value.
# #' @param thetazero A numeric vector of confidence criteria for response=0. Must be sorted in decreasing order.
# #' @param thetaone A numeric vector of confidence criteria for response=1. Must be sorted in increasing order.
# #' @param truth Stimulus value (0 or 1). For simulation, can be a vector of length `n`. If NULL, it's sampled. For density, it's the observed stimulus.
# #' @param response Observed response (0 or 1). For density calculation.
# #' @param confidence Observed confidence rating (integer >= 0). For density calculation.
# #' @param log Logical; if TRUE, returns the log-density.
# #'
# #' @details
# #' The SDT model assumes that for a given stimulus `truth` (0 or 1), sensory
# #' evidence `x` is drawn from a normal distribution. When `truth` is 1, the
# #' evidence is drawn from `N(dprime / 2, 1)`; when `truth` is 0, it's drawn
# #' from `N(-dprime / 2, 1)`.
# #' The `response` is 1 if `x > -c` and 0 otherwise. The decision criterion is thus `-c`.
# #' The confidence rating is determined by which region the evidence `x` falls into,
# #' defined by the confidence criteria `thetaone` and `thetazero`.
# #' A rating of 0 corresponds to the lowest confidence (i.e., evidence closest to the decision criterion).
# #'
# #' **Parameter Interpretation:**
# #' - `dprime`: The sensitivity (d'). It reflects the observer's ability to distinguish between the two stimulus categories. Higher values indicate better performance.
# #' - `c`: The response bias. It reflects the observer's tendency to favor one response over the other, independent of the stimulus. A positive value indicates a bias towards response 1, and a negative value indicates a bias towards response 0. A value of 0 is unbiased.
# #' - `thetazero` & `thetaone`: The confidence criteria. These thresholds partition the evidence space into different confidence levels for each response.
# #'   - For `response = 1`, `thetaone` are the criteria. They must be greater than the decision criterion (`-c`). For a given response, higher confidence ratings correspond to values of `x` further away from the decision criterion. For example, `confidence = 0` if `-c < x < thetaone[1]`, `confidence = 1` if `thetaone[1] < x < thetaone[2]`, and so on.
# #'   - For `response = 0`, `thetazero` are the criteria. They must be less than the decision criterion (`-c`). For example, `confidence = 0` if `thetazero[1] < x < -c`, `confidence = 1` if `thetazero[2] < x < thetazero[1]`, and so on.
# #'
# #' @references
# #' - Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.
# #'
# #' @examples
# #' # Simulate data with 2 criteria per response (3 confidence levels: 0, 1, 2)
# #' sim_data <- rconf_sdt(
# #'   n = 1000, dprime = 1, c = 0.2,
# #'   thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
# #' )
# #' table(sim_data$response, sim_data$confidence)
# #'
# #' # Calculate density for confidence=2 (highest) for response=1
# #' dconf_sdt(
# #'   truth = 1, response = 1, confidence = 2, dprime = 1, c = 0.2,
# #'   thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
# #' )
# #'
# #' @name conf_sdt
# NULL
#
# #' @rdname conf_sdt
# #' @export
# rconf_sdt <- function(n, dprime, c, thetazero, thetaone, truth = NULL) {
#   # --- Prepare and Validate Parameters ---
#   if (is.null(truth)) {
#     truth <- sample(c(0, 1), n, replace = TRUE)
#   }
#   params <- .prepare_conf_sdt(n = n, dprime = dprime, c = c, truth = truth)
#   if (length(c) > 1) stop("'c' must be a single value.")
#
#   decision_criterion <- -c
#
#   if (any(thetaone <= decision_criterion)) stop("All thetaone must be > -c.")
#   if (any(thetazero >= decision_criterion)) stop("All thetazero must be < -c.")
#   if (is.unsorted(thetaone)) stop("thetaone must be sorted in increasing order.")
#   if (is.unsorted(rev(thetazero))) stop("thetazero must be sorted in decreasing order.")
#
#   # --- Simulation ---
#   # Generate sensory evidence
#   mean_x <- ((params$truth * 2) - 1) * params$dprime / 2
#   x <- stats::rnorm(params$ndraws, mean = mean_x, sd = 1)
#
#   # Generate response (0 or 1)
#   response <- ifelse(x > decision_criterion, 1, 0)
#
#   # Generate confidence rating (0-indexed)
#   confidence <- numeric(params$ndraws)
#
#   idx_r1 <- which(response == 1)
#   if (length(idx_r1) > 0) {
#     # Confidence = 0 for -c < x < t1, 1 for t1 < x < t2, etc.
#     # findInterval is 1-based, so we subtract 1 to get 0-based ratings.
#     confidence[idx_r1] <- findInterval(x[idx_r1], vec = c(decision_criterion, thetaone)) - 1
#   }
#
#   idx_r0 <- which(response == 0)
#   if (length(idx_r0) > 0) {
#     # Confidence = 0 for tz1 < x < -c, 1 for tz2 < x < tz1, and so on.
#     # This is equivalent to counting how many criteria in thetazero x is smaller than.
#     confidence[idx_r0] <- rowSums(outer(x[idx_r0], thetazero, `<`))
#   }
#
#   # Return data frame
#   data.frame(truth = params$truth, response = response, confidence = confidence)
# }
#
# #' @rdname conf_sdt
# #' @export
# dconf_sdt <- function(truth, response, confidence, dprime, c, thetazero, thetaone, log = FALSE) {
#   # --- Prepare and Validate Parameters ---
#   params <- .prepare_conf_sdt(truth = truth, response = response, confidence = confidence, dprime = dprime, c = c)
#   if (length(c) > 1) stop("'c' must be a single value.")
#
#   decision_criterion <- -c
#
#   if (any(thetaone <= decision_criterion)) stop("All thetaone must be > -c.")
#   if (any(thetazero >= decision_criterion)) stop("All thetazero must be < -c.")
#   if (is.unsorted(thetaone)) stop("thetaone must be sorted in increasing order.")
#   if (is.unsorted(rev(thetazero))) stop("thetazero must be sorted in decreasing order.")
#
#   # --- Calculation ---
#   lower_bound <- numeric(params$ndraws)
#   upper_bound <- numeric(params$ndraws)
#
#   # --- response == 1 ---
#   idx_r1 <- which(params$response == 1)
#   if (length(idx_r1) > 0) {
#     conf_r1 <- params$confidence[idx_r1]
#     if (any(conf_r1 > length(thetaone))) {
#       stop("A confidence value for response=1 is higher than allowed by the number of thetaone criteria.")
#     }
#     boundaries <- c(decision_criterion, thetaone, Inf)
#     lower_bound[idx_r1] <- boundaries[conf_r1 + 1]
#     upper_bound[idx_r1] <- boundaries[conf_r1 + 2]
#   }
#
#   # --- response == 0 ---
#   idx_r0 <- which(params$response == 0)
#   if (length(idx_r0) > 0) {
#     conf_r0 <- params$confidence[idx_r0]
#     if (any(conf_r0 > length(thetazero))) {
#       stop("A confidence value for response=0 is higher than allowed by the number of thetazero criteria.")
#     }
#     # Upper bound candidates: -c, tzero1, tzero2, ...
#     upper_bound_cand <- c(decision_criterion, thetazero)
#     upper_bound[idx_r0] <- upper_bound_cand[conf_r0 + 1]
#     # Lower bound candidates: tzero1, tzero2, ..., -Inf
#     lower_bound_cand <- c(thetazero, -Inf)
#     lower_bound[idx_r0] <- lower_bound_cand[conf_r0 + 1]
#   }
#
#   # Calculate probability
#   mean_x <- ((params$truth * 2) - 1) * params$dprime / 2
#   prob <- stats::pnorm(upper_bound, mean = mean_x, sd = 1) - stats::pnorm(lower_bound, mean = mean_x, sd = 1)
#
#   # Return
#   if (log) {
#     return(log(prob))
#   } else {
#     return(prob)
#   }
# }
#
#
# #' @keywords internal
# .prepare_conf_sdt <- function(n = NULL, response = NULL, truth = NULL, confidence = NULL, dprime, c) {
#   # --- Determine Target Length ---
#   if (!is.null(n)) {
#     # For RNG (rconf_sdt)
#     if (length(n) != 1 || n <= 0 || n != floor(n)) {
#       stop("n must be a single positive integer.")
#     }
#     param_lengths <- c(length(dprime), length(truth))
#     m <- max(n, param_lengths)
#   } else {
#     # For Density/Likelihood (dconf_sdt)
#     if (is.null(response) || is.null(confidence) || is.null(truth)) stop("response, truth and confidence must be provided for dconf_sdt.")
#     if (any(!response %in% c(0, 1), na.rm = TRUE)) stop("response must contain only 0 or 1.")
#     if (any(!truth %in% c(0, 1), na.rm = TRUE)) stop("truth must contain only 0 or 1.")
#     if (any(confidence < 0, na.rm = TRUE)) stop("confidence must be >= 0.")
#
#     param_lengths <- c(length(truth), length(response), length(confidence), length(dprime))
#     m <- max(param_lengths)
#     if (m == 0) stop("At least one input vector must have non-zero length.")
#   }
#
#   # --- Recycle Parameters ---
#   params <- list(
#     dprime = rep_len(dprime, m),
#     c = c, # Keep as scalar, checked in main function
#     truth = rep_len(truth, m)
#   )
#
#   # --- Add other variables if provided ---
#   if (!is.null(response)) {
#     params$response <- rep_len(response, m)
#   }
#   if (!is.null(confidence)) {
#     params$confidence <- rep_len(confidence, m)
#   }
#
#   params$ndraws <- m
#   params
# }


# ---------------------------------------------------------------------------
# EXPERIMENTAL / NOT READY: the confidence SDT model is commented out until it
# is finished. Uncomment this file (and the matching entries in NAMESPACE and
# _pkgdown.yml) to re-enable it.
# ---------------------------------------------------------------------------

# # Stanvars ----------------------------------------------------------------
#
# #' @keywords internal
# .conf_sdt_lpdf <- function() {
#   "
# // Log-likelihood for the reparameterized confidence SDT model.
# // Y: observed confidence rating (0, 1, 2, ...).
# real conf_sdt_lpdf(real Y, real mu, real c, real tzeroone, real tzerotwo, real toneone, real tonetwo, int dec, int truth) {
#   real lower_bound;
#   real upper_bound;
#   real decision_criterion = -c;
#   real mean_x = (truth * 2 - 1) * mu / 2.0;
#   real eps = 1e-6; // Small epsilon to avoid numerical issues
#
#   // --- Reparameterization ---
#   // The input parameters are positive offsets. Calculate actual criteria.
#   real crit_toneone = decision_criterion + toneone;
#   real crit_tonetwo = crit_toneone + tonetwo;
#   real crit_tzeroone = decision_criterion - tzeroone;
#   real crit_tzerotwo = crit_tzeroone - tzerotwo;
#
#   // Input validation is now implicit in the parameterization.
#
#   if (dec == 1) {
#     if (Y <= eps) {
#       lower_bound = decision_criterion;
#       upper_bound = crit_toneone;
#     } else if (Y - 1.0 <= eps) {
#       lower_bound = crit_toneone;
#       upper_bound = crit_tonetwo;
#     } else if (Y - 2.0 <= eps) {
#       lower_bound = crit_tonetwo;
#       upper_bound = positive_infinity();
#     } else {
#       return negative_infinity(); // Invalid confidence level
#     }
#   } else if (dec == 0) {
#     if (Y <= eps) {
#       lower_bound = crit_tzeroone;
#       upper_bound = decision_criterion;
#     } else if (Y - 1.0 <= eps) {
#       lower_bound = crit_tzerotwo;
#       upper_bound = crit_tzeroone;
#     } else if (Y - 2.0 <= eps) {
#       lower_bound = negative_infinity();
#       upper_bound = crit_tzerotwo;
#     } else {
#       return negative_infinity(); // Invalid confidence level
#     }
#   } else {
#     return negative_infinity(); // Invalid response
#   }
#
#   // --- Log-likelihood Calculation ---
#   return log_diff_exp(normal_lcdf(upper_bound | mean_x, 1.0),
#                       normal_lcdf(lower_bound | mean_x, 1.0));
# }
# "
# }
#
# #' @rdname conf_sdt
# #' @export
# conf_sdt_stanvars <- function() {
#   brms::stanvar(scode = .conf_sdt_lpdf(), block = "functions")
# }
#
#
# #' @rdname conf_sdt
# #' @param link_mu Link function for the dprime parameter.
# #' @param link_c Link function for the c parameter.
# #' @param link_tzeroone,link_tzerotwo,link_toneone,link_tonetwo Link functions for the criteria offset parameters.
# #' @export
# conf_sdt_custom_family <- function(link_mu = "identity", link_c = "identity",
#                                    link_tzeroone = "log", link_tzerotwo = "log",
#                                    link_toneone = "log", link_tonetwo = "log") {
#   brms::custom_family(
#     name = "conf_sdt",
#     dpars = c("mu", "c", "tzeroone", "tzerotwo", "toneone", "tonetwo"),
#     links = c(link_mu, link_c, link_tzeroone, link_tzerotwo, link_toneone, link_tonetwo),
#     lb = c(NA, NA, 0, 0, 0, 0), # Lower bounds for positive offsets
#     ub = c(NA, NA, NA, NA, NA, NA),
#     type = "real",
#     vars = c("dec[n]", "vint1[n]")
#   )
# }
#
# # brms --------------------------------------------------------------------
#
# #' @rdname conf_sdt
# #' @export
# log_lik_conf_sdt <- function(i, prep) {
#   y <- prep$data$Y[i]
#   if (is.na(y)) {
#     return(NA_real_)
#   }
#
#   # Get parameters (offsets for criteria)
#   dprime <- brms::get_dpar(prep, "mu", i = i)
#   c <- brms::get_dpar(prep, "c", i = i)
#   tzeroone_off <- brms::get_dpar(prep, "tzeroone", i = i)
#   tzerotwo_off <- brms::get_dpar(prep, "tzerotwo", i = i)
#   toneone_off <- brms::get_dpar(prep, "toneone", i = i)
#   tonetwo_off <- brms::get_dpar(prep, "tonetwo", i = i)
#
#   truth <- prep$data$truth[i]
#   response <- prep$data$dec[i]
#
#   # Reconstruct actual criteria from offsets
#   crit_toneone <- -c + toneone_off
#   crit_tonetwo <- crit_toneone + tonetwo_off
#   crit_tzeroone <- -c - tzeroone_off
#   crit_tzerotwo <- crit_tzeroone - tzerotwo_off
#
#   # Call the R density function
#   dconf_sdt(
#     truth = truth,
#     response = response,
#     confidence = y,
#     dprime = dprime,
#     c = c,
#     thetazero = cbind(crit_tzeroone, crit_tzerotwo),
#     thetaone = cbind(crit_toneone, crit_tonetwo),
#     log = TRUE
#   )
# }
#
#
# #' @rdname conf_sdt
# #' @importFrom brms get_dpar
# #' @inheritParams rcogmod_lnr
# #' @export
# posterior_predict_conf_sdt <- function(i, prep, ...) {
#   # Get parameters (offsets for criteria)
#   dprime <- brms::get_dpar(prep, "mu", i = i)
#   c <- brms::get_dpar(prep, "c", i = i)
#   tzeroone_off <- brms::get_dpar(prep, "tzeroone", i = i)
#   tzerotwo_off <- brms::get_dpar(prep, "tzerotwo", i = i)
#   toneone_off <- brms::get_dpar(prep, "toneone", i = i)
#   tonetwo_off <- brms::get_dpar(prep, "tonetwo", i = i)
#
#   truth <- prep$data$truth[i]
#
#   # Reconstruct actual criteria from offsets
#   crit_toneone <- -c + toneone_off
#   crit_tonetwo <- crit_toneone + tonetwo_off
#   crit_tzeroone <- -c - tzeroone_off
#   crit_tzerotwo <- crit_tzeroone - tzerotwo_off
#
#   # Loop over draws because rconf_sdt is not vectorized over the structure of criteria
#   n_draws <- length(dprime)
#   predictions <- numeric(n_draws)
#   for (j in 1:n_draws) {
#     sim <- rconf_sdt(
#       n = 1,
#       dprime = dprime[j],
#       c = c[j],
#       thetazero = c(crit_tzeroone[j], crit_tzerotwo[j]),
#       thetaone = c(crit_toneone[j], crit_tonetwo[j]),
#       truth = truth
#     )
#     if (nrow(sim) > 0) {
#       predictions[j] <- sim$confidence
#     } else {
#       predictions[j] <- NA
#     }
#   }
#   predictions
# }
#
#
# #' @rdname conf_sdt
# #' @export
# posterior_epred_conf_sdt <- function(prep) {
#   # Get parameters (offsets for criteria)
#   dprime <- brms::get_dpar(prep, "mu")
#   c <- brms::get_dpar(prep, "c")
#   tzeroone_off <- brms::get_dpar(prep, "tzeroone")
#   tzerotwo_off <- brms::get_dpar(prep, "tzerotwo")
#   toneone_off <- brms::get_dpar(prep, "toneone")
#   tonetwo_off <- brms::get_dpar(prep, "tonetwo")
#
#   truth <- prep$data$truth
#   response <- prep$data$dec
#
#   # Reconstruct actual criteria from offsets
#   crit_toneone <- -c + toneone_off
#   crit_tonetwo <- crit_toneone + tonetwo_off
#   crit_tzeroone <- -c - tzeroone_off
#   crit_tzerotwo <- crit_tzeroone - tzerotwo_off
#
#   n_draws <- nrow(dprime)
#   n_obs <- ncol(dprime)
#   epred <- matrix(NA, nrow = n_draws, ncol = n_obs)
#   conf_levels <- 0:2
#
#   for (i in 1:n_obs) {
#     probs_k_list <- lapply(conf_levels, function(k) {
#       sapply(1:n_draws, function(j) {
#         dconf_sdt(
#           truth = truth[i],
#           response = response[i],
#           confidence = k,
#           dprime = dprime[j, i],
#           c = c[j, i],
#           thetazero = c(crit_tzeroone[j, i], crit_tzerotwo[j, i]),
#           thetaone = c(crit_toneone[j, i], crit_tonetwo[j, i])
#         )
#       })
#     })
#     probs_k <- do.call(cbind, probs_k_list)
#     epred[, i] <- probs_k %*% conf_levels
#   }
#
#   epred
# }
