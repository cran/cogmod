context("CHOCO")

test_that("rcogmod_choco empirical side-probabilities match theory", {
  skip_on_cran()
  set.seed(42)
  n <- 20000
  tol <- 0.05

  # Small grid of parameters, thinned to a covering subset - see helper-grid.R.
  grid <- covering_grid(
    p = c(0.2, 0.5, 0.8),
    pmid = c(0, 0.1),
    pex = c(0, 0.2),
    bex = c(0.3, 0.7),
    mid = c(0.3, 0.5, 0.7)
  )
  # fix symmetric confidences/precisions for simplicity
  confright <- confleft <- 0.5
  precright <- precleft <- 3

  for (.i in seq_len(nrow(grid))) {
    p <- grid$p[.i]
    pmid <- grid$pmid[.i]
    pex <- grid$pex[.i]
    bex <- grid$bex[.i]
    mid <- grid$mid[.i]

    # theoretical weights
    prob_not_mid <- 1 - pmid
    theo_left <- prob_not_mid * (1 - p)
    theo_mid <- pmid
    theo_right <- prob_not_mid * p

    x <- rcogmod_choco(n,
      p = p,
      confright = confright, precright = precright,
      confleft = confleft, precleft = precleft,
      pex = pex, bex = bex,
      pmid = pmid,
      mid = mid
    )

    emp_left <- mean(x < mid)
    emp_mid <- mean(abs(x - mid) < 1e-9)
    emp_right <- mean(x > mid)

    label_base <- paste0(
      "p=", p, ", pmid=", pmid,
      ", pex=", pex, ", bex=", bex,
      ", mid=", mid
    )

    expect_equal(emp_left, theo_left,
      tolerance = tol,
      label = paste(label_base, "— left mass")
    )
    expect_equal(emp_mid, theo_mid,
      tolerance = tol,
      label = paste(label_base, "— mid mass")
    )
    expect_equal(emp_right, theo_right,
      tolerance = tol,
      label = paste(label_base, "— right mass")
    )
  }
})


test_that("dcogmod_choco places correct point‐masses and integrates to 1", {
  set.seed(7)
  tol_mass <- 0.05
  tol_int <- 0.02

  # same grid
  grid <- covering_grid(
    p = c(0.2, 0.5, 0.8),
    pmid = c(0, 0.1),
    pex = c(0, 0.2),
    bex = c(0.3, 0.7),
    mid = c(0.3, 0.5, 0.7)
  )
  confright <- confleft <- 0.5
  precright <- precleft <- 3

  for (.i in seq_len(nrow(grid))) {
    p <- grid$p[.i]
    pmid <- grid$pmid[.i]
    pex <- grid$pex[.i]
    bex <- grid$bex[.i]
    mid <- grid$mid[.i]

    # Theoretical side‐weights
    prob_not_mid <- 1 - pmid
    prob_left <- prob_not_mid * (1 - p)
    prob_mid <- pmid
    prob_right <- prob_not_mid * p

    # Underlying cogmod_betagate masses
    mass0_left <- dcogmod_betagate(0,
      mu = 1 - confleft, phi = precleft,
      pex = pex * (1 - bex), bex = 0
    )
    mass1_right <- dcogmod_betagate(1,
      mu = confright, phi = precright,
      pex = pex * bex, bex = 1
    )

    # Theoretical point‐mass for full CHOCO
    theo_mass0 <- prob_left * mass0_left
    theo_massT <- prob_mid
    theo_mass1 <- prob_right * mass1_right

    # Evaluate dcogmod_choco at the three points
    got0 <- dcogmod_choco(0,
      p = p,
      confright = confright, precright = precright,
      confleft = confleft, precleft = precleft,
      pex = pex, bex = bex,
      pmid = pmid,
      mid = mid
    )
    gotT <- dcogmod_choco(mid,
      p = p,
      confright = confright, precright = precright,
      confleft = confleft, precleft = precleft,
      pex = pex, bex = bex,
      pmid = pmid,
      mid = mid
    )
    got1 <- dcogmod_choco(1,
      p = p,
      confright = confright, precright = precright,
      confleft = confleft, precleft = precleft,
      pex = pex, bex = bex,
      pmid = pmid,
      mid = mid
    )

    base_lbl <- paste0(
      "p=", p, ", pmid=", pmid,
      ", pex=", pex, ", bex=", bex,
      ", mid=", mid
    )

    expect_equal(got0, theo_mass0,
      tolerance = tol_mass,
      label = paste(base_lbl, "— point mass at 0")
    )
    expect_equal(gotT, theo_massT,
      tolerance = tol_mass,
      label = paste(base_lbl, "— point mass at mid")
    )
    expect_equal(got1, theo_mass1,
      tolerance = tol_mass,
      label = paste(base_lbl, "— point mass at 1")
    )

    # Quick numeric integration over (0, mid) and (mid,1)
    f_int <- function(f_lower, f_upper) {
      stats::integrate(
        function(xx) {
          dcogmod_choco(xx,
            p = p,
            confright = confright, precright = precright,
            confleft = confleft, precleft = precleft,
            pex = pex, bex = bex,
            pmid = pmid,
            mid = mid
          )
        },
        lower = f_lower, upper = f_upper,
        rel.tol = tol_int, subdivisions = 200
      )$value
    }

    cont_left <- if (mid > 0) f_int(0 + 1e-8, mid - 1e-8) else 0
    cont_right <- if (mid < 1) f_int(mid + 1e-8, 1 - 1e-8) else 0

    total_mass <- theo_mass0 + theo_massT + theo_mass1 + cont_left + cont_right
    expect_equal(total_mass, 1,
      tolerance = tol_int,
      label = paste(base_lbl, "— total integrates to 1")
    )
  }
})



context("CHOCO - brms")

test_that("CHOCO model can recover parameters with brms using variational inference", {
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan") # Needed by brms for custom families sometimes

  # --- 1. Set up Simulation ---
  true_mu <- 0.7
  true_confright <- 0.8
  true_precright <- 5
  true_confleft <- 0.6
  true_precleft <- 5
  true_pex <- 0.1
  true_bex <- 0.6
  true_pmid <- 0
  n_obs <- 5000

  # Generate synthetic data with known parameters
  set.seed(1234)
  df <- data.frame(
    y = rcogmod_choco(
      n = n_obs,
      p = true_mu,
      confright = true_confright,
      precright = true_precright,
      confleft = true_confleft,
      precleft = true_precleft,
      pex = true_pex,
      bex = true_bex,
      pmid = true_pmid,
      mid = 0.5
    )
  )

  # --- 2. Define and Fit brms Model ---
  f <- brms::bf(
    y ~ 1,
    confright ~ 1,
    precright ~ 1,
    confleft ~ 1,
    precleft ~ 1,
    pex ~ 1,
    bex ~ 1,
    pmid = 0,
    family = cogmod_choco()
  )

  m <- brms::brm(
    formula = f,
    data = df,
    init = 0,
    stanvars = cogmod_choco_stanvars(),
    backend = "cmdstanr",
    seed = 1234,
    refresh = 0,
    algorithm = "pathfinder"
  )

  # --- 3. Check Parameter Recovery (point estimate vs true value) ---
  post_summary <- brms::posterior_summary(m)

  # Mu
  mu_est <- brms::inv_logit_scaled(post_summary["b_Intercept", "Estimate"])
  expect_equal(mu_est, true_mu, tolerance = 0.15, label = "Recovered mu")

  # Confright
  confright_est <- brms::inv_logit_scaled(post_summary["b_confright_Intercept", "Estimate"])
  expect_equal(confright_est, true_confright, tolerance = 0.15, label = "Recovered confright")

  # Precright (softplus link)
  precright_est <- log(1 + exp(post_summary["b_precright_Intercept", "Estimate"]))
  expect_equal(precright_est, true_precright, tolerance = 0.5, label = "Recovered precright")

  # Confleft
  confleft_est <- brms::inv_logit_scaled(post_summary["b_confleft_Intercept", "Estimate"])
  expect_equal(confleft_est, true_confleft, tolerance = 0.15, label = "Recovered confleft")

  # Precleft (softplus link)
  precleft_est <- log(1 + exp(post_summary["b_precleft_Intercept", "Estimate"]))
  expect_equal(precleft_est, true_precleft, tolerance = 0.5, label = "Recovered precleft")

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


test_that("Stan cogmod_choco_lpdf matches R dcogmod_choco function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function if possible
  cogmod_choco_lpdf <- stan_fun("cogmod_choco")

  # --- Define parameter grids for testing ---
  y_values <- c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)
  mu_values <- c(0.2, 0.5, 0.8)
  confright_values <- c(0.3, 0.7)
  precright_values <- c(2, 5)
  confleft_values <- c(0.3, 0.7)
  precleft_values <- c(2, 5)
  pex_values <- c(0, 0.3, 0.7)
  bex_values <- c(0, 0.5, 1)
  pmid_values <- c(0, 0.1)

  for (mu in mu_values) {
    for (confright in confright_values) {
      for (precright in precright_values) {
        for (confleft in confleft_values) {
          for (precleft in precleft_values) {
            for (pex in pex_values) {
              for (bex in bex_values) {
                for (pmid in pmid_values) {
                  for (y in y_values) {
                    label <- sprintf(
                      "y=%.2f, mu=%.1f, confright=%.1f, precright=%.1f, confleft=%.1f, precleft=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                      y, mu, confright, precright, confleft, precleft, pex, bex, pmid
                    )

                    # Calculate log-density using Stan function
                    stan_log_lik <- cogmod_choco_lpdf(
                      y, mu, confright, precright, confleft, precleft, pex, bex, pmid
                    )

                    # Calculate log-density using R function
                    r_log_lik <- dcogmod_choco(
                      x = y,
                      p = mu,
                      confright = confright,
                      precright = precright,
                      confleft = confleft,
                      precleft = precleft,
                      pex = pex,
                      bex = bex,
                      pmid = pmid,
                      mid = 0.5,
                      log = TRUE
                    )

                    expect_equal(stan_log_lik, r_log_lik, tolerance = 1e-6,
                      label = paste("Log-likelihood comparison:", label))
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  # --- Test invalid parameter handling ---
  expect_error(
    dcogmod_choco(0.5, p = -0.1, confright = 0.5, precright = 5, confleft = 0.5, precleft = 5, pex = 0.1, bex = 0.5, pmid = 0.1, mid = 0.5),
    "p must be between 0 and 1"
  )
  expect_error(
    dcogmod_choco(0.5, p = 0.5, confright = -0.1, precright = 5, confleft = 0.5, precleft = 5, pex = 0.1, bex = 0.5, pmid = 0.1, mid = 0.5),
    "confright must be between 0 and 1"
  )
  expect_error(
    dcogmod_choco(0.5, p = 0.5, confright = 0.5, precright = -1, confleft = 0.5, precleft = 5, pex = 0.1, bex = 0.5, pmid = 0.1, mid = 0.5),
    "precright must be positive"
  )
  expect_error(
    dcogmod_choco(0.5, p = 0.5, confright = 0.5, precright = 5, confleft = -0.1, precleft = 5, pex = 0.1, bex = 0.5, pmid = 0.1, mid = 0.5),
    "confleft must be between 0 and 1"
  )
  expect_error(
    dcogmod_choco(0.5, p = 0.5, confright = 0.5, precright = 5, confleft = 0.5, precleft = -1, pex = 0.1, bex = 0.5, pmid = 0.1, mid = 0.5),
    "precleft must be positive"
  )
})
