# ORDERED BETA REGRESSION (https://github.com/saudiwin/ordbetareg_pack/blob/master/R/distribution.R)
# For reference


#' Generate Ordered Beta Variates
#'
#'
#' This function will generate ordered beta random variates given
#' values for the mean (`mu`), dispersion (`phi`) and cutpoints
#' governing the ratio of degenerate (discrete) to continuous
#' responses.
#'
#' @param n Number of variates to generate.
#' @param mu Value of the mean of the distribution.
#' Should be in the \(0,1\) interval (cannot be strictly equal to 0 or 1). If
#' length is greater than 1, should be of length `n`.
#' @param phi Value of the dispersion parameter. Should be strictly greater than 0. If
#' length is greater than 1, should be of length `n`.
#' @param cutzero,cutone Two numeric values for the cutpoints on the logit scale. Second value should
#' be strictly greater than the first value.
#' @return A vector of length `n` of variates from the ordered beta distribution.
#' 
#' @keywords internal
#' @noRd
.rordbeta <- function(n = 100,
                     mu = 0.5,
                     phi = 1,
                     cutzero = -1,
                     cutone = 1) {
  if (!all(mu > 0 & mu < 1))  stop("Please pass a numeric value for mu that is between 0 and 1.")
  if (!all(phi > 0)) stop("Please pass a numeric value for phi that is greater than 0.")
  if (!(length(mu) %in% c(1, n)))  stop("Please pass a vector for mu that is either length 1 or length N.")
  if (!(length(phi) %in% c(1, n))) stop("Please pass a vector for phi that is either length 1 or length N.")

  mu_ql <- stats::qlogis(mu)
  if (length(mu_ql) == 1) {
    mu_ql <- rep(mu_ql, n)
  }

  # probabilities for three possible categories (0, proportion, 1)
  low <- 1 - stats::plogis(mu_ql - cutzero)
  middle <- stats::plogis(mu_ql - cutzero) - stats::plogis(mu_ql - cutone)
  high <- stats::plogis(mu_ql - cutone)

  # we'll assume the same eta was used to generate outcomes

  out_beta <- stats::rbeta(n = n, mu * phi, (1 - mu) * phi)

  # now determine which one we get for each observation
  outcomes <- sapply(1:n, function(i) {
    sample(1:3, size = 1, prob = c(low[i], middle[i], high[i]))
  })

  # if we get underflow/overflow, assign the max/min based on double floating point precision
  out_beta[out_beta == 1] <- out_beta[out_beta == 1] - .Machine$double.eps
  out_beta[out_beta == 0] <- out_beta[out_beta == 0] + .Machine$double.eps

  # now combine binary (0/1) with proportion (beta)

  final_out <- sapply(1:n, function(i) {
    if (outcomes[i] == 1) {
      return(0)
    } else if (outcomes[i] == 2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })

  final_out
}






#' @keywords internal
#' @noRd
.dordbeta <- function(x = .9,
                     mu = 0.5,
                     phi = 1,
                     cutzero = -1,
                     cutone = 1,
                     log = FALSE) {
  if (!all(mu > 0 & mu < 1)) stop("Please pass a numeric value for mu that is between 0 and 1.")
  if (!all(phi > 0)) stop("Please pass a numeric value for phi that is greater than 0.")
  if (!(length(mu) %in% c(1, length(x)))) stop("Please pass a vector for mu that is either length 1 or length N.")
  if (!(length(phi) %in% c(1, length(x)))) stop("Please pass a vector for phi that is either length 1 or length N.")

  mu_ql <- stats::qlogis(mu)

  if (length(mu_ql) == 1) {
    mu_ql <- rep(mu_ql, length(x))
  }

  # probabilities for three possible categories (0, proportion, 1)
  low <- 1 - stats::plogis(mu_ql - cutzero)
  middle <- stats::plogis(mu_ql - cutzero) - stats::plogis(mu_ql - cutone)
  high <- stats::plogis(mu_ql - cutone)

  # we'll assume the same eta was used to generate outcomes
  out_beta <- stats::dbeta(x = x, mu * phi, (1 - mu) * phi)

  # now determine which one we get for each observation
  outcomes <- sapply(1:length(x), function(i) {
    if (x[i] == 0) {
      out <- low[i]
    } else if (x[i] == 1) {
      out <- high[i]
    } else {
      out <- middle[i] * out_beta[i]
    }

    return(out)
  })

  if (log) {
    return(log(outcomes))
  } else {
    return(outcomes)
  }
}





# --- STAN LPDF Function ---

#' @keywords internal
#' @noRd
.ordbeta_lpdf <- function() {
  "
real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {
  // cutone must be > cutzero
  if(y==0) {
      return log1m_inv_logit(mu - cutzero);
    } else if(y==1) {
      return log_inv_logit(mu  - cutone);
    } else {
      return log_diff_exp(log_inv_logit(mu - cutzero), log_inv_logit(mu - cutone)) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }
"
}




# Reparametrized Version --------------------------------------------------
# cutpoints are now derived from pex/bex on the identity scale
# phi = phi * 2

# .rordbeta2 <- function(n = 100,
#                      mu = 0.5,
#                      phi = 1,
#                      pex = 0.1,
#                      bex = 0.5) {
#   if (!all(mu > 0 & mu < 1))  stop("Please pass a numeric value for mu that is between 0 and 1.")
#   if (!all(phi > 0)) stop("Please pass a numeric value for phi that is greater than 0.")
#   if (!(length(mu) %in% c(1, n)))  stop("Please pass a vector for mu that is either length 1 or length N.")
#   if (!(length(phi) %in% c(1, n))) stop("Please pass a vector for phi that is either length 1 or length N.")

#   eps <- 1e-10

#   # Cutpoints
#   cutzero <- pex * (1 - bex)
#   cutone <- 1 - pex * bex
#   cutzerolog <- stats::qlogis(cutzero)
#   cutonelog <- stats::qlogis(cutone)

#   phi <- phi * 2
#   mu_ql <- stats::qlogis(mu)
#   if (length(mu_ql) == 1) {
#     mu_ql <- rep(mu_ql, n)
#   }

#   # probabilities for three possible categories (0, proportion, 1)
#   low <- 1 - stats::plogis(mu_ql - cutzerolog)
#   middle <- stats::plogis(mu_ql - cutzerolog) - stats::plogis(mu_ql - cutonelog)
#   high <- stats::plogis(mu_ql - cutonelog)

#   # we'll assume the same eta was used to generate outcomes
#   out_beta <- stats::rbeta(n = n, mu * phi, (1 - mu) * phi)

#   # now determine which one we get for each observation
#   outcomes <- sapply(1:n, function(i) {
#     sample(1:3, size = 1, prob = c(low[i], middle[i], high[i]))
#   })

#   # if we get underflow/overflow, assign the max/min based on double floating point precision
#   out_beta[out_beta == 1] <- out_beta[out_beta == 1] - eps
#   out_beta[out_beta == 0] <- out_beta[out_beta == 0] + eps

#   # now combine binary (0/1) with proportion (beta)
#   final_out <- sapply(1:n, function(i) {
#     if (outcomes[i] == 1) {
#       return(0)
#     } else if (outcomes[i] == 2) {
#       return(out_beta[i])
#     } else {
#       return(1)
#     }
#   })

#   final_out
# }

# .dordbeta2 <- function(x = .9,
#                      mu = 0.5,
#                      phi = 1,
#                      pex = 0.1,
#                      bex = 0.5,
#                      log = FALSE) {
#   if (!all(mu > 0 & mu < 1)) stop("Please pass a numeric value for mu that is between 0 and 1.")
#   if (!all(phi > 0)) stop("Please pass a numeric value for phi that is greater than 0.")
#   if (!(length(mu) %in% c(1, length(x)))) stop("Please pass a vector for mu that is either length 1 or length N.")
#   if (!(length(phi) %in% c(1, length(x)))) stop("Please pass a vector for phi that is either length 1 or length N.")

#   eps <- 1e-10

#   # Cutpoints
#   cutzero <- pex * (1 - bex)
#   cutone <- 1 - pex * bex
#   cutzerolog <- stats::qlogis(cutzero)
#   cutonelog <- stats::qlogis(cutone)

#   phi <- phi * 2
#   mu_ql <- stats::qlogis(mu)
#   if (length(mu_ql) == 1) {
#     mu_ql <- rep(mu_ql, length(x))
#   }

#   # probabilities for three possible categories (0, proportion, 1)
#   low <- 1 - stats::plogis(mu_ql - cutzerolog)
#   middle <- stats::plogis(mu_ql - cutzerolog) - stats::plogis(mu_ql - cutonelog)
#   high <- stats::plogis(mu_ql - cutonelog)

#   # we'll assume the same eta was used to generate outcomes
#   out_beta <- stats::dbeta(x = x, mu * phi, (1 - mu) * phi)

#   # now determine which one we get for each observation
#   outcomes <- sapply(1:length(x), function(i) {
#     if (x[i] == 0) {
#       out <- low[i]
#     } else if (x[i] == 1) {
#       out <- high[i]
#     } else {
#       out <- middle[i] * out_beta[i]
#     }

#     return(out)
#   })

#   if (log) {
#     return(log(outcomes))
#   } else {
#     return(outcomes)
#   }
# }


# # --- STAN LPDF Function ---
# .ordbeta2_lpdf <- function() {
# "
# real ord_beta_reg_lpdf(real y, real mu, real phi, real pex, real bex) {
#   // Cutpoints
#   real cutzero = pex * (1 - bex)
#   real cutone = 1 - pex * bex
#   real cutzerolog = logit(cutzero);
#   real cutonelog = logit(cutone);
#   real philog = log(phi * 2);

#   if(y==0) {
#       return log1m_inv_logit(mu - cutzerolog);
#     } else if(y==1) {
#       return log_inv_logit(mu  - cutonelog);
#     } else {
#       return log_diff_exp(log_inv_logit(mu - cutzerolog), log_inv_logit(mu - cutonelog)) +
#                 beta_lpdf(y|exp(log_inv_logit(mu) + philog),exp(log1m_inv_logit(mu) + philog));
#     }
#   }
# "
# }
