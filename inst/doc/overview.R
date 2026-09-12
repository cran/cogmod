## -----------------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## -----------------------------------------------------------------------------
library(cogmod)

set.seed(3)
x <- rcogmod_betagate(5000, mu = 0.6, phi = 4, pex = 0.15, bex = 0.4)
hist(x, breaks = 50, col = "#2196F3", border = NA,
     main = "Beta-Gate ratings", xlab = "Rating")

## -----------------------------------------------------------------------------
# library(brms)
# 
# f <- bf(rating ~ condition + (1 | participant), phi ~ 1, pex ~ 1, bex ~ 1,
#         family = cogmod_betagate())
# 
# m <- brm(
#   f,
#   data = df,
#   stanvars = cogmod_stanvars(f),
#   backend = "cmdstanr"
# )

## -----------------------------------------------------------------------------
citation("cogmod")

