# True response curve: Gaussian, min-max scaled over the range in Romania
true_response <- function(pars, v, rng, n = 100) {
  p <- pars[pars$variable == v, ]
  x <- seq(rng[1, v], rng[2, v], length.out = n)
  y <- dnorm(x, p$mean, p$sd)
  list(x = x, y = (y - min(y)) / (max(y) - min(y)))
}