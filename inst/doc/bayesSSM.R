## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(bayesSSM)
library(ggplot2)

## -----------------------------------------------------------------------------
set.seed(1405)
t_val <- 20
sigma_x <- 1
sigma_y <- 0.5

init_state <- rnorm(1, mean = 0, sd = 1)
x <- numeric(t_val)
y <- numeric(t_val)
x[1] <- sin(init_state) + rnorm(1, mean = 0, sd = sigma_x)
y[1] <- x[1] + rnorm(1, mean = 0, sd = sigma_y)
for (t in 2:t_val) {
  x[t] <- sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
  y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
}
x <- c(init_state, x)

## -----------------------------------------------------------------------------
ggplot() +
  geom_line(aes(x = 0:t_val, y = x), color = "blue", linewidth = 1) + # Latent
  geom_point(aes(x = 1:t_val, y = y), color = "red", size = 2) + # Observed
  labs(
    title = "Simulated Data: Latent State and Observations",
    x = "Time",
    y = "Value",
    caption = "Blue line: Latent state (x), Red points: Observed values (y)"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
init_fn <- function(num_particles) {
  rnorm(num_particles, mean = 0, sd = 1)
}

transition_fn <- function(particles, sigma_x) {
  sin(particles) + rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn <- function(y, particles, sigma_y) {
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

## -----------------------------------------------------------------------------
log_prior_sigma_x <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}

log_prior_sigma_y <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}

log_priors <- list(
  sigma_x = log_prior_sigma_x,
  sigma_y = log_prior_sigma_y
)

## -----------------------------------------------------------------------------
result <- pmmh(
  pf_wrapper = bootstrap_filter,
  y = y,
  m = 1000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  pilot_init_params = list(
    c(sigma_x = 0.4, sigma_y = 0.4),
    c(sigma_x = 0.8, sigma_y = 0.8)
  ),
  burn_in = 500,
  num_chains = 2,
  seed = 1405,
  param_transform = list(
    sigma_x = "log",
    sigma_y = "log"
  ),
  tune_control = default_tune_control(pilot_m = 100, pilot_burn_in = 10),
  verbose = TRUE
)

## -----------------------------------------------------------------------------
print(result)

## -----------------------------------------------------------------------------
chains <- result$theta_chain

## -----------------------------------------------------------------------------
ggplot(chains, aes(x = sigma_x, fill = factor(chain))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density plot of sigma_x chains",
    x = "Value",
    y = "Density",
    fill = "Chain"
  ) +
  theme_minimal()

