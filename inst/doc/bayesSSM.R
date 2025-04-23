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
t_val <- 50
phi_true <- 0.8
sigma_x_true <- 1
sigma_y_true <- 0.5

x <- numeric(t_val)
y <- numeric(t_val)
x[1] <- rnorm(1)
y[1] <- x[1] + sigma_y_true * rnorm(1)
for (t in 2:t_val) {
  x[t] <- phi_true * x[t - 1] + sin(x[t - 1]) + sigma_x_true * rnorm(1)
  y[t] <- x[t] + sigma_y_true * rnorm(1)
}

## -----------------------------------------------------------------------------
ggplot() +
  geom_line(aes(x = 1:t_val, y = x), color = "blue", linewidth = 1) + # Latent
  geom_point(aes(x = 1:t_val, y = y), color = "red", size = 2) + # Observed
  labs(
    title = "Simulated Data: Latent State and Observations",
    x = "Time",
    y = "Value",
    caption = "Blue line: Latent state (x), Red points: Observed values (y)"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
init_fn <- function(particles) {
  rnorm(particles, mean = 0, sd = 1)
}

transition_fn <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn <- function(y, particles, sigma_y) {
  dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}

## -----------------------------------------------------------------------------
log_prior_phi <- function(phi) {
  dunif(phi, min = 0, max = 1, log = TRUE)
}

log_prior_sigma_x <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}

log_prior_sigma_y <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}

log_priors <- list(
  phi = log_prior_phi,
  sigma_x = log_prior_sigma_x,
  sigma_y = log_prior_sigma_y
)

## -----------------------------------------------------------------------------
result <- pmmh(
  y = y,
  m = 1000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  pilot_init_params = list(
    c(phi = 0.4, sigma_x = 0.4, sigma_y = 0.4),
    c(phi = 0.8, sigma_x = 0.8, sigma_y = 0.8)
  ),
  burn_in = 500,
  num_chains = 2,
  seed = 1405,
  param_transform = list(
    phi = "identity",
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
ggplot(chains, aes(x = phi, fill = factor(chain))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density plot of phi chains",
    x = "Value",
    y = "Density",
    fill = "Chain"
  ) +
  theme_minimal()

