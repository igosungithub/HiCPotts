test_that("likelihood_combined calculates likelihood correctly for ZIP", {
  # Mock data
  z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
  y <- matrix(c(0, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- list()
  N <- 3
  theta <- 0.2
  size <- 10
  
  # Mock pred_combined function
  pred_combined <- function(params, z, x_vars, component, N) {
    # Return a simple function of params for testing
    rep(sum(params) + component, sum(z == component))
  }
  
  # Run the function
  likelihood <- likelihood_combined(pred_combined, params, z, y, x_vars, 1, theta, size, N, "ZIP")
  
  # Check that likelihood is numeric
  expect_type(likelihood, "double")
})

test_that("likelihood_combined handles Poisson distribution", {
  # Mock data
  z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
  y <- matrix(c(1, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- list()
  N <- 3
  theta <- 0.2
  size <- 10
  
  # Mock pred_combined function
  pred_combined <- function(params, z, x_vars, component, N) {
    rep(sum(params) + component, sum(z == component))
  }
  
  # Run the function
  likelihood <- likelihood_combined(pred_combined, params, z, y, x_vars, 1, theta, size, N, "Poisson")
  
  # Check that likelihood is numeric
  expect_type(likelihood, "double")
})

test_that("likelihood_combined handles invalid distribution", {
  # Mock data
  z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
  y <- matrix(c(1, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- list()
  N <- 3
  theta <- 0.2
  size <- 10
  
  # Mock pred_combined function
  pred_combined <- function(params, z, x_vars, component, N) {
    rep(sum(params) + component, sum(z == component))
  }
  
  # Expect an error for invalid distribution
  expect_error(likelihood_combined(pred_combined, params, z, y, x_vars, 1, theta, size, N, "InvalidDist"),
               "Invalid distribution specified")
})

test_that("likelihood_combined calculates likelihood for component 2 and NB", {
  # Mock data
  z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
  y <- matrix(c(1, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- list()
  N <- 3
  theta <- 0.2
  size <- 10
  
  # Mock pred_combined function
  pred_combined <- function(params, z, x_vars, component, N) {
    rep(sum(params) + component, sum(z == component))
  }
  
  # Run the function for component 2
  likelihood <- likelihood_combined(pred_combined, params, z, y, x_vars, 2, theta, size, N, "NB")
  
  # Check that likelihood is numeric
  expect_type(likelihood, "double")
})

test_that("likelihood_combined calculates likelihood for ZINB", {
  # Mock data
  z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
  y <- matrix(c(0, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- list()
  N <- 3
  theta <- 0.2
  size <- 10
  
  # Mock pred_combined function
  pred_combined <- function(params, z, x_vars, component, N) {
    rep(sum(params) + component, sum(z == component))
  }
  
  # Run the function for ZINB
  likelihood <- likelihood_combined(pred_combined, params, z, y, x_vars, 1, theta, size, N, "ZINB")
  
  # Check that likelihood is numeric
  expect_type(likelihood, "double")
})
