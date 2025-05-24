test_that("prior_combined calculates correct priors with data-driven priors for component 1", {
  # Mock data
  params <- c(5, 1, 2, 0, 1)
  component <- 1
  y <- rnorm(25)
  x_vars <- list(
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5))
  )
  z <- sample(1:3, 25, replace = TRUE)
  use_data_priors <- TRUE
  user_fixed_priors <- NULL

  # Call the function
  result <- prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors)

  # Validate result is numeric
  expect_type(result, "double")
})

test_that("prior_combined calculates correct priors with fixed priors for component 2", {
  # Mock data
  params <- c(3, 2, 4, 5, 1)
  component <- 2
  y <- rnorm(25)
  x_vars <- list(
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5))
  )
  z <- sample(1:3, 25, replace = TRUE)
  use_data_priors <- FALSE
  user_fixed_priors <- list(
    component1 = list(meany = 5, meanx1 = 1, meanx2 = 2, meanx3 = 3, meanx4 = 4, sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1),
    component2 = list(meany = 700, meanx1 = 2, meanx2 = 3, meanx3 = 4, meanx4 = 5, sdy = 2, sdx1 = 2, sdx2 = 2, sdx3 = 2, sdx4 = 2),
    component3 = list(meany = 10, meanx1 = 3, meanx2 = 4, meanx3 = 5, meanx4 = 6, sdy = 3, sdx1 = 3, sdx2 = 3, sdx3 = 3, sdx4 = 3)
  )

  # Call the function
  result <- prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors)

  # Validate result is numeric
  expect_type(result, "double")
})

test_that("prior_combined throws an error for invalid component", {
  # Mock data
  params <- c(5, 1, 2, 0, 1)
  component <- 4
  y <- rnorm(25)
  x_vars <- list(
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5))
  )
  z <- sample(1:3, 25, replace = TRUE)
  use_data_priors <- FALSE
  user_fixed_priors <- list(
    component1 = list(meany = 5, meanx1 = 1, meanx2 = 2, meanx3 = 3, meanx4 = 4, sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1)
  )

  # Expect error
  expect_error(
    prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors),
    "Invalid component specified"
  )
})

test_that("prior_combined handles small standard deviations correctly", {
  # Mock data
  params <- c(5, 1, 2, 0, 1)
  component <- 1
  y <- rnorm(25)
  x_vars <- list(
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5))
  )
  z <- sample(1:3, 25, replace = TRUE)
  use_data_priors <- TRUE
  user_fixed_priors <- NULL

  # Mock very small variances in x_vars
  x_vars[[1]][[1]] <- matrix(rnorm(25, sd = 1e-8), nrow = 5)

  # Call the function
  result <- prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors)

  # Validate result is numeric
  expect_type(result, "double")
  expect_false(is.nan(result)) # Ensure the result is not NaN
})

test_that("prior_combined returns a valid prior for all components", {
  # Mock data
  params <- c(5, 1, 2, 0, 1)
  y <- rnorm(25)
  x_vars <- list(
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5)),
    list(matrix(rnorm(25), nrow = 5))
  )
  z <- sample(1:3, 25, replace = TRUE)
  use_data_priors <- TRUE
  user_fixed_priors <- NULL

  # Test all components
  for (component in 1:3) {
    result <- prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors)
    expect_type(result, "double")
  }
})
