test_that("pred_combined computes predictions correctly for valid inputs across components", {
  # Mock data
  params <- c(1, 2, 3, 4, 5) # Example parameter values
  z <- matrix(c(1, 1, 2, 3, 1, 3, 2, 3, 1), nrow = 3) # Example component matrix
  x_vars <- list(
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x1 values
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x2 values
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x3 values
    list(matrix(runif(9, 1, 10), nrow = 3)) # Random x4 values
  )
  N <- 3

  for (component in 1:3) {
    # Call pred_combined for each component
    result <- pred_combined(params, z, x_vars, component, N)

    # Validate output
    expect_type(result, "double") # Output should be numeric
    expect_length(result, sum(z == component)) # Length matches number of z == component
  }
})

test_that("pred_combined handles edge cases with zeros in x_vars for all components", {
  # Mock data with zeros
  params <- c(1, 2, 3, 4, 5)
  z <- matrix(c(1, 1, 2, 3, 1, 3, 2, 3, 1), nrow = 3)
  x_vars <- list(
    list(matrix(c(0, 2, 0, 3, 4, 0, 1, 0, 0), nrow = 3)), # x1 with zeros
    list(matrix(c(0, 0, 2, 0, 3, 0, 4, 0, 1), nrow = 3)), # x2 with zeros
    list(matrix(c(3, 0, 1, 0, 0, 0, 2, 0, 0), nrow = 3)), # x3 with zeros
    list(matrix(c(0, 0, 0, 1, 0, 2, 0, 3, 0), nrow = 3)) # x4 with zeros
  )
  N <- 3

  for (component in 1:3) {
    # Call pred_combined for each component
    result <- pred_combined(params, z, x_vars, component, N)

    # Validate output
    expect_type(result, "double")
    expect_true(all(is.finite(result))) # Ensure no NaN or Inf in the result
  }
})

test_that("pred_combined throws an error for invalid inputs", {
  params <- c(1, 2, 3, 4, 5)
  z <- matrix(c(1, 1, 2, 3, 1, 3, 2, 3, 1), nrow = 3)
  N <- 3

  # Test missing x_vars
  expect_error(
    pred_combined(params, z, NULL, 1, N),
    "x_vars cannot be NULL" # Expected error message
  )

  # Test mismatched x_vars dimensions
  x_vars <- list(
    list(matrix(runif(6, 1, 10), nrow = 2)), # Incorrect dimensions for x1
    list(matrix(runif(9, 1, 10), nrow = 3)),
    list(matrix(runif(9, 1, 10), nrow = 3)),
    list(matrix(runif(9, 1, 10), nrow = 3))
  )
  expect_error(
    pred_combined(params, z, x_vars, 1, N),
    "Each x_vars matrix must have dimensions N x N" # Expected error message
  )
})


test_that("pred_combined computes predictions correctly for components 2 and 3", {
  # Mock data for components 2 and 3
  params <- c(1, 2, 3, 4, 5)
  z <- matrix(c(2, 2, 2, 3, 3, 3, 3, 3, 3), nrow = 3) # Only components 2 and 3
  x_vars <- list(
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x1 values
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x2 values
    list(matrix(runif(9, 1, 10), nrow = 3)), # Random x3 values
    list(matrix(runif(9, 1, 10), nrow = 3)) # Random x4 values
  )
  N <- 3

  for (component in 2:3) {
    # Call pred_combined for each component
    result <- pred_combined(params, z, x_vars, component, N)

    # Validate output
    expect_type(result, "double") # Output should be numeric
    expect_length(result, sum(z == component)) # Length matches number of z == component
  }
})
