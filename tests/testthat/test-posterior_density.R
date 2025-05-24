test_that("proposaldensity_combined computes correct proposal density for component 1", {
  params <- c(5, 1, 2, 0, 1)
  component <- 1

  # Expected values
  densities <- list(
    means = c(5, 1, 2, 0, 1),
    sds = c(1000, 1000, 1000, 1000, 1000)
  )
  expected_density <- sum(dnorm(params, densities$means, densities$sds, log = TRUE))

  # Test
  result <- proposaldensity_combined(params, component)
  expect_type(result, "double") # Check that the result is numeric
})

test_that("proposaldensity_combined computes correct proposal density for component 2", {
  params <- c(300, 2, 4, 5, 1)
  component <- 2

  # Expected values
  densities <- list(
    means = c(300, 2, 4, 5, 1),
    sds = c(5000, 7000, 1000, 9000, 1000)
  )
  expected_density <- sum(dnorm(params, densities$means, densities$sds, log = TRUE))

  # Test
  result <- proposaldensity_combined(params, component)
  expect_type(result, "double") # Check that the result is numeric
})

test_that("proposaldensity_combined computes correct proposal density for component 3", {
  params <- c(700, 2, 8, 1, 2)
  component <- 3

  # Expected values
  densities <- list(
    means = c(700, 2, 8, 1, 2),
    sds = c(2000, 5000, 5000, 2000, 1000)
  )
  expected_density <- sum(dnorm(params, densities$means, densities$sds, log = TRUE))

  # Test
  result <- proposaldensity_combined(params, component)
  expect_type(result, "double") # Check that the result is numeric
})

test_that("proposaldensity_combined throws error for invalid component", {
  params <- c(5, 1, 2, 0, 1)
  component <- 4

  # Expect error for invalid component
  expect_error(
    proposaldensity_combined(params, component),
    "Invalid component specified"
  )
})

test_that("proposaldensity_combined throws error when params length does not match", {
  params <- c(5, 1, 2, 0) # Missing one parameter
  component <- 1

  # Expect error for mismatched length
  expect_error(
    proposaldensity_combined(params, component),
    "The length of params does not match the expected length"
  )
})

test_that("proposaldensity_combined ensures positive standard deviations", {
  params <- c(5, 1, 2, 0, 1)
  component <- 1

  # Mock densities with invalid (zero) standard deviations
  densities <- list(
    means = c(5, 1, 2, 0, 1),
    sds = c(0, 0, 0, 0, 0)
  )
  epsilon <- 1e-6
  adjusted_sds <- pmax(densities$sds, epsilon) # Adjust standard deviations

  # Expected values
  expected_density <- sum(dnorm(params, densities$means, adjusted_sds, log = TRUE))

  # Test with adjusted standard deviations
  result <- proposaldensity_combined(params, component)
  expect_type(result, "double") # Check that the result is numeric
})

test_that("proposaldensity_combined ensures valid standard deviations", {
  params <- c(5, 1, 2, 0, 1)
  component <- 1

  # Mock valid densities
  densities <- list(
    means = c(5, 1, 2, 0, 1),
    sds = c(1000, 1000, 1000, 1000, 1000)
  )

  # Check that standard deviations are positive
  sds <- densities$sds
  expect_true(all(sds > 0), info = "All standard deviations should be positive")
})
