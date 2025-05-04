test_that("likelihood_gamma returns correct matrix dimensions", {
  # Mock inputs
  x <- c(0.1, 0.2, 0.3, 0.4)
  pair_neighbours_DA_x1 <- c(1, 2, 3, 4)
  N <- 2
  
  # Run the function
  result <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
  
  # Check dimensions of the result
  expect_equal(dim(result), c(N, N))
})

test_that("likelihood_gamma normalizes probabilities correctly", {
  # Mock inputs
  x <- c(1, 2, 3, 4)
  pair_neighbours_DA_x1 <- c(2, 3, 4, 5)
  N <- 2
  
  # Run the function
  result <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
  
  # Check that all elements are non-negative
  expect_true(all(result >= 0))
  
  # Check that the sum of all elements equals 1
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("likelihood_gamma handles numerical stability with large inputs", {
  # Mock inputs with large values
  x <- c(100, 200, 300, 400)
  pair_neighbours_DA_x1 <- c(1, 2, 3, 4)
  N <- 2
  
  # Run the function
  result <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
  
  # Check that the result contains valid probabilities
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("likelihood_gamma handles small inputs", {
  # Mock inputs with small values
  x <- c(1e-5, 1e-4, 1e-3, 1e-2)
  pair_neighbours_DA_x1 <- c(1, 1, 1, 1)
  N <- 2
  
  # Run the function
  result <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
  
  # Check that the result contains valid probabilities
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1, tolerance = 1e-6)
})

test_that("likelihood_gamma throws error for mismatched lengths", {
  # Mock inputs with mismatched lengths
  x <- c(1, 2, 3)
  pair_neighbours_DA_x1 <- c(1, 2, 3, 4)
  N <- 2
  
  # Expect an error
  expect_error(
    likelihood_gamma(x, pair_neighbours_DA_x1, N),
    "x and pair_neighbours_DA_x1 must have the same length"
  )
})
