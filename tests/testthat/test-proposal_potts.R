test_that("proposalfunction generates a valid beta-distributed value", {
  # Call the function
  result <- proposalfunction()

  # Check that the result is numeric
  expect_type(result, "double")

  # Check that the result is within the valid range for a Beta distribution (0, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("proposalfunction generates different values across multiple calls", {
  # Generate multiple values
  results <- replicate(25, proposalfunction())

  # Check that there is variability in the results
  expect_gt(length(unique(results)), 1)

  # Ensure all generated values are within the valid range
  expect_true(all(results >= 0 & results <= 1))
})
