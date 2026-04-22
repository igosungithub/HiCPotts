test_that("Returns sum of dnorm log-densities for valid component", {
  params <- c(5, 1, 2, 0, 1)
  means <- c(5, 1, 2, 0, 1); sds <- rep(1000, 5)
  expected <- sum(dnorm(params, means, sds, log = TRUE))
  expect_equal(proposaldensity_combined(params, 1), expected, tolerance = 1e-12)
})

test_that("Invalid component raises", {
  expect_error(proposaldensity_combined(c(1,2,3,4,5), 4),
               "Invalid component")
})

test_that("Wrong params length raises", {
  expect_error(proposaldensity_combined(c(1, 2, 3), 1),
               "does not match")
})
