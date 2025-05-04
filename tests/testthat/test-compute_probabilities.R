test_that("compute_HMRFHiC_probabilities works with valid inputs", {
  # Mock data
  data <- data.frame(
    start = c(1, 10, 20),
    end = c(5, 15, 30),
    interactions = c(10, 20, 30),
    GC = c(0.5, 0.8, 0.3),
    TES = c(0.2, 0.5, 0.7),
    ACC = c(0.9, 0.4, 0.6)
  )
  
  # Mock chain_betas with 5x5 matrices
  chain_betas <- list(
    list(
      chains = list(
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5)
      ),
      theta = runif(5, 0.1, 0.9),
      size = matrix(runif(15, 1, 10), nrow = 3)
    )
  )
  
  # Test with ZINB distribution
  result <- compute_HMRFHiC_probabilities(data, chain_betas, iterations = 5, dist = "ZINB")
  
  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_true(all(c("start", "end", "interactions", "prob1", "prob2", "prob3") %in% colnames(result)))
  
  # Check probabilities sum to 1 for each row
  prob_sums <- rowSums(result[, c("prob1", "prob2", "prob3")])
  expect_true(all(abs(prob_sums - 1) < 1e-6))
})

test_that("compute_HMRFHiC_probabilities handles missing required columns", {
  # Mock data with missing columns
  data <- data.frame(
    start = c(1, 10, 20),
    interactions = c(10, 20, 30),
    GC = c(0.5, 0.8, 0.3)
  )
  
  # Mock chain_betas with 5x5 matrices
  chain_betas <- list(
    list(
      chains = list(
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5)
      ),
      theta = runif(5, 0.1, 0.9),
      size = matrix(runif(15, 1, 10), nrow = 3)
    )
  )
  
  # Expect an error due to missing columns
  expect_error(compute_HMRFHiC_probabilities(data, chain_betas, iterations = 5, dist = "ZINB"),
               "The following required columns are missing")
})

test_that("compute_HMRFHiC_probabilities computes probabilities correctly for capped interactions", {
  # Mock data with high interactions
  data <- data.frame(
    start = c(1, 10, 20),
    end = c(5, 15, 30),
    interactions = c(100, 200, 300),  # Exceeds cap
    GC = c(0.5, 0.8, 0.3),
    TES = c(0.2, 0.5, 0.7),
    ACC = c(0.9, 0.4, 0.6)
  )
  
  # Mock chain_betas with 5x5 matrices
  chain_betas <- list(
    list(
      chains = list(
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5),
        matrix(runif(25, 0.1, 1), ncol = 5)
      ),
      theta = runif(5, 0.1, 0.9),
      size = matrix(runif(15, 1, 10), nrow = 3)
    )
  )
  
  # Test function
  result <- compute_HMRFHiC_probabilities(data, chain_betas, iterations = 5, dist = "ZINB")
  
  # Check that interactions are capped at 500
  expect_true(all(result$interactions <= 500))
})
