test_that("run_metropolis_MCMC_betas executes without errors and returns valid output", {
  # Example inputs
  N <- 5
  gamma_prior <- 0.3
  iterations <- 10
  x_vars <- list(
    distance = list(matrix(runif(N * N, 0, 1), nrow = N)),
    GC = list(matrix(runif(N * N, 0, 1), nrow = N)),
    TES = list(matrix(runif(N * N, 0, 1), nrow = N)),
    ACC = list(matrix(runif(N * N, 0, 1), nrow = N))
  )
  y <- matrix(rpois(N * N, lambda = 5), nrow = N)
  use_data_priors <- TRUE
  user_fixed_priors <- NULL
  dist <- "ZIP"
  epsilon <- NULL
  distance_metric <- "manhattan"
  size_start <- c(1, 1, 1)  # Correct type for Nullable<NumericMatrix>
  theta_start <- 0.5
  
  # Run the function
  result <- run_metropolis_MCMC_betas(
    N = N,
    gamma_prior = gamma_prior,
    iterations = iterations,
    x_vars = x_vars,
    y = y,
    use_data_priors = use_data_priors,
    dist = dist,
    theta_start <- 0.5,
    distance_metric = distance_metric
  )
  
  # Validate output structure
  expect_type(result, "list")
  expect_named(result, c("chains", "gamma", "theta", "size"))
  expect_type(result$chains, "list")
  expect_type(result$gamma, "double")
  expect_type(result$theta, "double")
  expect_type(result$size, "double")
})

test_that("Neighbours_combined calculates neighbors correctly", {
  # Example inputs
  potts_data <- matrix(sample(1:3, 25, replace = TRUE), nrow = 5)
  N <- nrow(potts_data)
  proposed_value <- matrix(sample(1:3, 25, replace = TRUE), nrow = 5)
  
  # Run the function
  result <- Neighbours_combined(potts_data, N, proposed_value)
  
  # Validate output
  expect_true(is.matrix(result), info = "Result should be a matrix")
  expect_equal(dim(result), dim(potts_data), info = "Result dimensions should match input dimensions")
  expect_true(all(result >= 0), info = "All values in the result should be non-negative")
})

test_that("pz_123 calculates probabilities correctly", {
  # Example inputs
  z <- matrix(sample(1:3, 25, replace = TRUE), nrow = 5)
  sum_neighbours <- matrix(runif(25, 0, 5), nrow = 5)
  y <- matrix(rpois(25, lambda = 5), nrow = 5)
  pred_combined <- function(params, z, x_vars, component, N) {
    return(runif(sum(z == component)))
  }
  chains <- list(
    matrix(rnorm(25), nrow = 5, ncol = 5),
    matrix(rnorm(25), nrow = 5, ncol = 5),
    matrix(rnorm(25), nrow = 5, ncol = 5)
  )
  chain_gamma <- rnorm(5)
  x_vars <- list(
    distance = list(matrix(runif(25, 0, 1), nrow = 5)),
    GC = list(matrix(runif(25, 0, 1), nrow = 5)),
    TES = list(matrix(runif(25, 0, 1), nrow = 5)),
    ACC = list(matrix(runif(25, 0, 1), nrow = 5))
  )
  theta <- 0.5
  size_chain <- matrix(runif(25, 0.5, 1.5), nrow = 5, ncol = 5)
  N <- nrow(z)
  iter <- 1
  dist <- "ZIP"
  
  # Run the function
  result <- pz_123(
    z = z,
    sum_neighbours = sum_neighbours,
    y = y,
    pred_combined = pred_combined,
    chains = chains,
    chain_gamma = chain_gamma,
    x_vars = x_vars,
    theta = theta,
    size_chain = size_chain,
    N = N,
    iter = iter,
    dist = dist
  )
  
  # Validate output
  expect_true(is.matrix(result), info = "Result should be a matrix")
  expect_equal(dim(result), dim(z), info = "Result dimensions should match input dimensions")
  expect_true(all(is.finite(result)), info = "Result should not contain NaN or Inf values")
})

