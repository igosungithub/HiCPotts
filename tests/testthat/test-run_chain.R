
test_that("run_metropolis_MCMC_betas executes MCMC chains correctly for valid inputs", {
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
  y <- matrix(rpois(N * N, lambda = 5), nrow = N)  # Wrap y in a list
  theta_start <- 0.5
  use_data_priors <- TRUE
  dist <- "ZIP"
  distance_metric <- "manhattan"
  
  # Run the function
  result <- run_metropolis_MCMC_betas(
    N = N,
    gamma_prior = gamma_prior,
    iterations = iterations,
    x_vars = x_vars,
    y = y,
    theta_start = theta_start,
    use_data_priors = use_data_priors,
    dist = dist,
    distance_metric = distance_metric
  )
  
  # Validate output structure
  expect_type(result, "list")
  expect_named(result, c("chains", "gamma", "theta","size"))
  expect_type(result$chains, "list")
  expect_type(result$gamma, "double")
  expect_type(result$theta, "double")
})


test_that("run_metropolis_MCMC_betas works for multiple distributions", {
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
  y <- matrix(rpois(N * N, lambda = 5), nrow = N)  # Wrap y in a list
  theta_start <- 0.5
  size_start <- c(1, 1, 1)
  use_data_priors <- TRUE
  distance_metric <- "manhattan"
  
  # Test the function with different distributions
  for (dist in c("Poisson", "ZIP", "NB", "ZINB")) {
    result <- run_metropolis_MCMC_betas(
      N = N,
      gamma_prior = gamma_prior,
      iterations = iterations,
      x_vars = x_vars,
      y = y,
      theta_start = theta_start,
      size_start = size_start,
      use_data_priors = use_data_priors,
      dist = dist,
      distance_metric = distance_metric
    )
    
    # Validate output
    expect_type(result, "list")
    expect_length(result, 4)  # Output should match the number of chains
    expect_named(result, c("chains", "gamma", "theta", "size"))
   }
})


test_that("run_metropolis_MCMC_betas works for user-fixed prior", {
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
  y <- matrix(rpois(N * N, lambda = 5), nrow = N)  # Wrap y in a list
  use_data_priors <- FALSE  # Disable data priors to test user-fixed priors
  dist <- "ZINB"
  epsilon <- NULL
  distance_metric <- "manhattan"
  size_start <- c(1, 1, 1)
  theta_start <- 0.5
  
  # User-fixed priors
  user_fixed_priors <- list(
    component1 = list(
      meany = 3, meanx1 = 5, meanx2 = 0.2, meanx3 = 0, meanx4 = 4,
      sdy = 1, sdx1 = 0.01, sdx2 = 0.01, sdx3 = 0.002, sdx4 = 0.005
    ),
    component2 = list(
      meany = 300, meanx1 = 8, meanx2 = 0.3, meanx3 = 3, meanx4 = 8,
      sdy = 0.1, sdx1 = 0.1, sdx2 = 0.01, sdx3 = 0.02, sdx4 = 0.02
    ),
    component3 = list(
      meany = 200, meanx1 = 10, meanx2 = 0.3, meanx3 = 1.6, meanx4 = 4.5,
      sdy = 0.1, sdx1 = 0.6, sdx2 = 0.1, sdx3 = 0.2, sdx4 = 0.2
    )
  )
  
  # Run the function
  result <- run_metropolis_MCMC_betas(
    N = N,
    gamma_prior = gamma_prior,
    iterations = iterations,
    x_vars = x_vars,
    y = y,
    theta_start = theta_start,
    size_start = size_start,
    use_data_priors = use_data_priors,
    user_fixed_priors = user_fixed_priors,
    dist = dist,
    distance_metric = distance_metric
  )
  
  # Validate output
  expect_type(result, "list")
  expect_length(result, 4)
  expect_named(result, c("chains", "gamma", "theta", "size"))
})


