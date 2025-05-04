test_that("posterior_combined computes general behavior for all components and distributions", {
  # Mock `pred_combined` to return valid predictions
  pred_combined <- function(params, z, x_vars, component, N) {
    return(rep(0, sum(z == component)))  # Return a vector of zeros matching the component count
  }
  
  # Mock dependencies
  likelihood_combined <- function(...) { return(10) }  # Mock likelihood for testing
  prior_combined <- function(...) { return(5) }  # Mock prior for testing
  size_prior <- function(size, component) { return(-2 * component) }  # Mock size prior for testing
  
  # Mock data
  params <- c(1, 2, 3, 4, 5)
  z <- matrix(sample(1:3, 5 * 5, replace = TRUE), nrow = 5)
  y <- matrix(rpois(5 * 5, lambda = 5), nrow = 5)
  x_vars <- list(
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5))
  )
  N <- 5
  use_data_priors <- TRUE
  user_fixed_priors <- list()
  
  distributions <- list(
    "Poisson" = list(expected = 13, theta = NULL, size = NULL),
    "ZIP" = list(expected = 8, theta = 0.5, size = NULL),
    "NB" = list(expected = 11, theta = NULL, size = 10),
    "ZINB" = list(expected = 7, theta = 0.5, size = 20)
  )
  
  for (dist in names(distributions)) {
    dist_params <- distributions[[dist]]
    theta <- dist_params$theta
    size <- dist_params$size
    
    for (component in 1:3) {
      expected <- dist_params$expected - 2 * component
      
      # Run the function
      result <- posterior_combined(
        pred_combined, params, z, y, x_vars, component,
        theta, N, use_data_priors, user_fixed_priors, dist, size
      )
      
      # Validate the result type
      expect_type(result, "double")  # Check that the result is numeric
      expect_length(result, 1)  # Ensure it returns a single value
      
   }
  }
})


test_that("posterior_combined throws error for invalid distribution", {
  pred_combined <- function(params, z, x_vars, component, N) { 
    return(rep(0, sum(z == component))) 
  }
  params <- c(1, 2, 3, 4, 5)
  z <- matrix(sample(1:3, 25, replace = TRUE), nrow = 5)
  y <- matrix(rpois(25, lambda = 5), nrow = 5)
  x_vars <- list(
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5)),
    list(matrix(runif(25), nrow = 5))
  )
  N <- 5
  use_data_priors <- TRUE
  user_fixed_priors <- list()
  theta <- NULL
  size <- NULL
  
  expect_error(
    posterior_combined(
      pred_combined, params, z, y, x_vars, component = 1, 
      theta, N, use_data_priors, user_fixed_priors, "InvalidDist", size
    ),
    "Invalid distribution specified"
  )
})

