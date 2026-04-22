test_that("posterior_combined = likelihood + prior (+ size prior for NB/ZINB)", {
  N <- 3
  z <- matrix(1, N, N)
  y <- matrix(2, N, N)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  params <- rep(0, 5)
  
  upri <- list(component1 = list(
    meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
    sdy   = 1, sdx1   = 1, sdx2   = 1, sdx3   = 1, sdx4   = 1
  ))
  
  # Poisson branch: posterior = likelihood + prior only
  got_pois <- posterior_combined(
    pred_combined = pred_combined,
    params = params, z = z, y = y, x_vars = x_vars,
    component = 1, theta = NULL, N = N,
    use_data_priors = FALSE, user_fixed_priors = upri,
    dist = "Poisson", size = NULL
  )
  ll <- likelihood_combined(
    pred_combined, params, z, y, x_vars, 1, theta = NULL, size = NULL,
    N, "Poisson"
  )
  lp <- prior_combined(params, 1, y, x_vars, z, FALSE, upri)
  expect_equal(got_pois, ll + lp, tolerance = 1e-12)
  
  # NB branch: +size prior
  got_nb <- posterior_combined(
    pred_combined = pred_combined,
    params = params, z = z, y = y, x_vars = x_vars,
    component = 1, theta = NULL, N = N,
    use_data_priors = FALSE, user_fixed_priors = upri,
    dist = "NB", size = 2
  )
  ll_nb <- likelihood_combined(
    pred_combined, params, z, y, x_vars, 1, theta = NULL, size = 2,
    N, "NB"
  )
  sp <- size_prior(2, 1)
  expect_equal(got_nb, ll_nb + lp + sp, tolerance = 1e-12)
})

test_that("NB/ZINB with invalid size raises an error", {
  N <- 2
  z <- matrix(1, N, N); y <- matrix(0, N, N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  expect_error(
    posterior_combined(pred_combined, rep(0, 5), z, y, x_vars,
                       1, theta = NULL, N = N,
                       use_data_priors = TRUE, user_fixed_priors = NULL,
                       dist = "NB", size = -1),
    "size"
  )
})