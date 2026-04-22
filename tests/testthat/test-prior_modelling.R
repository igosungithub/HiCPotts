test_that("prior_combined is DETERMINISTIC (regression test)", {
  set.seed(4291)
  N <- 4
  y      <- rpois(N * N, 5)
  x_vars <- replicate(4, list(matrix(runif(N * N), N, N)), simplify = FALSE)
  z      <- sample(1:3, N * N, replace = TRUE)
  params <- c(1, 0.2, -0.1, 0.3, 0.05)
  
  v1 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  v2 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  v3 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  expect_identical(v1, v2)
  expect_identical(v1, v3)
  expect_true(is.finite(v1))
})

test_that("prior_combined with user_fixed_priors equals sum(dnorm(..., log=TRUE))", {
  params <- c(1, 2, 3, 4, 5)
  upri <- list(
    component1 = list(meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
                      sdy   = 1, sdx1   = 1, sdx2   = 1, sdx3   = 1, sdx4   = 1),
    component2 = list(meany = 1, meanx1 = 2, meanx2 = 3, meanx3 = 4, meanx4 = 5,
                      sdy   = 2, sdx1   = 2, sdx2   = 2, sdx3   = 2, sdx4   = 2),
    component3 = list(meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
                      sdy   = 1, sdx1   = 1, sdx2   = 1, sdx3   = 1, sdx4   = 1)
  )
  expected <- with(upri$component2,
                   dnorm(1, meany, sdy, log = TRUE) +
                     dnorm(2, meanx1, sdx1, log = TRUE) +
                     dnorm(3, meanx2, sdx2, log = TRUE) +
                     dnorm(4, meanx3, sdx3, log = TRUE) +
                     dnorm(5, meanx4, sdx4, log = TRUE))
  
  N <- 3
  y      <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z      <- rep(2L, N * N)
  
  got <- prior_combined(params, 2, y, x_vars, z, FALSE, upri)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("Missing user_fixed_priors component raises an informative error", {
  params <- rep(0, 5)
  N <- 2
  y <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z <- rep(1L, N * N)
  
  expect_error(
    prior_combined(params, 2, y, x_vars, z, FALSE, user_fixed_priors = NULL),
    "user_fixed_priors"
  )
  
  upri <- list(component1 = list(
    meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1
  ))
  expect_error(
    prior_combined(params, 2, y, x_vars, z, FALSE, upri),
    "component2"
  )
})

test_that("Empty component is handled via fallback (not a crash)", {
  N <- 3
  y <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z <- rep(1L, N * N)                    # no 2s at all
  params <- rep(0, 5)
  
  expect_silent(val <- prior_combined(params, 2, y, x_vars, z, TRUE, NULL))
  expect_true(is.finite(val))
})