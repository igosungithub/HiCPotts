test_that("pred_combined implements a + b*log1p(x1) + ... exactly", {
  N <- 3
  params <- c(0.5, 1.0, 2.0, -1.0, 0.25)
  z <- matrix(1, N, N)
  x1 <- matrix(seq_len(N * N), N, N)
  x2 <- matrix(seq_len(N * N) * 2, N, N)
  x3 <- matrix(seq_len(N * N) * 3, N, N)
  x4 <- matrix(seq_len(N * N) * 4, N, N)
  x_vars <- list(list(x1), list(x2), list(x3), list(x4))
  
  got <- pred_combined(params, z, x_vars, 1, N)
  expected <- params[1] +
    params[2] * log1p(as.vector(x1)) +
    params[3] * log1p(as.vector(x2)) +
    params[4] * log1p(as.vector(x3)) +
    params[5] * log1p(as.vector(x4))
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("Empty component returns numeric(0)", {
  N <- 2
  z <- matrix(1, N, N)
  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  expect_equal(pred_combined(params, z, x_vars, 2, N), numeric(0))
})

test_that("Input validation fires correctly", {
  N <- 3
  params <- rep(0, 5)
  z <- matrix(1, N, N)
  
  expect_error(pred_combined(params, z, NULL, 1, N), "x_vars cannot be NULL")
  expect_error(pred_combined(params, z, list(matrix(0, N, N)), 1, N),
               "list of four")
  x_vars_bad <- replicate(4, list(matrix(0, 2, 2)), simplify = FALSE)  # wrong N
  expect_error(pred_combined(params, z, x_vars_bad, 1, N),
               "N x N")
  expect_error(pred_combined(1:3, z,
                             replicate(4, list(matrix(0,N,N)), simplify = FALSE),
                             1, N),
               "length 5")
})