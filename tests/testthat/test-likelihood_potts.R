test_that("Returns an NxN matrix", {
  for (N in c(2, 3, 5)) {
    x <- rep(0.5, N * N)
    n <- rnorm(N * N)
    m <- likelihood_gamma(x, n, N)
    expect_equal(dim(m), c(N, N))
    expect_true(is.numeric(m))
    expect_true(all(is.finite(m)))
  }
})

test_that("Element-wise formula: m[i,j] == x[i,j] * pair_neighbours[i,j]", {
  N <- 3
  x <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  n <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  expected <- matrix(x * n, nrow = N, ncol = N)
  got <- likelihood_gamma(x, n, N)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("Scalar x broadcasts to the lattice", {
  N <- 4
  n <- sample(0:4, N * N, replace = TRUE)
  m <- likelihood_gamma(0.3, n, N)
  expect_equal(dim(m), c(N, N))
  expect_equal(as.vector(m), 0.3 * n, tolerance = 1e-12)
})

test_that("Agrees with the formal Potts log-potential gamma * neighbours", {
  N <- 5
  set.seed(1L)
  gamma_val <- 0.7
  neigh     <- matrix(sample(0:4, N * N, replace = TRUE), N, N)
  m <- likelihood_gamma(gamma_val, as.vector(neigh), N)
  expect_equal(as.vector(m), gamma_val * as.vector(neigh),
               tolerance = 1e-12)
})

test_that("Zero neighbours produce zero log-potential", {
  N <- 3
  m <- likelihood_gamma(5.0, rep(0, N * N), N)
  expect_equal(m, matrix(0, N, N), tolerance = 1e-12)
})

test_that("Gamma = 0 produces zero log-potential regardless of neighbours", {
  N <- 3
  n <- sample(0:4, N * N, replace = TRUE)
  m <- likelihood_gamma(0, n, N)
  expect_equal(m, matrix(0, N, N), tolerance = 1e-12)
})

test_that("Numerically stable for large gamma (no softmax overflow)", {
  N <- 3
  n <- sample(0:4, N * N, replace = TRUE)
  expect_silent(m <- likelihood_gamma(300, n, N))
  expect_true(all(is.finite(m)))
  expect_equal(as.vector(m), 300 * n, tolerance = 1e-12)
})

test_that("Length mismatch and bad inputs raise informative errors", {
  expect_error(likelihood_gamma(1:3, 1:4, 2),
               "same length")
  expect_error(likelihood_gamma(1:3, 1:3, 2),
               "N\\*N")
  expect_error(likelihood_gamma(rep(0.5, 4), rep(0, 4), 0),
               "positive integer")
})