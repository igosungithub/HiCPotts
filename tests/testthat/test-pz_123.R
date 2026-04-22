test_that("pz_123 returns an NxN finite matrix for a ZIP configuration", {
  set.seed(3L)
  N <- 4
  z <- matrix(sample(1:3, N * N, replace = TRUE), N, N)
  sn <- matrix(sample(0:4, N * N, replace = TRUE), N, N)
  y <- matrix(rpois(N * N, 3), N, N)
  pred <- function(params, z, x_vars, component, N) {
    rep(log(3), sum(z == component))
  }
  chains <- list(matrix(0, nrow = 5, ncol = 5),
                 matrix(0, nrow = 5, ncol = 5),
                 matrix(0, nrow = 5, ncol = 5))
  chain_gamma <- rep(0.3, 20)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  names(x_vars) <- c("distance","GC","TES","ACC")
  size_chain <- matrix(5, 3, 20)
  
  got <- pz_123(z = z, sum_neighbours = sn, y = y,
                pred_combined = pred, chains = chains,
                chain_gamma = chain_gamma, x_vars = x_vars,
                theta = 0.25, size_chain = size_chain,
                N = N, iter = 1, dist = "ZIP")
  expect_equal(dim(got), c(N, N))
  expect_true(all(is.finite(got)))
})

test_that("ZIP zero cell uses exp(-mu), not exp(-eta) (regression test)", {
  ## Build a 1-cell lattice (N=1). z=1 so ZIP zero branch applies.
  ## With eta = log(3) and theta = 0.2 the unweighted probability should be
  ##   theta + (1-theta)*exp(-3)  = 0.2 + 0.8 * exp(-3)
  ## NOT
  ##   theta + (1-theta)*exp(-log(3)) = 0.2 + 0.8/3
  N <- 1
  z <- matrix(1, 1, 1)
  sn <- matrix(0, 1, 1)            # exp(gamma * 0) = 1
  y <- matrix(0, 1, 1)
  pred <- function(params, z, x_vars, component, N) rep(log(3), 1)
  chains <- list(matrix(0, nrow = 5, ncol = 5),
                 matrix(0, nrow = 5, ncol = 5),
                 matrix(0, nrow = 5, ncol = 5))
  chain_gamma <- rep(0.3, 5)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  names(x_vars) <- c("distance","GC","TES","ACC")
  size_chain <- matrix(5, 3, 5)
  theta <- 0.2
  
  got <- as.vector(pz_123(z, sn, y, pred, chains, chain_gamma, x_vars,
                          theta, size_chain, N, iter = 1, dist = "ZIP"))
  expected <- log(theta + (1 - theta) * exp(-3))
  wrong    <- log(theta + (1 - theta) * exp(-log(3)))
  expect_equal(got, expected, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(got, wrong)))
})