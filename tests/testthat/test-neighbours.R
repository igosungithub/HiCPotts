## Tests for Neighbours_combined --------------------------------------
## Regression test: the function should count matching neighbours in {0..4}.

test_that("uniform lattice: every interior cell has 4 matching neighbours", {
  N <- 5
  z <- matrix(1, N, N)                  # all same state
  nb <- Neighbours_combined(z, N)
  # corners have 2 neighbours on the grid, edges have 3, interior has 4
  expected <- matrix(4, N, N)
  expected[1, ]            <- 3
  expected[N, ]            <- 3
  expected[, 1]            <- 3
  expected[, N]            <- 3
  expected[1, 1] <- expected[1, N] <- expected[N, 1] <- expected[N, N] <- 2
  expect_equal(nb, expected)
})

test_that("checkerboard lattice: every cell has 0 matching neighbours", {
  N <- 4
  z <- outer(1:N, 1:N, function(i, j) ifelse((i + j) %% 2 == 0, 1, 2))
  nb <- Neighbours_combined(z, N)
  expect_equal(nb, matrix(0, N, N))
})

test_that("known 3x3 configuration: counts match hand calculation", {
  ## Layout:    1 2 1
  ##            2 1 2
  ##            1 2 1
  ## Centre (2,2) has neighbours {up=2, down=2, left=2, right=2}, ref=1 -> 0
  ## Corner (1,1) has neighbours {down=2, right=2}, ref=1 -> 0
  ## Edge cell (1,2) has neighbours {down=1, left=1, right=1}, ref=2 -> 0
  N <- 3
  z <- matrix(c(1,2,1, 2,1,2, 1,2,1), N, N, byrow = TRUE)
  nb <- Neighbours_combined(z, N)
  expect_equal(nb, matrix(0, N, N))
})

test_that("all counts are in {0, 1, 2, 3, 4} for random lattices", {
  set.seed(1L)
  for (rep in 1:10) {
    N <- sample(3:8, 1)
    z <- matrix(sample(1:3, N * N, replace = TRUE), N, N)
    nb <- Neighbours_combined(z, N)
    expect_true(all(nb %in% 0:4),
                info = sprintf("rep %d failed: max = %g, min = %g",
                               rep, max(nb), min(nb)))
  }
})

test_that("proposed_value mode: counts reference matches of proposed against current", {
  ## current:   1 1           proposed:   1 2
  ##            1 1                       1 1
  ## At (1,2): ref = 2, neighbours (down=1, left=1) -> 0 matches
  ## At (1,1): ref = 1, neighbours (down=1, right=1) -> 2 matches
  ## At (2,1): ref = 1, neighbours (up=1,   right=1) -> 2 matches
  ## At (2,2): ref = 1, neighbours (up=2,   left=1)  -> 1 match
  current  <- matrix(1, 2, 2)
  proposed <- matrix(c(1, 2, 1, 1), 2, 2, byrow = TRUE)
  nb <- Neighbours_combined(current, 2, proposed_value = proposed)
  expect_equal(nb, matrix(c(2, 0, 2, 2), 2, 2, byrow = TRUE))
})