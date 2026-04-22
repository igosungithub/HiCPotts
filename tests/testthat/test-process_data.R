test_that("Returns N x N matrices in the expected shape", {
  set.seed(1L)
  N <- 4
  df <- data.frame(
    start        = rep(seq_len(N),    each = N),
    end          = rep(seq_len(N) * 2, N),
    interactions = rpois(N * N, 5),
    GC           = runif(N * N),
    TES          = rpois(N * N, 1),
    ACC          = runif(N * N)
  )
  out <- process_data(df, N = N, standardization_y = FALSE)
  expect_named(out, c("x_vars", "y"))
  expect_named(out$x_vars, c("distance","GC","TES","ACC"))
  expect_equal(dim(out$x_vars$distance[[1]]), c(N, N))
  expect_equal(dim(out$y[[1]]), c(N, N))
})

test_that("Raises when nrow(data) is not a multiple of N^2 by default (regression)", {
  df <- data.frame(
    start        = 1:3,    # only 3 rows
    end          = 4:6,
    interactions = 1:3,
    GC           = c(.1,.2,.3),
    TES          = 1:3,
    ACC          = c(.1,.2,.3)
  )
  expect_error(process_data(df, N = 2), "not a multiple")
})

test_that("pad_with_zero=TRUE zero-pads", {
  df <- data.frame(
    start        = 1:3, end = 4:6, interactions = 1:3,
    GC = c(.1,.2,.3), TES = 1:3, ACC = c(.1,.2,.3)
  )
  out <- process_data(df, N = 2, pad_with_zero = TRUE)
  expect_equal(dim(out$y[[1]]), c(2, 2))
})

test_that("standardization_y emits a warning and stays in [1, scale_max]", {
  N <- 3
  df <- data.frame(
    start        = rep(1:N, each = N),
    end          = rep(1:N * 2, N),
    interactions = 1:(N * N),
    GC = runif(N * N), TES = rep(1, N * N), ACC = runif(N * N)
  )
  expect_warning(
    out <- process_data(df, N = N,
                        scale_max = 100,
                        standardization_y = TRUE),
    "min-max"
  )
  vals <- as.vector(out$y[[1]])
  expect_true(all(vals >= 1 & vals <= 100))
})

test_that("Missing columns error with required-column message", {
  df <- data.frame(start = 1, end = 1)
  expect_error(process_data(df, N = 1), "required columns")
})