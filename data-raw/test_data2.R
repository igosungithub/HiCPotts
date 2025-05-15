## code to prepare `DATASET` dataset goes here

# Load required libraries
library(usethis)
library(readr)

# Path to the raw data file (adjust the path if needed)
data_file <- "inst/extdata/test_data.csv"

# Read the raw data
raw_data <- read_csv(data_file)
N=36

DATASET <- function(raw_data, N, scale_max=500, standardization_y = TRUE) {
  # Check required columns in the data
  required_columns <- c("start", "end", "interactions", "GC", "ACC", "TES")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the data:", paste(missing_columns, collapse = ", ")))
  }
  
  # Check for NA values in the data
  if (anyNA(data)) {
    stop("The data contains NA values. Please handle missing data before processing.")
  }
  
  # Process interactions
  y_dat <- data$interactions
  if (standardization_y) {
    min_val <- min(y_dat)
    max_val <- max(y_dat)
    scaled_data <- round((y_dat - min_val) / (max_val - min_val) * (scale_max - 1) + 1)  # Scaling to range [1, 1000]
  } else {
    scaled_data <- y_dat
  }
  
  # Process independent variables
  x111 <- abs(data$end - data$start)
  x222 <- data$GC
  x333 <- data$TES
  x444 <- data$ACC
  
  # Convert data to symmetric matrices
  as <- split(x111, ceiling(seq_along(x111) / (N * N)))
  ab <- split(x222, ceiling(seq_along(x222) / (N * N)))
  ac <- split(x333, ceiling(seq_along(x333) / (N * N)))
  ad <- split(x444, ceiling(seq_along(x444) / (N * N)))
  yy <- split(scaled_data, ceiling(seq_along(scaled_data) / (N * N)))
  
  
  expand <- function(ss, N) {
    lapply(ss, function(x) {
      length_diff = max(N^2 - length(x), 0)
      matrix(c(x, rep(0, length_diff)), N, N)
    })
  }
  
  
  
  ass <- expand(as, N)
  abb <- expand(ab, N)
  acc <- expand(ac, N)
  add <- expand(ad, N)
  y1 <- expand(yy, N)
  
  # Prepare outputs
  x_vars <- list(
    distance = ass,
    GC = abb,
    TES = acc,
    ACC = add
  )
  y_sim1 <- y1
  
  # Check for NAs in the output
  if (anyNA(x_vars)) {
    warning("NA values found in x_vars.")
  }
  if (anyNA(y_sim1)) {
    warning("NA values found in y_sim1.")
  }
  
  # Return outputs
  list(x_vars = x_vars, y = y_sim1)
}


# Save the processed dataset in the `data/` folder
usethis::use_data(DATASET, overwrite = TRUE)
