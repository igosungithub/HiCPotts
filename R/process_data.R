#' @title CSV file Data Processing for Hi-C Interaction Matrices and Covariates
#'
#' @description
#' The \code{process_data} function takes a data frame with interaction information and associated covariates
#' (genomic bins, GC content, transposable elements(TES), and accessibility) and converts it into
#' a structured list of \eqn{N \times N} matrices suitable for modeling. It also provides options for scaling interaction counts and checking for
#' required columns and missing values.
#'
#' @usage
#' process_data(data, N, scale_max = 500, standardization_y = TRUE)
#'
#' @param data A \code{data.frame} containing the following required columns:
#'   \itemize{
#'     \item \code{start}: Numeric start coordinates for locus i.
#'     \item \code{end}: Numeric end coordinates for locus j.
#'     \item \code{interactions}: Numeric vector of observed interaction counts between loci i and j.
#'     \item \code{GC}: GC content measure for the given loci.
#'     \item \code{ACC}: Accessibility score for the given loci.
#'     \item \code{TES}: A measure related to transposable elements for the given loci.
#'   }
#'
#' @param N An integer specifying the dimension of the resulting \eqn{N \times N} matrices. The data
#'   provided should correspond to \eqn{N^2} interactions or multiples thereof, as it will be reshaped
#'   into one or more \eqn{N \times N} matrices.
#'
#' @param scale_max A numeric value indicating the maximum scaling factor for interaction counts
#'   when \code{standardization_y = TRUE}. Defaults to 500. After scaling, interaction counts
#'   range between [1, \code{scale_max}].
#'
#' @param standardization_y A logical value. If \code{TRUE}, interaction counts are scaled to the
#'   range [1, \code{scale_max}] by applying a min-max normalization and rounding. If \code{FALSE},
#'   the raw interaction counts are used as provided.
#'
#' @details
#' This function processes a long-format data frame where each row represents an interaction between
#' two loci (denoted by \code{start} and \code{end}) and associated covariate values.
#' Key steps include:
#' \enumerate{
#'   \item Validating that required columns are present and checking for missing values.
#'   \item Optionally scaling the interaction counts to a fixed range to reduce skew or prepare data
#'         for models sensitive to the magnitude of counts.
#'   \item Computing genomic distances as \eqn{|end - start|}.
#'   \item Reshaping the vectorized interaction and covariate data into \eqn{N \times N} matrices.
#'         If the number of rows in \code{data} is greater than \eqn{N^2}, it splits them into multiple
#'         \eqn{N \times N} matrices (one for each segment of length \eqn{N^2}).
#'   \item Storing the resulting data in a list structure suitable for downstream analyses.
#' }
#'
#' The returned \code{x_vars} list contains multiple \eqn{N \times N} matrices for each covariate:
#' \itemize{
#'   \item \code{x_vars$distance}: A list of one or more \eqn{N \times N} matrices containing genomic distances.
#'   \item \code{x_vars$GC}: A list of \eqn{N \times N} matrices for GC content.
#'   \item \code{x_vars$TES}: A list of \eqn{N \times N} matrices for TES values.
#'   \item \code{x_vars$ACC}: A list of \eqn{N \times N} matrices for accessibility scores.
#' }
#'
#' The \code{y} element in the returned list contains the scaled (or raw) interaction counts structured
#' as a list of \eqn{N \times N} matrices.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{\code{x_vars}}{A named list containing lists of \eqn{N \times N} matrices for each covariate:
#'   \code{distance}, \code{GC}, \code{TES}, \code{ACC}.}
#'
#'   \item{\code{y}}{A list of \eqn{N \times N} matrices representing the (scaled) interaction counts
#'   corresponding to the covariates in \code{x_vars}.}
#' }
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   start = rep(1:10, each = 10),
#'   end = rep(11:20, times = 10),
#'   interactions = rpois(100, 5),
#'   GC = runif(100, 0, 1),
#'   TES = runif(100, 0, 1),
#'   ACC = runif(100, 0, 1)
#' )
#' processed <- process_data(df, N = 10, scale_max = 500, standardization_y = TRUE)
#' x_vars <- processed$x_vars
#' y_matrices <- processed$y
#' str(x_vars) # Show structure of covariates
#' str(y_matrices) # Show structure of interaction matrices
#'
#' #\donttest{
#' # Extended example with larger dataset
#' # Suppose we have a data frame 'large_df' corresponding to a 20x20 interaction matrix
#' large_df <- data.frame(
#'   start = rep(1:20, each = 20),
#'   end = rep(21:40, times = 20),
#'   interactions = rpois(400, 5),
#'   GC = runif(400, 0, 1),
#'   TES = runif(400, 0, 1),
#'   ACC = runif(400, 0, 1)
#' )
#' processed <- process_data(large_df, N = 20, scale_max = 500, standardization_y = TRUE)
#' x_vars <- processed[[1]]
#' y_matrices <- processed[[2]]
#' str(x_vars)
#' str(y_matrices)
#' # See vignette("HMRFHiC_vignette") for detailed examples with real Hi-C data.
#' #}
#'
#' @seealso
#' \code{\link{run_chain_betas}} for downstream MCMC inference
#' and modeling steps once data is processed into matrix form.
#'
#' @export
#
#
process_data <- function(data, N, scale_max = 500, standardization_y = TRUE) {
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
    scaled_data <- round((y_dat - min_val) / (max_val - min_val) * (scale_max - 1) + 1) # Scaling to range [1, 1000]
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
      length_diff <- max(N^2 - length(x), 0)
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
