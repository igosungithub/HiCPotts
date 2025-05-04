#' @title Prediction Function for Combined Model
#'
#' @description
#' This function computes the predicted linear predictor (\eqn{\lambda}) values for a specified model component.
#' It is designed to correct for various sources of bias by incorporating multiple covariates through a 
#' log-transformed relationship. The resulting predictions serves as inputs to downstream likelihood 
#' or posterior calculations in the Bayesian framework.
#'
#' @usage
#' pred_combined(params, z, x_vars, component, N)
#'
#' @param params A numeric vector of parameters for the model. Typically, \code{params} = \eqn{(a, b, c, d, e)}, 
#'   where \eqn{a} is the intercept and \eqn{b, c, d, e} are coefficients for the log-transformed covariates.
#'
#' @param z A matrix used to indicate component membership of each observation. 
#'   The function uses \code{z == component} to select the subset of observations belonging to the specified component.
#'
#' @param x_vars A list of covariates, each element itself a list or vector of length \code{N}, 
#'   corresponding to different sources of bias or predictors. For example:
#'   \itemize{
#'     \item \code{x_vars[[1]]}: The first covariate (wrapped in a list), accessible via \code{x_vars[[1]][[1]]}.
#'     \item \code{x_vars[[2]]}, \code{x_vars[[3]]}, \code{x_vars[[4]]}: Additional covariates in similar nested-list structures.
#'   }
#'   Each covariate vector must be of length \code{N}, and the values are expected to be non-negative since a \eqn{\log(x+1)} transformation is applied.
#'
#' @param component An integer specifying which component of the mixture or hierarchical model to use 
#'   when subsetting the data. For the three components, \code{component = 1} 
#'   will extract the subset of data where \code{z == 1}.
#'
#' @param N An integer specifying the total number of observations. This should match the length of the 
#'   \code{y} vector and each covariate vector.
#'
#' @details
#' The predicted values (\eqn{\lambda}) are computed as:
#'
#' \deqn{\lambda_i = a + b \cdot \log(x1_i + 1) + c \cdot \log(x2_i + 1) + d \cdot \log(x3_i + 1) + e \cdot \log(x4_i + 1)}
#'
#' for each observation \eqn{i} in the specified component. This log-transform often stabilizes variance and 
#' normalizes skewed distributions of covariates.
#'
#' The resulting \eqn{\lambda} values serves as linear predictors (in a Poisson, NB, ZIP, or ZINB model) 
#' and as inputs into the likelihood and posterior computation steps.
#'
#' @return
#' A numeric vector containing the predicted values (\eqn{\lambda}) for observations belonging to the specified component.
#'
#' @examples
#' \dontrun{
#' # Assume we have:
#' # N = 100 observations
#' # z indicates component membership for each observation
#' # x_vars is a list of four covariates, each nested in a list: x_vars[[i]][[1]]
#' N <- 100
#' z <- sample(1:3, N, replace = TRUE)
#' x_vars <- list(
#'   list(runif(N, 0, 10)), # first covariate
#'   list(runif(N, 0, 5)),  # second covariate
#'   list(runif(N, 1, 20)), # third covariate
#'   list(runif(N, 0, 2))   # fourth covariate
#' )
#' params <- c(a = 0.5, b = 0.1, c = -0.05, d = 0.2, e = 0.3)
#'
#' # Get predicted values for component 1:
#' preds_comp1 <- pred_combined(params, z, x_vars, component = 1, N = N)
#' head(preds_comp1)
#' }
#'
#' @export
#'
pred_combined <- function(params, z, x_vars, component, N) {
  
  # Validate inputs
  if (is.null(x_vars)) {
    stop("x_vars cannot be NULL")
  }
  if (!is.list(x_vars) || length(x_vars) != 4) {
    stop("x_vars must be a list of four matrices")
  }
  if (any(sapply(x_vars, function(x) !is.matrix(x[[1]]) || dim(x[[1]])[1] != N || dim(x[[1]])[2] != N))) {
    stop("Each x_vars matrix must have dimensions N x N")
  }
  
  # Extract parameters
  a <- params[1]
  b <- params[2]
  c <- params[3]
  d <- params[4]
  e <- params[5]

  # Create logical mask and subset matrices
  logical_mask <- z == component
  x1_sub <- x_vars[[1]][[1]][logical_mask]
  x2_sub <- x_vars[[2]][[1]][logical_mask]
  x3_sub <- x_vars[[3]][[1]][logical_mask]
  x4_sub <- x_vars[[4]][[1]][logical_mask]

  # Compute pred for any component
  pred <- a + b * log(x1_sub+1) + c * log(x2_sub+1) + d * log(x3_sub+1) + e * log(x4_sub+1)
  return(pred)
}


