#' @title Likelihood Function for a Specified Model Component and Distribution
#'
#' @description
#' This function computes the likelihood contribution of a specific model component using a chosen distribution. 
#' It is designed to work within a mixture-modeling or hierarchical modeling framework, where each interaction 
#' (such as genomic interactions in a Hi-C experiment) can be modeled by a combination of regression parameters and 
#' potentially zero-inflation or overdispersion parameters.
#'
#' @usage likelihood_combined(params, z, x_vars, component, N, pred_combined, y, dist, theta = NULL, size = NULL)
#'
#' @param params A numeric vector of parameters associated with the model’s linear predictors. 
#'   Typically includes intercepts and regression coefficients related to the covariates in \code{x_vars}.
#'
#' @param z A matrix or array representing the latent state indicators or other structural variables 
#'   that influence the model. \code{z} helps map each observation to a particular mixture component.
#'
#' @param x_vars A list of covariates used as predictors in the linear model. Each element in the list 
#'   corresponds to a covariate vector, and all must be of length \code{N}, the number of observations.
#'
#' @param component An integer or factor specifying which component of the mixture model is currently 
#'   being evaluated. In the 3-component mixture, this is \code{1}, \code{2}, or \code{3}.
#'
#' @param N An integer specifying the number of observations. This 
#'   should match the length of the response variable \code{y} and the covariates in \code{x_vars}.
#'
#' @param pred_combined A numeric vector of predictor values obtained from a prediction function. This 
#'   typically represents the linear predictor \eqn{\lambda} for the given component.
#'
#' @param y A numeric vector of observed interaction counts (the response variable). Each element 
#'   corresponds to one observation (e.g., interaction count between a pair of genomic loci).
#'
#' @param dist A character string specifying the distribution to be used for modeling the interaction counts. 
#'   Options may include:
#'   \itemize{
#'     \item \code{"Poisson"}: Poisson distribution
#'     \item \code{"NB"}: Negative Binomial distribution
#'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
#'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
#'   }
#'
#' @param theta (Optional) A numeric value for the zero-inflation parameter. Required for ZIP and ZINB models. 
#'   This parameter controls the probability of excess zeros not explained by the Poisson or NB component.
#'
#' @param size (Optional) A numeric value for the size (overdispersion) parameter. Required for NB and ZINB models. 
#'   This parameter captures variance that exceeds that of a Poisson distribution.
#'
#' @details
#' This function calculates the likelihood for a single component given a set of parameters and covariates. 
#' The steps typically involved are:
#' \enumerate{
#'   \item Compute the linear predictor \eqn{\lambda} from the supplied parameters and covariates. 
#'         The \code{pred_combined} represents \eqn{\log(\lambda)}.
#'   \item Depending on \code{dist}, convert the linear predictor into a mean parameter (e.g., \eqn{\lambda} for Poisson or NB).
#'   \item Compute the likelihood of each observed count \code{y[i]} under the chosen distribution with 
#'         the given parameters (\eqn{\lambda}, \eqn{\theta}, \eqn{size}, etc.).
#'   \item For zero-inflated models (ZIP, ZINB), the likelihood incorporates the probability of an extra zero.
#'   \item Return the computed likelihood values. This is used internally to 
#'         evaluate and update model parameters during the MCMC steps.
#' }
#'
#' The function is a building block in a larger modeling framework (MCMC inference for mixture models). 
#' While end-users might not call it directly, it enables flexible specification of distributions and 
#' model components for complex hierarchical models.
#'
#' @return
#' The function returns the computed likelihood under the specified model setup. 
#' It returns a numeric vector of likelihood contributions for each observation.
#'
#' @examples
#' \dontrun{
#' # Example setup
#' N <- 100
#' y <- rpois(N, lambda = 5)  # synthetic response data
#' x_vars <- list(dist = runif(N, 0, 1))
#' params <- c(intercept = 1.5, beta = 0.3)
#' z <- matrix(1, nrow = N, ncol = 3)  # placeholder latent structure
#' pred_combined <- exp(params["intercept"] + params["beta"] * x_vars$dist)
#'
#' # Compute likelihood under a Poisson model, component 1
#' ll_values <- likelihood_combined(
#' pred_combined = pred_combined,
#'   params = params,
#'   z = z,
#'   y = y,
#'   x_vars = x_vars,
#'   component = 1,
#'   theta,
#'   size,
#'   N = N,
#'   dist = "Poisson"
#' )
#'
#' sum(ll_values) # sum of likelihood contributions
#' }
#'
#' @seealso
#' \code{\link{dpois}}, \code{\link{dnbinom}} for related probability mass functions.
#'
#' @importFrom stats dgamma dnbinom dnorm dpois rbeta rgamma rnorm
#' 
#' 
#' 
#' @export
#'
likelihood_combined <- function(pred_combined, params, z, y, x_vars, component, theta, size, N, dist) {
  # Subset the data based on the component
  yc = y[z == component]
  lambda=pred_combined(params, z, x_vars, component, N)
  # Calculate the likelihood based on the specified distribution
  if (component == 1) {
    if (dist == "ZIP") {
      singlelikelihoods = ifelse(yc == 0,
                                 log(theta + (1 - theta) * exp(-lambda)),
                                 log(1 - theta) + dpois(yc, lambda = exp(lambda), log = TRUE))
    } else if (dist == "Poisson") {
      singlelikelihoods = dpois(yc, lambda = exp(lambda), log = TRUE)
    } else if (dist == "NB") {
      singlelikelihoods = dnbinom(yc, size = size, mu = exp(lambda), log = TRUE)
    } else if (dist == "ZINB") {
      singlelikelihoods = ifelse(yc == 0,
                                 log(theta + (1 - theta) * exp(-lambda)),
                                 log(1 - theta) + dnbinom(yc, size = size, mu = exp(lambda), log = TRUE))
    } else {
      stop("Invalid distribution specified.")
    }
  } else if(component==2 || component==3) {
    # For components other than 1, restrict to Poisson and NB distributions
    if (dist == "Poisson" || dist == "ZIP") {
      singlelikelihoods = dpois(yc[component], lambda = exp(lambda)[component], log = TRUE)
    } else if (dist == "NB" || dist == "ZINB") {
      singlelikelihoods = dnbinom(yc[component], size = size, mu = exp(lambda)[component], log = TRUE)
    } else {
      stop("Invalid distribution specified for component > 1.")
    }
  }
  sumll = sum(singlelikelihoods)
  return(sumll)
}
