#' @title Size Prior for the Negative Binomial-Type Distributions
#'
#' @description
#' This function computes the prior contribution of the \code{size} parameter (also known as the dispersion parameter) 
#' for a specified component in models using a Negative Binomial (NB) or Zero-Inflated Negative Binomial (ZINB) distribution. 
#' The prior is assumed to follow a Gamma distribution. 
#'
#' @usage
#' size_prior(size_value, component)
#'
#' @param size_value A numeric value representing the \code{size} (dispersion) parameter of the NB or ZINB distribution.
#'   This parameter controls the variance of the distribution, with larger values implying less overdispersion.
#'
#' @param component An integer specifying which component of the mixture model to consider (e.g., \code{1}, \code{2}, or \code{3}). 
#'   Different components can have different prior distributions for the size parameter.
#'
#' @details
#' In NB and ZINB distributions, the \code{size} parameter controls the variance relative to the mean. 
#' Bayesian inference often places a prior on this parameter to regularize its estimation. 
#' By assigning a Gamma prior, we leverage its conjugacy properties and flexibility in capturing a range of plausible 
#' dispersion values.
#'
#' This function assigns a component-specific Gamma prior. Each component can have distinct shape (and possibly rate) parameters, 
#' allowing the model to accommodate different dispersion characteristics across mixture components.
#'
#' The returned value is the log-density of the Gamma prior at the given \code{size_value}.
#'
#' @return
#' A numeric value representing the log of the prior density for the \code{size_value} parameter under the Gamma distribution 
#' defined for the specified component.
#'
#' @examples
#' # Example: Compute the size prior for size_value = 2.5 in component 1
#' log_prior_comp1 <- size_prior(size_value = 2.5, component = 1)
#' log_prior_comp1
#'
#' # Example: Compute the size prior for size_value = 10 in component 3
#' log_prior_comp3 <- size_prior(size_value = 10, component = 3)
#' log_prior_comp3
#'
#' @export
#' 
# Gamma prior for size parameters
size_prior <- function(size_value, component) {
  if (component < 1 || component > 3) {
    stop("Invalid component specified. Must be 1, 2, or 3.")
  }
  
  # Define parameters for each component
  shape_params <- c(7.5, 2, 5)
  rate_param <- 1

  # Compute log-probability using the appropriate shape parameter
  dgamma(size_value, shape = shape_params[component], rate = rate_param, log = TRUE)
}


