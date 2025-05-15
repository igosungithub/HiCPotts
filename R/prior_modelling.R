#' @title Prior Function for Model Components
#'
#' @description
#' This function computes the prior contribution to the log-posterior for the parameters associated with a specified model component. 
#' It can use either data-driven priors (estimated from the subset of data corresponding to that component) or user-specified fixed priors. 
#' By integrating these priors into the Bayesian inference framework, the function helps ensure that parameter estimates remain 
#' grounded in reasonable values and incorporate domain knowledge or data-derived distributions.
#'
#' @usage
#' prior_combined(params, component, y, x_vars, z, use_data_priors = TRUE, user_fixed_priors)
#'
#' @param params A numeric vector of parameters for the model. Typically, \code{params} will include a set of parameters 
#'   associated with the regression coefficients for the model. For example, \code{params} might be 
#'   \eqn{(a, b, c, d, e)} representing coefficients for the covariates in \code{x_vars}.
#'
#' @param component An integer specifying which component of the model to consider. The model is a mixture with multiple 
#'   components, this argument indicates which component’s parameters and priors are being evaluated.
#'
#' @param y A numeric vector of observed responses (e.g., interaction counts). This is used when estimating data-driven 
#'   priors. The vector should be of length \code{N}, where \code{N} is the size of the system.
#'
#' @param x_vars A list of covariate information. Each element of \code{x_vars} is typically another list or vector 
#'   representing a single covariate. These covariates are used to compute data-driven priors by extracting summary statistics
#'   from the data corresponding to the specified component.
#'
#' @param z A matrix or similar structure used to assign observations to components. If \code{z == component}, that observation 
#'   belongs to the given component. The function uses this to subset \code{y} and \code{x_vars} for calculating data-driven priors.
#'
#' @param use_data_priors A logical value indicating whether to use data-driven priors (\code{TRUE}) or user-specified fixed priors (\code{FALSE}). 
#'   When \code{TRUE}, the function computes means and standard deviations from the subset of the data assigned to the specified component.
#'
#' @param user_fixed_priors A list of user-defined priors for each component. This argument is only used when 
#'   \code{use_data_priors = FALSE}. The list should contain entries named \code{component1}, \code{component2}, and \code{component3},
#'   each defining \code{meany}, \code{meanx1}, \code{meanx2}, \code{meanx3}, \code{meanx4}, 
#'   and corresponding \code{sdy}, \code{sdx1}, \code{sdx2}, \code{sdx3}, \code{sdx4}.
#'
#' @details
#' This function supports two modes:
#' \enumerate{
#'   \item \strong{Data-Driven Priors:} When \code{use_data_priors = TRUE}, the function subsets the data for the specified \code{component}, 
#'         calculates summary statistics, and then draws from those statistics to form priors for each parameter. 
#'         This approach allows the priors to adapt to the characteristics of the data associated with that component.
#'
#'   \item \strong{User-Specified Fixed Priors:} When \code{use_data_priors = FALSE}, the function uses priors provided by the user 
#'         through \code{user_fixed_priors}, bypassing any data-driven estimation.
#' }
#'
#' In both modes, the priors are assumed to be normal distributions for each parameter. The function returns the sum of the log of these priors.
#'
#' @return
#' A single numeric value representing the sum of the log-priors for all parameters in \code{params}. 
#' This value can be integrated into a Bayesian inference framework to update parameter estimates based on the chosen priors.
#'
#' @examples
#' \dontrun{
#' # Example setup:
#' N <- 100
#' y <- rpois(N, lambda = 5)
#' x_vars <- list(
#'   list(runif(N)),  # x1
#'   list(rnorm(N)),  # x2
#'   list(rnorm(N)),  # x3
#'   list(rnorm(N))   # x4
#' )
#' z <- sample(1:3, N, replace = TRUE)
#' params <- c(a = 1, b = 0.5, c = -0.2, d = 2, e = 0.3)
#'
#' # Using data-driven priors for component 1
#' prior_val_data <- prior_combined(
#'   params = params,
#'   component = 1,
#'   y = y,
#'   x_vars = x_vars,
#'   z = z,
#'   use_data_priors = TRUE,
#'   user_fixed_priors = NULL
#' )
#'
#' # Using user-specified fixed priors for component 2
#' user_fixed_priors <- list(
#'   component1 = list(meany=5, meanx1=0, meanx2=0, meanx3=0, meanx4=0,
#'   sdy=1, sdx1=1, sdx2=1, sdx3=1, sdx4=1),\cr
#'   component2 = list(meany=10, meanx1=1, meanx2=1, meanx3=1, meanx4=1,
#'   sdy=2, sdx1=2, sdx2=2, sdx3=2, sdx4=2),\cr
#'   component3 = list(meany=2, meanx1=2, meanx2=2, meanx3=2, meanx4=2,
#'   sdy=1, sdx1=1, sdx2=1, sdx3=1, sdx4=1)
#' )
#'
#' prior_val_fixed <- prior_combined(
#'   params = params,
#'   component = 2,
#'   y = y,
#'   x_vars = x_vars,
#'   z = z,
#'   use_data_priors = FALSE,
#'   user_fixed_priors = user_fixed_priors
#' )
#'
#' }
#'
#' @export
#' 
#prior for the sources of biases
prior_combined <- function(params, component, y, x_vars, z, use_data_priors = TRUE, user_fixed_priors) {
  # Extract parameters from the 'params' vector
  a = params[1]
  b = params[2]
  c = params[3]
  d = params[4]
  e = params[5]

  # Set a small positive value to avoid using zero or negative standard deviations
  epsilon = 1e-6

  # Subset data based on component if using data-driven priors
  if (use_data_priors) {
    logical_mask = z == component
    x1 = x_vars[[1]][[1]][logical_mask]
    x2 = x_vars[[2]][[1]][logical_mask]
    x3 = x_vars[[3]][[1]][logical_mask]
    x4 = x_vars[[4]][[1]][logical_mask]
    yc = y[logical_mask]

    # Data-driven priors: Calculate means and standard deviations from the data
    # Different settings for meany based on the component
    if (component == 1) {
      inversesdy = rgamma(1,10,1000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 5, sdy/10)  # Use mean of 2 for component 1

    } else if (component == 2) {
      inversesdy = rgamma(1,500,2000000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 700, sdy/10)  # Use mean of 1 for other components

    } else{
      inversesdy = rgamma(1,10,10000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 10, sdy/10)

    }

    inversesdx1 = rgamma(1, (length(x1)-1)/2, sum((x1 - mean(x1))^2)/2)
    sdx1 = 1 / inversesdx1

    meanx1 = rnorm(1, mean(x1), sdx1)

    inversesdx2 = rgamma(1, (length(x2)-1)/2, sum((x2 - mean(x2))^2)/2)
    sdx2 = 1 / inversesdx2

    meanx2 = rnorm(1, mean(x2), sdx2)

    inversesdx3 = rgamma(1, (length(x3)-1)/2, sum((x3 - mean(x3))^2)/2)
    sdx3 = 1 / inversesdx3

    meanx3 = rnorm(1, mean(x3), sdx3)

    inversesdx4 = rgamma(1, (length(x4)-1)/2, sum((x4 - mean(x4))^2)/2)
    sdx4 = 1 / inversesdx4

    meanx4 = rnorm(1, mean(x4), sdx4)
  } else {
    # Fixed priors: Use user-provided fixed priors for the component
    if (component == 1) {
      priors <- user_fixed_priors$component1
    } else if (component == 2) {
      priors <- user_fixed_priors$component2
    } else if (component == 3) {
      priors <- user_fixed_priors$component3
    } else {
      stop("Invalid component specified!")
    }

    # Assign user-provided priors
    meany <- priors$meany
    meanx1 <- priors$meanx1
    meanx2 <- priors$meanx2
    meanx3 <- priors$meanx3
    meanx4 <- priors$meanx4
    sdy <- priors$sdy
    sdx1 <- priors$sdx1
    sdx2 <- priors$sdx2
    sdx3 <- priors$sdx3
    sdx4 <- priors$sdx4
  }

  # Calculate the log priors for each parameter based on the means and standard deviations
  a_prior = dnorm(a, meany, sdy, log = TRUE)
  b_prior = dnorm(b, meanx1, sdx1, log = TRUE)
  c_prior = dnorm(c, meanx2, sdx2, log = TRUE)
  d_prior = dnorm(d, meanx3, sdx3, log = TRUE)
  e_prior = dnorm(e, meanx4, sdx4, log = TRUE)

  # Return the sum of log priors for all parameters
  return(a_prior + b_prior + c_prior + d_prior + e_prior)
}
