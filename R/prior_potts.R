#' @title Prior Value for the Potts Model Interaction Parameter
#'
#' @description
#' This function generates a prior value for the interaction parameter in a Potts model from a Beta distribution.
#' The Potts model is a spatial statistical model where the interaction parameter influences the tendency of
#' neighboring sites on a lattice to take on similar states. By drawing the interaction parameter from a Beta
#' distribution, we impose a prior belief on the range and likely values of this parameter.
#'
#' @usage
#' gamma_prior_value()
#'
#' @details
#' The function samples a single random value from a \code{Beta(10, 5)} distribution. This distribution places
#' more mass towards the higher end of the \[0, 1\] interval (since Beta(10, 5) is skewed towards 1), indicating
#' a prior belief that the interaction parameter is likely to encourage some degree of spatial clustering.
#'
#'
#' @return
#' A numeric value between \[0,1\] representing the prior draw for the interaction parameter. This is a random draw
#' from the specified Beta distribution.
#'
#' @examples
#' #\donttest{
#' # Generate a single prior value for the Potts model interaction parameter
#' prior_val <- gamma_prior_value()
#' prior_val
#' #
#' #}
#'
#' @export
#
gamma_prior_value <- function() {
  rbeta(1, 10, 5)
}
