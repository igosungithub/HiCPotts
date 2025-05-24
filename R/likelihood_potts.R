#' @title Likelihood Function for the Interaction Parameter in a Potts Model
#'
#' @description
#' This function computes the matrix of probabilities (or likelihood values) associated with the interaction parameter
#' in the Potts model. The Potts model is a generalization of the Ising model and is commonly used in spatial statistics,
#' image analysis, and other fields where data can be represented on a lattice. The interaction parameter controls
#' how neighboring sites on the lattice influence each other.
#'
#' @usage
#' likelihood_gamma(x, pair_neighbours_DA_x1, N)
#'
#' @param x A numeric vector representing the interaction parameter(s) for the Potts model.
#'   This could be a single parameter repeated for each pair of neighbors.
#'
#' @param pair_neighbours_DA_x1 A numeric vector of the same length as \code{x}, representing the pairwise interaction
#'   structure between neighbors in the Potts model. Each element indicates whether a pair of sites is assigned
#'   as neighbors.
#'
#' @param N An integer specifying the dimension of the Potts lattice. The returned matrix will be \code{N} by \code{N},
#'   representing the spatial layout of the Potts model.
#'
#' @details
#' The Potts model describes configurations of states on a lattice, where the interaction parameter \code{x}
#' influences how often neighboring sites are in the same state. The likelihood or probability of a
#' configuration is determined by exponentials of the interaction terms between neighboring sites.
#'
#' This function:
#' \enumerate{
#'   \item Multiplies the interaction parameter vector \code{x} by the vector \code{pair_neighbours_DA_x1}, which encodes
#'         neighbor relationships or their contributions.
#'   \item Applies a stability adjustment by subtracting the maximum value (\code{max_val}) from each term
#'         \eqn{x * pair\_neighbours\_DA\_x1} to avoid numerical overflow.
#'   \item Exponentiates these adjusted values and normalizes them to obtain a set of values that sum to 1,
#'         interpreting them as probabilities or likelihood contributions.
#'   \item Reshapes these normalized values into an \code{N x N} matrix, \code{potts_DA}, which represents
#'         the spatial layout of likelihoods (or probabilities) under the Potts model.
#' }
#'
#'
#'
#' @return A numeric \code{N x N} matrix of normalized probabilities corresponding to the Potts model likelihood
#' for the given interaction parameter. Each element of the matrix represents the likelihood contribution of
#' that site/pair under the model parameterized by \code{x}.
#'
#' @examples
#' #\donttest{
#' # Suppose we have an N=4 lattice and a vector x and pair_neighbours_DA_x1 of length 16
#' N <- 4
#' x <- rep(0.5, N * N) # interaction parameter repeated for each pair
#' pair_neighbours_DA_x1 <- rnorm(N * N) # random neighbor influences
#'
#' # Compute the Potts DA matrix
#' potts_matrix <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
#' print(potts_matrix)
#' #}
#'
#' @export
#'
likelihood_gamma <- function(x, pair_neighbours_DA_x1, N) {
  # Validate input dimensions
  if (length(x) != length(pair_neighbours_DA_x1)) {
    stop("x and pair_neighbours_DA_x1 must have the same length.")
  }

  # Step 1: Calculate max value for numerical stability
  max_val <- max(x * pair_neighbours_DA_x1)

  # Step 2: Calculate adjusted exponents
  exponent_diff <- x * pair_neighbours_DA_x1 - max_val
  exp_values <- exp(exponent_diff)

  # Step 3: Sum the adjusted exponentiated values
  sum_exp_values <- sum(exp_values)

  # Step 4: Calculate normalized probabilities
  a <- exp_values / sum_exp_values

  # Step 5: Create Potts DA matrix using the normalized probabilities
  potts_DA <- matrix(a, nrow = N, ncol = N)

  return(potts_DA)
}
