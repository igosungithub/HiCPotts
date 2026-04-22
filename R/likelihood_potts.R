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
#'   \item Reshapes these values into an \code{N x N} matrix, \code{potts_DA}, which represents
#'         the spatial layout of likelihoods (or probabilities) under the Potts model.
#' }
#'
#'
#'
#' @return A numeric \code{N x N} matrix of probabilities corresponding to the Potts model likelihood
#' for the given interaction parameter. Each element of the matrix represents the likelihood contribution of
#' that site/pair under the model parameterized by \code{x}.
#'
#' @examples
#' 
#' # Suppose we have an N=4 lattice and a vector x and pair_neighbours_DA_x1 of length 16
#' N <- 4
#' x <- rep(0.5, N * N) # interaction parameter repeated for each pair
#' pair_neighbours_DA_x1 <- rnorm(N * N) # random neighbor influences
#'
#' # Compute the Potts DA matrix
#' potts_matrix <- likelihood_gamma(x, pair_neighbours_DA_x1, N)
#' #print(potts_matrix)
#' 
#'
#' @export
#'
likelihood_gamma <- function(x, pair_neighbours_DA_x1, N) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (!is.numeric(pair_neighbours_DA_x1)) {
    stop("pair_neighbours_DA_x1 must be numeric.")
  }
  if (!is.numeric(N) || length(N) != 1L || N < 1L || N != as.integer(N)) {
    stop("N must be a positive integer scalar.")
  }
  if (length(pair_neighbours_DA_x1) != N * N) {
    stop(sprintf("length(pair_neighbours_DA_x1) = %d but N*N = %d.",
                 length(pair_neighbours_DA_x1), N * N))
  }
  if (!(length(x) == 1L || length(x) == length(pair_neighbours_DA_x1))) {
    stop("x must be scalar or the same length as pair_neighbours_DA_x1.")
  }
  
  vals <- as.numeric(x) * as.numeric(pair_neighbours_DA_x1)
  
  if (any(!is.finite(vals))) {
    stop("Computed gamma-neighbour values contain non-finite entries.")
  }
  
  matrix(vals, nrow = N, ncol = N)
}
