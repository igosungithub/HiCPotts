#' @title Posterior Function for the Potts Model
#'
#' @description
#' Computes log-posterior probabilities for each site in an \eqn{N \times N} lattice under the Potts model,
#' combined with an interaction count model (e.g., genomic interactions). Incorporates mixture components,
#' neighbor relationships, and distributional assumptions (Poisson, Negative Binomial, Zero-Inflated Poisson,
#' Zero-Inflated Negative Binomial).
#'
#' @usage
#' pz_123(z, sum_neighbours, y, pred_combined, chains, chain_gamma,
#'        x_vars, theta, size_chain, N, iter, dist)
#'
#' @param z A numeric \eqn{N \times N} matrix representing the component assignments for each site on the lattice.
#'   Each element is a component label (e.g., 1, 2, or 3), indicating which mixture component that interaction currently belongs to.
#'
#' @param sum_neighbours An \eqn{N \times N} numeric matrix representing the number of neighboring sites
#'   (up, down, left, right) that agree in state with each site. Typically computed by the function
#'   \code{Neighbours_combined()}.
#'
#' @param y An \eqn{N \times N} numeric matrix of observed interaction counts (e.g., genomic interaction data).
#'   Each element corresponds to the count observed at that site.
#'
#' @param pred_combined An R function (passed from R to C++) that computes predicted values (linear predictors)
#'   for the specified component. This is often a reference to the \code{pred_combined()} function in R.
#'
#' @param chains A \code{list} of parameter vectors (MCMC chains) for the model. Each element corresponds to
#'   a component and contains the parameters associated with that component.
#'
#' @param chain_gamma A numeric vector (using ABC algorithm) for the interaction parameter \eqn{\gamma} that governs
#'   neighbor effects in the Potts model. The \code{iter}-th element of this vector is used in the calculation.
#'
#' @param x_vars A list of covariates (sources of biases), each a nested structure as used by
#'   \code{pred_combined()}. These covariates are used to compute component-specific predicted values (\eqn{\lambda}).
#'
#' @param theta A numeric scalar representing the zero-inflation parameter. Used when the distribution is
#'   Zero-Inflated Poisson (ZIP) or Zero-Inflated Negative Binomial (ZINB), typically for the first component.
#'
#' @param size_chain An \eqn{3 \times M} matrix representing the MCMC draws of the size (overdispersion) parameter
#'   for Negative Binomial or Zero-Inflated Negative Binomial distributions. One row per component, and columns
#'   correspond to iterations. The value at \code{size_chain[comp - 1, iter]} is used for the specified iteration.
#'
#' @param N An integer specifying the dimension of the lattice (\eqn{N \times N}).
#'
#' @param iter An integer indicating which iteration of the MCMC chains (parameters, \eqn{\gamma}, size) to use
#'   in the calculation.
#'
#' @param dist A character string specifying the distribution:
#'   \itemize{
#'     \item \code{"Poisson"}: Poisson distribution
#'     \item \code{"NB"}: Negative Binomial distribution
#'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
#'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
#'   }
#'
#' @details
#' This function interfaces with C++ code to calculate log-posterior probabilities for an \eqn{N \times N} lattice
#' under the Potts model. The model accounts for:
#' \enumerate{
#'   \item Component assignments (\code{z}) for each site.
#'   \item Observed interaction counts (\code{y}) from genomic data.
#'   \item Neighbor effects (\code{sum_neighbours}) weighted by \eqn{\gamma} (\code{chain_gamma}).
#'   \item Covariates (\code{x_vars}) used to compute predicted values (\eqn{\lambda}) via \code{pred_combined}.
#'   \item Distributional assumptions (\code{dist}: Poisson, NB, ZIP, or ZINB), with zero-inflation (\code{theta})
#'     for ZIP/ZINB (first component only) and overdispersion (\code{size_chain}) for NB/ZINB.
#' }
#' The \code{params} list bundles all necessary parameters, which are passed to the C++ function
#' \code{_HMRFHiC_potts_posterior}. Probabilities are computed for each site, adjusted by neighbor effects
#' using \eqn{\exp(\gamma \times \text{sum_neighbours})}, and returned as log-posteriors.
#'
#' @return
#' A numeric \eqn{N \times N} matrix of log-posterior probabilities for each site for the specified component.
#'
#' @examples
#' #\donttest{
#'  z <- matrix(sample(1:3, 25, replace = TRUE), nrow = 5)
#'  sum_neighbours <- matrix(runif(25, 0, 5), nrow = 5)
#'  y <- matrix(rpois(25, lambda = 5), nrow = 5)
#'  pred_combined <- function(params, z, x_vars, component, N) {
#'    return(runif(sum(z == component)))
#'  }
#'  chains <- list(
#'    matrix(rnorm(25), nrow = 5, ncol = 5),
#'    matrix(rnorm(25), nrow = 5, ncol = 5),
#'    matrix(rnorm(25), nrow = 5, ncol = 5)
#'  )
#'  chain_gamma <- rnorm(5)
#'  x_vars <- list(
#'    distance = list(matrix(runif(25, 0, 1), nrow = 5)),
#'    GC = list(matrix(runif(25, 0, 1), nrow = 5)),
#'    TES = list(matrix(runif(25, 0, 1), nrow = 5)),
#'    ACC = list(matrix(runif(25, 0, 1), nrow = 5))
#'  )
#'  theta <- 0.5
#'  size_chain <- matrix(runif(25, 0.5, 1.5), nrow = 5, ncol = 5)
#'  N <- nrow(z)
#'  iter <- 1
#'  dist <- "ZIP"
#'  
#'  # Run the function
#'  pz_vals <- pz_123(
#'    z = z,
#'    sum_neighbours = sum_neighbours,
#'    y = y,
#'    pred_combined = pred_combined,
#'    chains = chains,
#'    chain_gamma = chain_gamma,
#'    x_vars = x_vars,
#'    theta = theta,
#'    size_chain = size_chain,
#'    N = N,
#'    iter = iter,
#'    dist = dist
#'  )
#' 
#' print(pz_vals)
#' #}
#'
#' @useDynLib HMRFHiC, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
pz_123 <- function(z, sum_neighbours, y, pred_combined, chains, chain_gamma, x_vars, theta, size_chain, N, iter, dist) {
  .Call(
    "_HMRFHiC_pz_123",
    PACKAGE = "HMRFHiC",
    z, sum_neighbours, y, pred_combined, chains, chain_gamma, x_vars, theta, size_chain, N, iter, dist
  )
}
