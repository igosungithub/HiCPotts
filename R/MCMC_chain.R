#' @title Run Metropolis-Hastings MCMC for HMRFHiC Model
#'
#' @description
#' Runs a Markov Chain Monte Carlo (MCMC) algorithm using Metropolis-Hastings to sample parameters for the Hidden Markov Random Field (HMRF) model in HMRFHiC.
#' Integrates proposal distributions, acceptance/rejection steps based on posterior calculations, and Approximate Bayesian Computation (ABC) updates for the \eqn{\gamma} parameter.
#' Supports data-driven or user-defined priors for covariates.
#'
#'   
#'
#' @param N Integer specifying the dimension of the lattice (\eqn{N \times N}).
#' @param gamma_prior Numeric. Initial (prior) value for the interaction parameter \eqn{\gamma} in the Potts model. Set at 0.3 as default.
#' @param iterations Integer. Number of MCMC iterations to run.
#' @param x_vars A list of covariates used as predictors. Must contain named elements "distance", "GC", "TES", and "ACC", each a list of \eqn{N \times N} matrices.
#' @param y An \eqn{N \times N} numeric matrix of observed interaction counts, with elements corresponding to locus pairs (i,j).
#' @param theta_start (Optional) Numeric initial value for \eqn{\theta}, the zero-inflation parameter in ZIP/ZINB models. Defaults to 0.5 if \code{NULL} and \code{dist} is ZIP/ZINB.
#' @param size_start (Optional) Numeric vector of length 3 specifying initial size (overdispersion) parameters for each component if \code{dist} is "NB" or "ZINB".
#' @param use_data_priors Logical. If \code{TRUE}, uses data-driven priors for each component. If \code{FALSE}, \code{user_fixed_priors} must be provided.
#' @param user_fixed_priors (Optional) A list of user-specified priors for model components if \code{use_data_priors = FALSE}. Each component’s priors should include means and standard deviations.
#' @param dist Character string specifying the distribution family: \code{"Poisson"}, \code{"NB"} (Negative Binomial), \code{"ZIP"} (Zero-Inflated Poisson, default), or \code{"ZINB"} (Zero-Inflated Negative Binomial).
#' @param epsilon (Optional) Numeric tolerance for the ABC update of \eqn{\gamma}. If \code{NULL}, computed dynamically from data.
#' @param distance_metric Character string specifying the distance metric for ABC: \code{"manhattan"} (default) or \code{"euclidean"}.
#' 
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item \strong{Initialization:}
#'     Initializes the component assignment matrix \eqn{z} based on \code{y} and sets starting values for parameters (\eqn{\gamma}, \eqn{\theta}, and \eqn{\text{size}} if applicable).
#'   \item \strong{MCMC Updates per Iteration:}
#'     \itemize{
#'       \item Updates latent assignments \eqn{z} by proposing new states and computing posterior probabilities using internal functions.
#'       \item Updates model parameters (betas) for each component via Metropolis-Hastings:
#'         \itemize{
#'           \item Proposes new values from normal distributions.
#'           \item Computes posterior probabilities and accepts/rejects proposals.
#'         }
#'       \item Updates overdispersion \code{size} or zero-inflation \eqn{\theta} parameters if applicable.
#'       \item Updates \eqn{\gamma} via ABC:
#'         \itemize{
#'           \item Proposes a new \eqn{\gamma} and simulates synthetic data.
#'           \item Compares to observed data using the specified \code{distance_metric} and \eqn{\epsilon}.
#'           \item Accepts/rejects based on whether the distance is below \eqn{\epsilon}.
#'         }
#'     }
#'   \item \strong{Adaptive Tuning:}
#'     Tunes proposal standard deviations to achieve a target acceptance rate.
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{chains}: List of three \eqn{(iterations+1) \times 5} matrices, each storing parameter chains for one of the three components.
#'   \item \code{gamma}: Numeric vector of length \eqn{iterations+1}, the chain of \eqn{\gamma} values.
#'   \item \code{theta}: Numeric vector of length \eqn{iterations+1}, the chain of \eqn{\theta} values (unchanged if not ZIP/ZINB).
#'   \item \code{size}: \eqn{3 \times (iterations+1)} matrix of size (overdispersion) parameters if \code{dist} is NB or ZINB; otherwise unused.
#' }
#'
#' @examples
#' \donttest{
#' N <- 10
#' gamma_prior <- 0.5
#' iterations <- 100
#' x_vars <- list(
#'   distance = list(matrix(runif(N*N, 0, 10), nrow=N)),
#'   GC = list(matrix(runif(N*N, 0, 1), nrow=N)),
#'   TES = list(matrix(runif(N*N, 0, 2), nrow=N)),
#'   ACC = list(matrix(runif(N*N, 0, 5), nrow=N))
#' )
#' y <- matrix(rpois(N*N, lambda=5), nrow=N)
#' results <- run_metropolis_MCMC_betas(
#'   N = N,
#'   gamma_prior = gamma_prior,
#'   iterations = iterations,
#'   x_vars = x_vars,
#'   y = y,
#'   use_data_priors = TRUE,
#'   dist = "ZIP"
#' )
#' plot(results$gamma, type="l")
#' }
#'
#' @useDynLib HMRFHiC, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
run_metropolis_MCMC_betas <- function(
    N, gamma_prior, iterations, x_vars, y, theta_start = NULL,
    size_start = NULL, use_data_priors, user_fixed_priors = NULL,
    dist = "ZIP", epsilon = NULL, distance_metric = "manhattan"
) {
  # Handle mc_cores in R (e.g., set up parallelization if needed)
  if (mc_cores > 1) {
    warning("Parallelization with mc_cores > 1 is not yet implemented; using 1 core")
    mc_cores <- 1
  }
  
  # Call C++ function, excluding mc_cores
  .Call('_HMRFHiC_run_metropolis_MCMC_betas', PACKAGE = 'HMRFHiC',
        N, gamma_prior, iterations, x_vars, y, use_data_priors,
        user_fixed_priors, dist, epsilon, distance_metric, size_start, theta_start)
}