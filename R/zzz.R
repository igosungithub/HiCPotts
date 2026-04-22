#' @useDynLib HiCPotts, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

#' @name pz_123
#' @title Internal: conditional component probabilities
#' @description Internal helper called by the C++ MCMC backend. Not for direct use.
#' @param z,sum_neighbours,y,pred_combined,chains,chain_gamma,x_vars,theta,size_chain,N,iter,dist Internal.
#' @return An N by N numeric matrix.
#' @keywords internal
#' @export
NULL

#' @name run_metropolis_MCMC_betas
#' @title Internal: Metropolis-Hastings MCMC backend
#' @description Internal C++ MCMC driver. Not for direct use.
#' @param N,gamma_prior,iterations,x_vars,y,use_data_priors,user_fixed_priors,dist,epsilon,distance_metric,size_start,theta_start Internal.
#' @return A list of MCMC results.
#' @keywords internal
#' @export
NULL