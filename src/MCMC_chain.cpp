#include <RcppArmadillo.h>
using namespace Rcpp;


// Helper function to update theta based on z_current and y
double update_theta(const NumericMatrix &z_current, const NumericMatrix &y) {
  int n1 = 0, n0 = 0;
  for (int i = 0; i < z_current.nrow(); i++) {
    for (int j = 0; j < z_current.ncol(); j++) {
      if(y(i ,j) == 0){
        if (z_current(i, j) == 1) ++n1;
          else ++n0;
      }    
    }
  }
  return R::rbeta(n1 + 1, n0 + 1);
}



// #' @title MCMC Chain Function for the HMRFHiC Model
// #'
// #' @description
// #' The function integrates multiple steps, including proposal distributions, acceptance/rejection steps based on posterior calculations,
// #' and Approximate Bayesian Computation (ABC) style updates for the \eqn{\gamma} parameter. Covariates and prior specifications can be either
// #' data-driven or user-defined. The output includes the chains of sampled parameters across the specified number of iterations.
// #'
// #'
// #' 
// #' @usage
// #' run_metropolis_MCMC_betas(N, gamma_prior, iterations, x_vars, y, theta_start = NULL,
// #'                           size_start = NULL, use_data_priors, user_fixed_priors = NULL,
// #'                           dist = "ZIP", epsilon = NULL,
// #'                           distance_metric = "manhattan", mc_cores
// #'                            )
// #'
// #' @param N Integer specifying the dimension of the lattice (\eqn{N \times N}).
// #'
// #' @param gamma_prior Numeric. The initial (prior) value for the interaction parameter \eqn{\gamma} 
// #'   in the Potts model.
// #'
// #' @param iterations Integer. The number of MCMC iterations to run.
// #'
// #' @param x_vars A list of covariates used as predictors in the model. This list should contain named elements 
// #'   corresponding to "distance", "GC", "TES", and "ACC", each containing a list of matrices of values for all \eqn{N \times N}.
// #'
// #' @param y An \eqn{N \times N} numeric matrix of observed interaction counts. Each element corresponds to interacting locus (i,j).
// #'
// #' @param use_data_priors Logical. If \code{TRUE}, data-driven priors will be used for each component. 
// #'   If \code{FALSE}, then \code{user_fixed_priors} must be provided.
// #'
// #' @param user_fixed_priors (Optional) A list of user-specified priors for the model components if 
// #'   \code{use_data_priors = FALSE}. Each components priors should be specified as a list of means and standard 
// #'   deviations for the parameters.
// #'
// #' @param dist A character string specifying the distribution family to use:
// #'   \itemize{
// #'     \item \code{"Poisson"}: Poisson distribution
// #'     \item \code{"NB"}: Negative Binomial distribution
// #'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
// #'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
// #'   }
// #'   The default is \code{"ZIP"}.
// #'
// #' @param epsilon (Optional) A numeric value specifying a tolerance (epsilon) for the ABC (Approximate Bayesian Computation) step 
// #'   used to update \eqn{\gamma}. If not provided, the function computes epsilon dynamically based on the data.
// #'
// #' @param distance_metric A character string specifying the distance metric used in the ABC step. 
// #'   Supported options:
// #'   \itemize{
// #'     \item \code{"manhattan"}: Use absolute differences in means.
// #'     \item \code{"euclidean"}: Use squared differences in means.
// #'   }
// #'
// #' @param size_start (Optional) Required if \code{dist} is "NB" or "ZINB". A numeric vector of length 3 specifying 
// #'   initial values for the size (overdispersion) parameter for each of the 3 components.
// #'
// #' @param theta_start (Optional) A numeric value providing an initial value for \eqn{\theta}, the zero-inflation 
// #'   parameter in ZIP/ZINB models. If not provided and \code{dist} is ZIP or ZINB, \eqn{\theta} is initialized to 0.5.
// #'
// #' @details
// #' The algorithm proceeds as follows:
// #' \enumerate{
// #'   \item \strong{Initialization:} 
// #'     - The component assignment matrix \eqn{z} is initialized based on \code{y}.
// #'     - Parameters (\eqn{\gamma}, \eqn{\theta}, and \eqn{\text{size}} if applicable) are initialized from given or default values.
// #'
// #'   \item \strong{MCMC Updates per Iteration:}
// #'     \itemize{
// #'       \item Update the latent component assignments \eqn{z} by proposing new states and using \code{pz_123} and \code{Neighbours_combined} 
// #'             to compute posterior probabilities.
// #'       \item Update the model parameters (betas) for each component using a Metropolis-Hastings step:
// #'         - Propose new parameter values from normal distributions.
// #'         - Compute posterior probabilities with \code{posterior_combined} and accept or reject the proposal.
// #'       \item If using NB or ZINB, similarly update the \code{size} (overdispersion) parameter via a Metropolis-Hastings step.
// #'       \item If using ZIP or ZINB, update the zero-inflation parameter \eqn{\theta} using \code{update_theta}.
// #'       \item Update \eqn{\gamma} using an ABC step:
// #'         - Propose a new \eqn{\gamma}.
// #'         - Generate synthetic data from the model with current parameters.
// #'         - Compare synthetic data and observed data using the chosen distance metric and \eqn{\epsilon} threshold.
// #'         - Accept or reject the new \eqn{\gamma} based on whether the distance is less than \eqn{\epsilon}.
// #'     }
// #'
// #'   \item \strong{Adaptive Tuning:}
// #'     - The proposal standard deviations for the parameter proposals are adaptively tuned to achieve a target acceptance rate.
// #'
// #' @return
// #' A \code{list} containing:
// #' \itemize{
// #'   \item \code{chains}: A list of three \eqn{(iterations+1) \times 5} matrices, each storing the parameter chains 
// #'         for one of the three components.
// #'   \item \code{gamma}: A numeric vector of length \eqn{iterations+1}, the chain of \eqn{\gamma} values.
// #'   \item \code{theta}: A numeric vector of length \eqn{iterations+1}, the chain of \eqn{\theta} values 
// #'         (or initialized and unchanged if not ZIP/ZINB).
// #'   \item \code{size}: A \eqn{3 \times (iterations+1)} matrix of size (overdispersion) parameters 
// #'         if \code{dist} is NB or ZINB; otherwise, this matrix might be unused.
// #' }
// #'
// #' @examples
// #' 
// #' N <- 10
// #' gamma_prior <- 0.5
// #' iterations <- 10
// #' x_vars <- list(
// #'   distance = list(matrix(runif(N*N, 0, 10), nrow=N)),
// #'   GC = list(matrix(runif(N*N, 0, 1), nrow=N)),
// #'   TES = list(matrix(runif(N*N, 0, 2), nrow=N)),
// #'   ACC = list(matrix(runif(N*N, 0, 5), nrow=N))
// #' )
// #' y <- matrix(rpois(N*N, lambda=5), nrow=N)
// #'
// #' # Using data-driven priors, ZIP distribution, no user_fixed_priors:
// #' results <- run_metropolis_MCMC_betas(
// #'   N = N,
// #'   gamma_prior = gamma_prior,
// #'   iterations = iterations,
// #'   x_vars = x_vars,
// #'   y = y,
// #'   use_data_priors = TRUE,
// #'   dist = "ZIP"
// #' )
// #'
// #' # Inspect gamma chain:
// #' plot(results$gamma, type='l')
// #' 
// #' @aliases run_metropolis_MCMC_betas
// #'
// #' @export
//
// [[Rcpp::export]]
List run_metropolis_MCMC_betas(int N, double gamma_prior, int iterations,
                               List x_vars, NumericMatrix y,
                               bool use_data_priors, Nullable<List> user_fixed_priors = R_NilValue,
                               std::string dist = "ZIP", Nullable<double> epsilon = R_NilValue,
                               std::string distance_metric = "manhattan",
                               Nullable<NumericVector> size_start = R_NilValue, Nullable<double> theta_start = R_NilValue) {
   // Map x1, x2, x3, and x4 to distance, GC, TES, and ACC
  List x1 = x_vars["distance"];
  List x2 = x_vars["GC"];
  List x3 = x_vars["TES"];
  List x4 = x_vars["ACC"];

  // Check if user_fixed_priors is needed
  if (!use_data_priors && user_fixed_priors.isNull()) {
    stop("Error: When use_data_priors is set to false, user_fixed_priors must be provided.");
  }

  // Initialize chains and other variables
  List chains = List::create(NumericMatrix(iterations + 1, 5),
                             NumericMatrix(iterations + 1, 5),
                             NumericMatrix(iterations + 1, 5));
  NumericVector chain_gamma(iterations + 1, NA_REAL);
  chain_gamma[0] = gamma_prior; // Initialize with gamma_prior
  NumericVector theta(iterations + 1, NA_REAL);
  // Initialize theta with theta_start if provided and if dist is ZIP or ZINB
  if (theta_start.isNotNull() && (dist == "ZIP" || dist == "ZINB")) {
    theta[0] = as<double>(theta_start);
  } else {
    theta[0] = 0.5; // Default initialization for theta if not using ZIP or ZINB
  }


  // Initialize acceptance counts and standard deviations for components 1, 2, and 3
  std::vector<int> acceptance_counts = {0, 0, 0};
  List sd_values = List::create(NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5),
                                NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5),
                                NumericVector::create(1.0, 0.5, 0.5, 0.5, 0.5));

  // Initialize size_chain as a matrix for each component
  NumericMatrix size_chain(3, iterations + 1);  // 3 components, each with a chain of size `iterations + 1`  
  if (dist == "NB" || dist == "ZINB") {
    if (size_start.isNotNull()) {
      NumericVector init_sizes = Rcpp::as<NumericVector>(size_start.get());
      if (init_sizes.size() != 3) stop("size_start should have 3 initial values for 3 components");
      for (int c = 0; c < 3; c++) size_chain(c, 0) = init_sizes[c];
    } else {
      stop("size_start must be provided for Negative Binomial or ZINB distribution.");
    }
  }

  // Adaptive tuning parameters for MCMC
  double target_acceptance_rate = 0.22;
  int adaptation_interval = 50;
  double adaptation_scaling = 1.2;


  //Epsilon for the ABC 
  bool epsilon_is_fixed = epsilon.isNotNull();
  double epsilon_value  = epsilon_is_fixed ? as<double>(epsilon) : 0.0;


  // Initialize z with zeros and setup function calls to R
  Function Neighbours_combined("Neighbours_combined");
  Function pred_combined("pred_combined");
  Function pz_123("pz_123");
  Function gamma_prior_value("gamma_prior_value");
  Function posterior_combined("posterior_combined");
  Function quantile("quantile");

  // Initialize first z values based on y
  NumericMatrix z_current (N,N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (y(i, j) == 0) {
        z_current(i, j) = 1;
      } else {
        z_current(i, j) = std::floor(R::runif(1, 4));
      }
    }
  }

  // ==========================================================================
  // MAIN MCMC LOOP
  // ==========================================================================

  for (int iter = 0; iter < iterations; iter++) {
   if ((iter + 1) % 50 == 0) {
     Rcout << "Iteration: " << iter + 1 << std::endl;
    }

    // Initialize Pros matrices and fill them based on y and z_current

    NumericMatrix Pros1(N, N), Pros2(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (y(i, j) == 0) {
          Pros1(i, j) = 1;
          Pros2(i, j) = 1;
        } else {
          do {
            Pros1(i, j) = std::floor(R::runif(1, 4));
          } while (Pros1(i, j) == z_current(i, j));

          do {
            Pros2(i, j) = std::floor(R::runif(1, 4));
          } while (Pros2(i, j) == z_current(i, j) || Pros2(i, j) == Pros1(i, j));
        }
      }
    }

    NumericMatrix sum1 = Neighbours_combined(z_current, N);
    NumericMatrix sum2 = Neighbours_combined(z_current, N, Pros1);
    NumericMatrix sum3 = Neighbours_combined(z_current, N, Pros2);

    NumericMatrix P_initials = pz_123(z_current, sum1, y, pred_combined, chains, chain_gamma, x_vars, theta[iter], size_chain, N, iter, dist);
    NumericMatrix P_proposed1 = pz_123(Pros1, sum2, y, pred_combined, chains, chain_gamma, x_vars, theta[iter], size_chain, N, iter, dist);
    NumericMatrix P_proposed2 = pz_123(Pros2, sum3, y, pred_combined, chains, chain_gamma, x_vars, theta[iter], size_chain, N, iter, dist);

    // Compute probabilities for sampling next z
    NumericMatrix z_next(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
       // double psum = exp(P_initials(i, j)) + exp(P_proposed1(i, j)) + exp(P_proposed2(i, j));
        if (y(i, j) == 0) {
          z_next(i, j) = 1;
          continue;
        } 
        // Use log-sum-exp for numerical stability
        double lp0 = P_initials(i, j);
        double lp1 = P_proposed1(i, j);
        double lp2 = P_proposed2(i, j);
        double m   = std::max(lp0, std::max(lp1, lp2));
        double e0 = std::exp(lp0 - m);
        double e1 = std::exp(lp1 - m);
        double e2 = std::exp(lp2 - m);
        double s  = e0 + e1 + e2;
        double p0 = e0 / s;
        double p1 = e1 / s;
        double u  = R::runif(0, 1);
         if      (u < p0)         z_next(i, j) = z_current(i, j);
        else if (u < p0 + p1)    z_next(i, j) = Pros1(i, j);
        else                     z_next(i, j) = Pros2(i, j);
      }
    }

    // MCMC proposal update for betas, per component ---------------

    for (int comp = 1; comp <= 3; comp++) {
      NumericMatrix chain_matrix = chains[comp - 1];
      NumericVector proposal(5);
      NumericVector sd_component = as<NumericVector>(sd_values[comp - 1]);

      for (int j = 0; j < 5; j++) {
        proposal[j] = R::rnorm(chain_matrix(iter, j), sd_component[j]);
      }

      double current_size = size_chain(comp - 1, iter);
      double posterior_current = as<double>(posterior_combined(pred_combined, chain_matrix(iter, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, current_size));
      double posterior_proposal = as<double>(posterior_combined(pred_combined, proposal, z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, current_size));

      // Call proposaldensity_combined to get log proposal densities for current and proposed states
       // SYMMETRIC RW proposal => proposal-density ratio cancels (FIX)
      double log_alpha = posterior_proposal - posterior_current;
      if (std::log(R::runif(0, 1)) < log_alpha) {
        for (int k = 0; k < 5; k++)
          chain_matrix(iter + 1, k) = proposal[k];
        acceptance_counts[comp - 1]++;
      } else {
        for (int k = 0; k < 5; k++)
          chain_matrix(iter + 1, k) = chain_matrix(iter, k);
      }

      // Update size for NB or ZINB if needed
      if (dist == "NB" || dist == "ZINB") {
        double size_proposal = R::rnorm(current_size, 0.5);  // Propose a new size for this component
        if (size_proposal > 0) {
          double posterior_current_size = as<double>(posterior_combined(pred_combined, chain_matrix(iter + 1, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, current_size));
          double posterior_proposal_size = as<double>(posterior_combined(pred_combined, chain_matrix(iter + 1, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, size_proposal));

            double log_alpha_s = posterior_proposal_size - posterior_current_size;
          if (std::log(R::runif(0, 1)) < log_alpha_s)
            size_chain(comp - 1, iter + 1) = size_proposal;
          else
            size_chain(comp - 1, iter + 1) = current_size;
        } else {
          size_chain(comp - 1, iter + 1) = current_size;
        }
      }
    }

    // Adaptive tuning of sd_values
    if (((iter+1) % adaptation_interval == 0) && ((iter + 1) <= iterations / 2)) {
      for (int comp = 1; comp <= 3; comp++) {
        double acceptance_rate = static_cast<double>(acceptance_counts[comp - 1]) / adaptation_interval;
        NumericVector sd_component = as<NumericVector>(sd_values[comp - 1]);
        if (acceptance_rate < target_acceptance_rate) {
          sd_component = sd_component / adaptation_scaling;
        } else if (acceptance_rate > target_acceptance_rate) {
          sd_component = sd_component * adaptation_scaling;
        }
        sd_values[comp - 1] = sd_component;
        acceptance_counts[comp - 1] = 0;
      }
    }

        // --- ABC for updating gamma ---
      // Step 1: Propose a new gamma
    double gamma_proposed = as<double>(gamma_prior_value());

    NumericMatrix gamma_matrix(N, N); // Initialize an empty N x N matrix
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        gamma_matrix(i, j) = gamma_proposed; // Fill all elements with gamma_current
      }
    }

    // Extract beta values (5 parameters per component) for the current iteration
    NumericMatrix beta_chain_1 = chains[0]; // First component
    NumericMatrix beta_chain_2 = chains[1]; // Second component
    NumericMatrix beta_chain_3 = chains[2]; // Third component
    
    // Get the beta parameters for the current iteration
    NumericVector beta_current_1 = beta_chain_1(iter+1, _); // 5 parameters
    NumericVector beta_current_2 = beta_chain_2(iter+1, _); // 5 parameters
    NumericVector beta_current_3 = beta_chain_3(iter+1, _); // 5 parameters



    // Step 2: Compute lambda using likelihood_gamma
    NumericMatrix neighbours = Neighbours_combined(z_next, N);

    // likelihood_gamma now returns the raw Potts log-potential term gamma * neighbours
    Function likelihood_gamma("likelihood_gamma");
    NumericMatrix base_lambda = likelihood_gamma(gamma_matrix, neighbours, N);

    // Covariate matrices
    NumericMatrix x11 = as<NumericMatrix>(x1[0]);
    NumericMatrix x22 = as<NumericMatrix>(x2[0]);
    NumericMatrix x33 = as<NumericMatrix>(x3[0]);
    NumericMatrix x44 = as<NumericMatrix>(x4[0]);

    // Step 3a: build lambda_matrix first
    NumericMatrix lambda_matrix(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        NumericVector b;
        int zz = (int) z_next(i, j);

        if (zz == 1) {
          b = beta_current_1;
        } else if (zz == 2) {
          b = beta_current_2;
        } else {
          b = beta_current_3;
        }

        double eta =
          base_lambda(i, j) +
          b[0] +
          b[1] * std::log1p(x11(i, j)) +
          b[2] * std::log1p(x22(i, j)) +
          b[3] * std::log1p(x33(i, j)) +
          b[4] * std::log1p(x44(i, j));

        double mu = std::exp(eta);
        if (!std::isfinite(mu) || mu <= 0.0) mu = 1e-6;
        lambda_matrix(i, j) = mu;
      }
    }

    // Step 3b: simulate synthetic data
    NumericMatrix synthetic_data(N, N);
    double next_theta = theta[iter];

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double mu = lambda_matrix(i, j);
        int zz = (int) z_next(i, j);
        double sz = size_chain(zz - 1, iter);
        if (!(sz > 0.0) || std::isnan(sz)) sz = 1.0;

        if (dist == "Poisson") {
          synthetic_data(i, j) = R::rpois(mu);

        } else if (dist == "ZIP") {
          if (zz == 1) {
            double u = R::runif(0, 1);
            synthetic_data(i, j) = (u < next_theta) ? 0.0 : R::rpois(mu);
          } else {
            synthetic_data(i, j) = R::rpois(mu);
          }

        } else if (dist == "NB") {
          double prob = sz / (mu + sz);
          if (!(prob > 0.0 && prob < 1.0)) prob = 0.5;
          synthetic_data(i, j) = R::rnbinom(sz, prob);

        } else if (dist == "ZINB") {
          double prob = sz / (mu + sz);
          if (!(prob > 0.0 && prob < 1.0)) prob = 0.5;

          if (zz == 1) {
            double u = R::runif(0, 1);
            synthetic_data(i, j) = (u < next_theta) ? 0.0 : R::rnbinom(sz, prob);
          } else {
            synthetic_data(i, j) = R::rnbinom(sz, prob);
          }

        } else {
          stop("Unsupported distribution specified.");
        }
      }
    }


    // Step 4: Compute distance metric (absolute difference of means)
    // Distance: MEAN ABSOLUTE difference (Manhattan) or root-mean-square
    // (Euclidean) between observed and synthetic data. Same scale as y.
    double dist_stat = 0.0;
    double ssq = 0.0, sabs = 0.0;
    int ncells = N * N;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double d = y(i, j) - synthetic_data(i, j);
        sabs += std::abs(d);
        ssq  += d * d;
      }
    }
    if (distance_metric == "euclidean") dist_stat = std::sqrt(ssq / ncells);
    else                                dist_stat = sabs / ncells;    // "manhattan"
 
    // Epsilon: data-driven default is 10% quantile of per-cell |y - syn|
    // (same scale as dist_stat). User-supplied epsilon overrides.
    if (!epsilon_is_fixed) {
      NumericVector diffs(ncells);
      int k = 0;
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
          diffs[k++] = std::abs(y(i, j) - synthetic_data(i, j));
      epsilon_value = as<double>(quantile(diffs, 0.1));
      if (!(epsilon_value > 0.0)) epsilon_value = 1e-6;
    }
 
    chain_gamma[iter + 1] = (dist_stat < epsilon_value) ? gamma_proposed
                                                        : chain_gamma[iter];
    
     // Update theta if using ZIP or ZINB distribution
    if (dist == "ZIP" || dist == "ZINB")
      theta[iter + 1] = update_theta(z_next, y);
      else
      theta[iter + 1] = theta[iter];
     z_current = z_next;

  }

  return List::create(Named("chains") = chains, Named("gamma") = chain_gamma, Named("theta") = theta, Named("size") = size_chain);
}





