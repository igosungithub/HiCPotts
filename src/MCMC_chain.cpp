#include <RcppArmadillo.h>
using namespace Rcpp;


// Helper function to update theta based on z_current and y
double update_theta(const NumericMatrix &z_current, const NumericMatrix &y) {
  int n1 = 0, n0 = 0;
  for (int i = 0; i < z_current.nrow(); i++) {
    for (int j = 0; j < z_current.ncol(); j++) {
      if (z_current(i, j) == 1) {
        n1++;
      } else {
        n0++;
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
// #' \dontrun{
// #' N <- 10
// #' gamma_prior <- 0.5
// #' iterations <- 1000
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
// #' }
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
                               Nullable<NumericMatrix> size_start = R_NilValue, Nullable<double> theta_start = R_NilValue) {
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
  List sd_values = List::create(NumericVector::create(70, 50, 50, 50, 10),
                                NumericVector::create(70, 50, 50, 50, 10),
                                NumericVector::create(70, 20, 70, 70, 70));

  // Initialize size_chain as a matrix for each component
  NumericMatrix size_chain(3, iterations + 1);  // 3 components, each with a chain of size `iterations + 1`  
  if (dist == "NB" || dist == "ZINB") {
    if (size_start.isNotNull()) {
      NumericVector init_sizes = Rcpp::as<NumericVector>(size_start.get());
      if (init_sizes.size() != 3) stop("size_start should have 3 initial values for 3 components");
      for (int comp = 1; comp <= 3; comp++) {
        size_chain(comp - 1, 0) = init_sizes[comp - 1];  // Initialize each components size_chain with size_start
      }
    } else {
      stop("Error: size_start must be provided for Negative Binomial or ZINB distribution.");
    }
  }

  // Adaptive tuning parameters for MCMC
  double target_acceptance_rate = 0.2;
  int adaptation_interval = 5;
  double adaptation_scaling = 1.5;


  //Epsilon for the ABC 
    double epsilon_value;
    if (epsilon.isNotNull()) {
        // Use the user-provided epsilon
        epsilon_value = as<double>(epsilon);
    }else {
        // Placeholder for dynamic epsilon calculation
        epsilon_value = 0.0; // Will be set dynamically in the ABC section based on quantile
    }


  // Initialize z with zeros and setup function calls to R
  List z;
  z.push_back(NumericMatrix(N, N));
  Function Neighbours_combined("Neighbours_combined");
  Function pred_combined("pred_combined");
  Function pz_123("pz_123");
  Function gamma_prior_value("gamma_prior_value");
  Function posterior_combined("posterior_combined");
  Function proposaldensity_combined("proposaldensity_combined");
  Function likelihood_gamma("likelihood_gamma");
  Function quantile("quantile");

  // Initialize first z values based on y
  NumericMatrix z_current = as<NumericMatrix>(z[0]);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (y(i, j) == 0) {
        z_current(i, j) = 1;
      } else {
        z_current(i, j) = floor(R::runif(1, 4));
      }
    }
  }
  z[0] = z_current;

  for (int iter = 0; iter < iterations; iter++) {
    Rcout << "Iteration: " << iter + 1 << std::endl;

    // Initialize Pros matrices and fill them based on y and z_current

    NumericMatrix Pros1(N, N), Pros2(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (y(i, j) == 0) {
          Pros1(i, j) = 1;
          Pros2(i, j) = 1;
        } else {
          do {
            Pros1(i, j) = floor(R::runif(1, 4));
          } while (Pros1(i, j) == z_current(i, j));

          do {
            Pros2(i, j) = floor(R::runif(1, 4));
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
        double psum = exp(P_initials(i, j)) + exp(P_proposed1(i, j)) + exp(P_proposed2(i, j));
        if (y(i, j) == 0) {
          z_next(i, j) = 1;
        } else {
          double prob1 = exp(P_initials(i, j)) / psum;
          double prob2 = exp(P_proposed1(i, j)) / psum;
          double random_val = R::runif(0, 1);
          z_next(i, j) = (random_val < prob1) ? 1 : ((random_val < prob1 + prob2) ? 2 : 3);
        }
      }
    }
    z.push_back(z_next);

    // MCMC proposal

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
      double proposal_density_current = as<double>(proposaldensity_combined(chain_matrix(iter, _), comp));
      double proposal_density_proposal = as<double>(proposaldensity_combined(proposal, comp));


      double probab = (posterior_proposal + proposal_density_current)- (posterior_current + proposal_density_proposal);
      if (log(R::runif(0, 1)) < probab) {
        for (int j = 0; j < 5; j++) {
          chain_matrix(iter + 1, j) = proposal[j];
        }
        acceptance_counts[comp - 1]++;
      } else {
        for (int j = 0; j < 5; j++) {
          chain_matrix(iter + 1, j) = chain_matrix(iter, j);
        }
      }

      // Update size for NB or ZINB if needed
      if (dist == "NB" || dist == "ZINB") {
        double size_proposal = R::rnorm(current_size, 3.0);  // Propose a new size for this component
        if (size_proposal > 0) {
          double posterior_current_size = as<double>(posterior_combined(pred_combined, chain_matrix(iter, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, current_size));
          double posterior_proposal_size = as<double>(posterior_combined(pred_combined, chain_matrix(iter, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors, dist, size_proposal));

          double size_acceptance_prob = posterior_proposal_size - posterior_current_size;
          if (log(R::runif(0, 1)) < size_acceptance_prob) {
            size_chain(comp - 1, iter + 1) = size_proposal;
          } else {
            size_chain(comp - 1, iter + 1) = current_size;
          }
        }
      }
    }

    // Adaptive tuning of sd_values
    if (iter % adaptation_interval == 0 && iter <= iterations / 2) {
      for (int comp = 1; comp <= 3; comp++) {
        double acceptance_rate = (double)acceptance_counts[comp - 1] / adaptation_interval;
        NumericVector sd_component = as<NumericVector>(sd_values[comp - 1]);
        if (acceptance_rate < target_acceptance_rate) {
          sd_component = sd_component * (1.0 / adaptation_scaling);
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
        gamma_matrix(i, j) = chain_gamma[iter]; // Fill all elements with gamma_current
      }
    }

    // Extract beta values (5 parameters per component) for the current iteration
    NumericMatrix beta_chain_1 = chains[0]; // First component
    NumericMatrix beta_chain_2 = chains[1]; // Second component
    NumericMatrix beta_chain_3 = chains[2]; // Third component
    
    // Get the beta parameters for the current iteration
    NumericVector beta_current_1 = beta_chain_1(iter, _); // 5 parameters
    NumericVector beta_current_2 = beta_chain_2(iter, _); // 5 parameters
    NumericVector beta_current_3 = beta_chain_3(iter, _); // 5 parameters



   // Step 2: Compute lambda using likelihood_gamma
    NumericMatrix neighbours = Neighbours_combined(z_current, N);
    NumericMatrix lambda_matrix = likelihood_gamma(gamma_matrix, neighbours, N);


// Debugging: Output lambda before beta adjustments
//Rcout << "Lambda[50,50] before beta adjustment: " << lambda_matrix(50, 50) << std::endl;

    // Incorporate beta parameters into lambda_matrix
    // Use x_vars["distance"], x_vars["GC"], x_vars["TES"], x_vars["ACC"] for component 1
    // Cast each element of the list to NumericMatrix
        NumericMatrix x11 = as<NumericMatrix>(x1[0]);
        NumericMatrix x22 = as<NumericMatrix>(x2[0]);
        NumericMatrix x33 = as<NumericMatrix>(x3[0]);
        NumericMatrix x44 = as<NumericMatrix>(x4[0]);
    
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
      //  double lambda = lambda_matrix(i, j); // Start with initial lambda from likelihood_gamma
         if (z_current(i, j) == 1) {
            // Example: Adjust lambda using all 5 beta parameters for each component
            lambda_matrix(i, j) += beta_current_1[0] +
                                   (beta_current_1[1] * std::log(1 + x11(i,j))) +
                                   (beta_current_1[2] * std::log(1 + x22(i,j)))+
                                   (beta_current_1[3] * std::log(1 + x33(i,j))) +
                                   (beta_current_1[4] * std::log(1 + x44(i,j)));
                                   
                                 
            // Additional terms from other components, if needed
          } else if (z_current(i, j) == 2) {
            lambda_matrix(i, j) += beta_current_2[0]  +
                                   (beta_current_2[1] * std::log(1 + x11(i,j))) +
                                   (beta_current_2[2] * std::log(1 + x22(i,j))) +
                                   (beta_current_2[3] * std::log(1 + x33(i,j))) +
                                   (beta_current_2[4] * std::log(1 + x44(i,j)));
    
          } else if (z_current(i, j) == 3) {
            lambda_matrix(i, j) += beta_current_3[0] +
                                   (beta_current_3[1] * std::log(1 + x11(i,j))) +
                                   (beta_current_3[2] * std::log(1 + x22(i,j))) +
                                   (beta_current_3[3] * std::log(1 + x33(i,j))) +
                                   (beta_current_3[4] * std::log(1 + x44(i,j)));
          }
           // Update the lambda_matrix with the adjusted value
       // lambda_matrix(i, j) = lambda;
      }
    }

// Debugging: Output lambda after beta adjustments
//Rcout << "Lambda[50,50] after beta adjustment: " << lambda_matrix(50, 50) << std::endl;


    // Step 3: Simulate synthetic data
    NumericMatrix synthetic_data(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
        
        lambda_matrix(i, j) = exp(lambda_matrix(i, j));
        
        if (lambda_matrix(i, j) < 0) {
    //Rcout << "Invalid lambda: " << lambda_matrix(i, j) << " at (" << i << ", " << j << ")" << std::endl;
    lambda_matrix(i, j) = 1e-6; // Replace with a fallback
        }
        int component = z_current(i, j) - 1; // Assuming values are 1, 2, 3
        double size = size_chain(component, iter);

        // Validate size
        if (size <= 0 || std::isnan(size)) {
          //  Rcout << "Invalid size at (" << i << ", " << j << "): " << size << std::endl;
            size = 1e-6; // Fallback value
        }


            if (dist == "Poisson") {
                synthetic_data(i, j) = R::rpois(lambda_matrix(i, j));
            } else if (dist == "ZIP") {
                double lambda = lambda_matrix(i, j);
                double u = R::runif(0, 1);
                synthetic_data(i, j) = (u < theta[iter]) ? 0 : R::rpois(lambda);
            } else if (dist == "NB") {
                //double lambda = lambda_matrix(i, j);
                //int component = z_current(i, j) - 1; // Assuming z_current values are 1, 2, or 3
      
                double prob = size / (lambda_matrix(i, j) + size);
            if (prob <= 0 || prob >= 1 || std::isnan(prob)) {
               // Rcout << "Invalid prob at (" << i << ", " << j << "): " << prob << std::endl;
                prob = 0.5; // Fallback value
            }
                synthetic_data(i, j) = R::rnbinom(size, prob);
            } else if (dist == "ZINB") {
                //double lambda = lambda_matrix(i, j);
                double u = R::runif(0, 1); 
                
                double prob = size / (lambda_matrix(i, j) + size);
            if (prob <= 0 || prob >= 1 || std::isnan(prob)) {
                //Rcout << "Invalid prob at (" << i << ", " << j << "): " << prob << std::endl;
                prob = 0.5; // Fallback value
            }
                synthetic_data(i, j) = (u < theta[iter]) ? 0 : R::rnbinom(size, prob); 
            } else {
                stop("Unsupported distribution specified.");
            }
            if (std::isnan(synthetic_data(i, j))) {
               //Rcout << "NaN in synthetic_data at (" << i << ", " << j << ")" << std::endl;
            }
        }
    }
  



    // Step 4: Compute distance metric (absolute difference of means)
      double mean_y = mean(y);
      double mean_synthetic_data = mean(synthetic_data);
      double distance = (distance_metric == "manhattan") ? 
                  std::abs(mean_y - mean_synthetic_data) : 
                  std::sqrt(std::pow(mean_y - mean_synthetic_data, 2));


    // Step 6: Calculate epsilon dynamically if not provided
    if (epsilon.isNull()) {
    NumericVector raw_differences;
       // Vectors to store non-zero values
      //NumericMatrix raw_differences(N, N);
        // Compute the raw differences (no absolute value) between y and synthetic_data
        for (int i = 0; i < N; i++) {
         for (int j = 0; j < N; j++) {
            if (!(y(i, j) == 0 && synthetic_data(i, j) == 0)) {
                raw_differences.push_back(std::exp(std::abs(y(i, j) - synthetic_data(i, j))));
            }
            if (R_IsNA(y(i, j)) || R_IsNA(synthetic_data(i, j))) {
    //Rcout << "NA detected at position (" << i << ", " << j << ")" << std::endl;
}

          }
        }

           // Ensure raw differences are non-empty
         if (raw_differences.size() == 0) {
            stop("Raw differences array is empty. Cannot compute epsilon.");
          }
        // Flatten raw_differences into a NumericVector for quantile calculation
        NumericVector flat_differences = as<NumericVector>(raw_differences);
        // Compute the 1% quantile
        // Use R quantile function
            epsilon_value = as<double>(quantile(flat_differences, 0.2)); // 1% quantile

          //Rcout << "Iteration " << iter + 1 << ": raw_differences[0] = " << raw_differences(0) << std::endl;

          //  Rcout << "Iteration " << iter + 1 
          //  << ": Distance = " << distance 
          //  << ", Epsilon = " << epsilon_value 
          //  << ", Gamma Proposed = " << gamma_proposed 
          //  << ", Gamma Accepted = " << (distance < epsilon_value) 
          //  << std::endl;
     }



    // Step 7: Accept or reject gamma based on distance
    
    if (distance < epsilon_value) {
      chain_gamma[iter + 1] = gamma_proposed;
    } else {
      chain_gamma[iter + 1] = chain_gamma[iter];
    }
    
     // Update theta if using ZIP or ZINB distribution
    if (dist == "ZIP" || dist == "ZINB") {
      theta[iter + 1] = update_theta(z_next, y);
    }

  }

  return List::create(Named("chains") = chains, Named("gamma") = chain_gamma, Named("theta") = theta, Named("size") = size_chain);
}





