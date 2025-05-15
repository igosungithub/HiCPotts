#include <RcppArmadillo.h>
using namespace Rcpp;


// Helper function to create shifted matrices
NumericMatrix shift_matrix(const NumericMatrix& data, std::string shift_direction, int N) {
  NumericMatrix shifted(N, N);

  if (shift_direction == "down") {
    for (int i = 1; i < N; i++) {
      for (int j = 0; j < N; j++) {
        shifted(i, j) = data(i - 1, j);
      }
    }
  } else if (shift_direction == "up") {
    for (int i = 0; i < N - 1; i++) {
      for (int j = 0; j < N; j++) {
        shifted(i, j) = data(i + 1, j);
      }
    }
  } else if (shift_direction == "right") {
    for (int i = 0; i < N; i++) {
      for (int j = 1; j < N; j++) {
        shifted(i, j) = data(i, j - 1);
      }
    }
  } else if (shift_direction == "left") {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N - 1; j++) {
        shifted(i, j) = data(i, j + 1);
      }
    }
  }

  return shifted;
}


// Helper function to compute neighbours
NumericMatrix compute_neighbours(const NumericMatrix& data1, const NumericMatrix& compare_matrix, int N) {
  NumericMatrix result_matrix(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (data1(i, j) == compare_matrix(i, j)) {
        result_matrix(i, j) = 1;
      }
    }
  }
  return result_matrix;
}


//' @name Neighbours_combined
//' @title Neighbours Function for the Potts Model
//'
//' @description
//' Computes the number of agreeing neighbors for each site in a Potts model lattice.
//' The Potts model is a generalization of the Ising model, where each site on a lattice can 
//' take on multiple states. By counting how many of a site's four neighbors share the same 
//' state, this function provides a measure of clustering within the lattice.
//'
//' @usage
//' Neighbours_combined(potts_data, N, proposed_value = NULL)
//'
//' @param potts_data A numeric \eqn{N \times N} matrix representing the current configuration 
//'   of the Potts model. Each element corresponds to the state of a particular site.
//'
//' @param N An integer specifying the dimension of the \code{potts_data} lattice (i.e., the 
//'   lattice has \code{N} rows and \code{N} columns).
//'
//' @param proposed_value (Optional) A numeric \eqn{N \times N} matrix representing a proposed 
//'   configuration of the Potts model. If provided, the function computes neighbor relationships 
//'   relative to this proposed configuration; otherwise, it uses a shifted version of \code{potts_data} 
//'   to determine neighbors.
//'
//' @details
//' The function checks each site and compares it with its four directional neighbors 
//' (up, down, left, right). For each neighbor that has the same state, a value of 1 is assigned. 
//' Summing these values for each site results in a measure of how "similar" the immediate neighborhood 
//' is to that site.
//'
//' Internally, the function:
//' \enumerate{
//'   \item Uses \code{shift_matrix} to create versions of the lattice shifted in each of the four directions.
//'   \item Uses \code{compute_neighbours} to determine if neighboring positions match in state.
//'   \item Sums the contributions from all four directions to form the \code{Neighbours_total} matrix.
//' }
//'
//' By providing an optional \code{proposed_value} matrix, this function can be integrated into an 
//' MCMC algorithm where new configurations are proposed and evaluated 
//' based on their local coherence.
//'
//' @return
//' A numeric \eqn{N \times N} matrix, \code{Neighbours_total}, where each element is the count of 
//' how many of the four neighbors match the state of that site.
//'
//' @examples
//' \dontrun{
//' # Suppose we have a 5x5 Potts model lattice:
//' potts_data <- matrix(sample(1:3, 25, replace = TRUE), ncol = 5)
//' N <- 5
//'
//' # Compute the neighbors without a proposed configuration:
//' neigh <- Neighbours_combined(potts_data, N)
//' neigh
//'
//' # Suppose we propose a new configuration:
//' proposed <- matrix(sample(1:3, 25, replace = TRUE), ncol = 5)
//'
//' # Evaluate neighbors for the proposed configuration:
//' neigh_proposed <- Neighbours_combined(potts_data, N, proposed_value = proposed)
//' neigh_proposed
//' }
//' @keywords internal
//'
//' @export
//
// [[Rcpp::export]]
NumericMatrix Neighbours_combined(NumericMatrix potts_data, int N, Nullable<NumericMatrix> proposed_value = R_NilValue) {
  bool is_proposed = proposed_value.isNotNull();
  NumericMatrix compare_matrix1(N, N);

  // Determine the shifted matrices and comparison matrix
  NumericMatrix mydata1 = shift_matrix(potts_data, is_proposed ? "up" : "down", N);
  if (is_proposed) {
    compare_matrix1 = as<NumericMatrix>(proposed_value);
  } else {
    compare_matrix1 = mydata1;
  }

  // Compute neighbour relationships
  NumericMatrix neighbour1 = compute_neighbours(mydata1, compare_matrix1, N);
  NumericMatrix neighbour2 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "down" : "up", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);
  NumericMatrix neighbour3 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "left" : "right", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);
  NumericMatrix neighbour4 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "right" : "left", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);

  // Calculating the total neighbours
  NumericMatrix Neighbours_total(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Neighbours_total(i, j) = neighbour1(i, j) + neighbour2(i, j) + neighbour3(i, j) + neighbour4(i, j);
    }
  }

  return as<NumericMatrix>(Neighbours_total);

}
