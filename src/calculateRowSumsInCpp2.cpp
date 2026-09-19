#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Function to calculate row sums for subsets of columns and count occurrences of 1 and 2
// for each iteration separately
// [[Rcpp::export]]
IntegerMatrix countOccurrencesInCpp(NumericMatrix mat, IntegerMatrix cols) {
  int nrow = mat.nrow();
  int ncalc = cols.nrow();
  IntegerMatrix counts(ncalc, 2); // Result matrix for counts of occurrences of 1 and 2
  
  // Iterate over calculations
  for (int k = 0; k < ncalc; ++k) {
    // Get column indices for current calculation
    IntegerVector cur_cols = cols(k, _);
    int ncol = cur_cols.size();
    
    // Initialize counts for current iteration
    int count1 = 0;
    int count2 = 0;
    
    // Iterate over rows
    for (int i = 0; i < nrow; ++i) {
      double sum = 0;
      // Iterate over selected columns
      for (int j = 0; j < ncol; ++j) {
        sum += mat(i, cur_cols[j] - 1); // Adjust for 0-based indexing
      }
      // Count occurrences of 1 and 2
      if (sum == 1) {
        count1++;
      } else if (sum == 2) {
        count2++;
      }
    }
    
    // Store counts for current iteration in result matrix
    counts(k, 0) = count1;
    counts(k, 1) = count2;
  }
  
  return counts;
}



// Accepts a Matrix::dgCMatrix
// [[Rcpp::export]]
IntegerMatrix countOccurrencesSparseCpp(
    S4 mat,
    IntegerMatrix cols
) {
  IntegerVector dim = mat.slot("Dim");
  IntegerVector p   = mat.slot("p");
  IntegerVector row = mat.slot("i");
  NumericVector x   = mat.slot("x");

  const int nrow  = dim[0];
  const int ncol  = dim[1];
  const int ncalc = cols.nrow();
  const int nsel  = cols.ncol();

  IntegerMatrix counts(ncalc, 2);

  // Reused between calculations
  std::vector<double> sums(nrow, 0.0);
  std::vector<unsigned char> seen(nrow, 0);
  std::vector<int> touched;
  touched.reserve(nrow);

  for (int k = 0; k < ncalc; ++k) {
    touched.clear();

    for (int j = 0; j < nsel; ++j) {
      int col_r = cols(k, j);

      if (IntegerVector::is_na(col_r) ||
          col_r < 1 || col_r > ncol) {
        stop("Invalid column index in row %d of cols", k + 1);
      }

      int col = col_r - 1; // R to C++ indexing

      // Stored entries in this sparse column
      for (int q = p[col]; q < p[col + 1]; ++q) {
        int r = row[q];

        if (!seen[r]) {
          seen[r] = 1;
          touched.push_back(r);
        }

        sums[r] += x[q];
      }
    }

    int count1 = 0;
    int count2 = 0;

    // Untouched rows have sum zero and cannot count as 1 or 2
    for (int r : touched) {
      if (sums[r] == 1.0) {
        ++count1;
      } else if (sums[r] == 2.0) {
        ++count2;
      }

      // Reset only rows used in this calculation
      sums[r] = 0.0;
      seen[r] = 0;
    }

    counts(k, 0) = count1;
    counts(k, 1) = count2;
  }

  return counts;
}