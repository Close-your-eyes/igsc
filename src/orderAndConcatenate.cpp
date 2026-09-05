#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

// Function to order entries in each row of a matrix alphabetically and concatenate them with underscores
// [[Rcpp::export]]
std::vector<std::string> orderAndConcatenateStrings(Rcpp::CharacterMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  std::vector<std::string> result(nrow);
  
  for (int i = 0; i < nrow; ++i) {
    std::vector<std::string> row_entries;
    for (int j = 0; j < ncol; ++j) {
      std::string entry = Rcpp::as<std::string>(mat(i, j));
      row_entries.push_back(entry);
    }
    std::sort(row_entries.begin(), row_entries.end());
    std::string concatenated_row;
    for (const auto& entry : row_entries) {
      concatenated_row += entry + "_";
    }
    concatenated_row.pop_back(); // Remove the last underscore
    result[i] = concatenated_row;
  }
  
  return result;
}
