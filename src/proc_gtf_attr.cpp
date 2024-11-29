#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List process_attr_col_rcpp(std::vector<std::string> x) {
  List result(x.size()); // Initialize result as an R List with the size of input vector

  for (size_t i = 0; i < x.size(); ++i) {
    std::vector<std::string> current_result; // Holds processed parts for the current string
    std::string current_str = x[i];
    
    // Remove the last semicolon if it exists
    //if (!current_str.empty() && current_str.back() == ';') {
     // current_str.pop_back(); // Remove the last character
    //}

    size_t start = 0, end;

    // Split by "; "
    while ((end = current_str.find("; ", start)) != std::string::npos) {
      std::string attr = current_str.substr(start, end - start);
      start = end + 2; // Advance past "; "


      // Split at the first space only
      size_t first_space = attr.find(' ');
      if (first_space != std::string::npos) {
        // Part before the space
        std::string key = attr.substr(0, first_space);
        // Part after the space
        std::string value = attr.substr(first_space + 1);

        // Remove double quotes
        key.erase(std::remove(key.begin(), key.end(), '"'), key.end());
        value.erase(std::remove(value.begin(), value.end(), '"'), value.end());

        // Push results
        current_result.push_back(key);
        current_result.push_back(value);
      } else {
        // If no space is found, treat the whole string as a key
        attr.erase(std::remove(attr.begin(), attr.end(), '"'), attr.end());

        // Trim whitespaces
        size_t first = attr.find_first_not_of(" \t");
        size_t last = attr.find_last_not_of(" \t");
        if (first != std::string::npos)
          attr = attr.substr(first, last - first + 1);

        current_result.push_back(attr);
      }
    }

    // Store results for the current string
    result[i] = current_result;
  }

  return result; // Return the R List
}

