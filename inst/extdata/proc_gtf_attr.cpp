#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List process_attr_col_rcpp(std::vector<std::string> x) {
  List result(x.size());

  for (size_t i = 0; i < x.size(); ++i) {
    std::vector<std::string> current_result;
    std::string current_str = x[i];
    size_t start = 0, end;

    // Split by semicolon
    while ((end = current_str.find(';', start)) != std::string::npos) {
      std::string attr = current_str.substr(start, end - start);
      start = end + 1;

      // Remove double quotes
      attr.erase(std::remove(attr.begin(), attr.end(), '"'), attr.end());

      // Trim spaces
      size_t first = attr.find_first_not_of(" \t");
      size_t last = attr.find_last_not_of(" \t");
      if (first != std::string::npos)
        attr = attr.substr(first, last - first + 1);

      // Split by space
      size_t sub_start = 0, sub_end;
      while ((sub_end = attr.find(' ', sub_start)) != std::string::npos) {
        current_result.push_back(attr.substr(sub_start, sub_end - sub_start));
        sub_start = sub_end + 1;
      }
      current_result.push_back(attr.substr(sub_start));
    }

    // Store results for current string
    result[i] = current_result;
  }

  return result;
}
