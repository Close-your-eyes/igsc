#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector mutate_value_cpp(
    CharacterVector value,
    CharacterVector ref_col,
    std::string match_symbol,
    std::string mismatch_symbol,
    bool keep_match_gaps,
    std::string nonref_mismatch_as
) {
  int n = value.size();
  CharacterVector out(n);
  
  for (int i = 0; i < n; i++) {
    std::string v = as<std::string>(value[i]);
    std::string r = as<std::string>(ref_col[i]);
    
    if (v == r) {
      // match branch
      if (v != "-") {
        out[i] = match_symbol;
      } else {
        out[i] = keep_match_gaps ? "-" : match_symbol;
      }
    } else {
      // mismatch branch
      if (nonref_mismatch_as == "base") {
        out[i] = v;
      } else {
        if (v != "-") {
          out[i] = mismatch_symbol;
        } else {
          out[i] = keep_match_gaps ? "-" : mismatch_symbol;
        }
      }
    }
  }
  
  return out;
}
