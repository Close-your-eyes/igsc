#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame compare_nonref_cpp(
    DataFrame df,
    std::string ref,
    std::string pos_col = "position",
    std::string seq_col = "seq",
    std::string name_col = "seq.name",
    std::string match_symbol = ".",
    std::string mismatch_symbol = "x",
    bool keep_gaps = true,
    std::string nonref_mismatch_as = "base"
) {

  IntegerVector pos = df[pos_col];
  CharacterVector seq = df[seq_col];
  CharacterVector names = df[name_col];

  int n = seq.size();

  // build reference lookup: position -> reference base
  std::unordered_map<int, std::string> ref_map;

  for (int i = 0; i < n; i++) {
    std::string nm = Rcpp::as<std::string>(names[i]);

    if (nm == ref) {
      ref_map[pos[i]] = Rcpp::as<std::string>(seq[i]);
    }
  }

  // mutate non-reference rows
  for (int i = 0; i < n; i++) {

    std::string nm = Rcpp::as<std::string>(names[i]);

    // skip reference sequence itself
    if (nm == ref)
      continue;

    int p = pos[i];

    // skip if no reference value at this position
    if (ref_map.find(p) == ref_map.end())
      continue;

    std::string v1 = Rcpp::as<std::string>(seq[i]);
    std::string v2 = ref_map[p];

    if (v1 == v2) {

      // match branch
      if (v1 != "-" || !keep_gaps) {
        seq[i] = match_symbol;
      }

    } else {

      // mismatch branch
      if (nonref_mismatch_as == "mismatch_symbol") {

        if (v1 != "-" || !keep_gaps) {
          seq[i] = mismatch_symbol;
        }

      }
    }
  }

  df[seq_col] = seq;

  return df;
}