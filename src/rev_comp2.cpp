#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Function to compute the complement of a single DNA string
std::string complement(const std::string& dna) {
  std::string result;
  result.reserve(dna.size());

  for (char nucleotide : dna) {
    switch (nucleotide) {
    case 'A': result += 'T'; break;
    case 'T': result += 'A'; break;
    case 'C': result += 'G'; break;
    case 'G': result += 'C'; break;
    default: result += 'N'; // Handle invalid characters
    }
  }
  return result;
}

// Function to compute the reverse of a single string
std::string reverse(const std::string& dna) {
  std::string result = dna;
  std::reverse(result.begin(), result.end());
  return result;
}

// [[Rcpp::export]]
std::vector<std::string> revcomp_rcpp2(std::vector<std::string> dna_strings, std::string mode = "both") {
  std::vector<std::string> result;
  result.reserve(dna_strings.size());

  for (const std::string& dna : dna_strings) {
    if (mode == "reverse") {
      result.push_back(reverse(dna));
    } else if (mode == "complement") {
      result.push_back(complement(dna));
    } else if (mode == "both") {
      result.push_back(reverse(complement(dna)));
    } else {
      stop("Invalid mode. Use 'reverse', 'complement', or 'both'.");
    }
  }

  return result;
}
