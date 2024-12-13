#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Function to compute the reverse complement of a single DNA string
std::string reverse_complement(const std::string& dna) {
    std::string complement;
    complement.reserve(dna.size());

    for (char nucleotide : dna) {
        switch (nucleotide) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            default: complement += 'N'; // Handle invalid characters
        }
    }

    // Reverse the complement string
    std::reverse(complement.begin(), complement.end());
    return complement;
}

// [[Rcpp::export]]
std::vector<std::string> revcomp_rcpp(std::vector<std::string> dna_strings) {
    std::vector<std::string> result;
    result.reserve(dna_strings.size());

    for (const std::string& dna : dna_strings) {
        result.push_back(reverse_complement(dna));
    }

    return result;
}