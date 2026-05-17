#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List qual_stats_cpp(List pq) {

    int n = pq.size();

    IntegerVector minQual(n);
    NumericVector meanQual(n);
    IntegerVector n_belowQ30(n);

    // keep qualities as a list of integer vectors
    List readQualNum(n);

    for (int i = 0; i < n; i++) {

        IntegerVector x = pq[i];
        int len = x.size();

        // store original vector directly
        readQualNum[i] = x;

        if (len == 0) {
            minQual[i] = NA_INTEGER;
            meanQual[i] = NA_REAL;
            n_belowQ30[i] = NA_INTEGER;
            continue;
        }

        int mn = x[0];
        double sm = 0.0;
        int below = 0;

        for (int j = 0; j < len; j++) {

            int val = x[j];

            if (val < mn)
                mn = val;

            sm += val;

            if (val < 30)
                below++;
        }

        minQual[i] = mn;
        meanQual[i] = sm / len;
        n_belowQ30[i] = below;
    }

    return List::create(
        _["readQualNum"] = readQualNum,
        _["minQual"] = minQual,
        _["meanQual"] = meanQual,
        _["n_belowQ30"] = n_belowQ30
    );
}