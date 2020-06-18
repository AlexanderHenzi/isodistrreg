#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pavaCorrect(NumericMatrix y) {
    /* Apply the Pool-Adjacent Violators Algorithm (PAVA) to the rows of a
     * numeric matrix y. This code is adapted from the R code by Lutz Duembgen
     * (2008) at
     * http://www.imsv.unibe.ch/about_us/files/lutz_duembgen/software/index_eng.html
     */
    int n = y.nrow();
    int m = y.ncol();
    NumericMatrix out(n, m);
    NumericVector weight(m);
    IntegerVector index(m);
    int ci = 0, j = 0, i;
    double nw;
    for (int k = 0; k < n; k++) {
        index[ci] = 0;
        weight[ci] = 1;
        out(k, ci) = y(k, 0);
        
        while (j < m - 1) {
            j += 1;
            ci += 1;
            index[ci] = j;
            weight[ci] = 1;
            out(k, ci) = y(k, j);
            while (ci >= 1 && out(k, ci) <= out(k, ci - 1)) {
                nw = weight[ci - 1] + weight[ci];
                out(k, ci - 1) = out(k, ci - 1) +
                    (weight[ci] / nw) * (out(k, ci) - out(k, ci - 1));
                weight[ci - 1] = nw;
                ci -= 1;
            }
        }
        
        while (j >= 0) {
            for (i = index[ci]; i <= j; i++) {
                out(k, i) = out(k, ci);
            }
            j = index[ci] - 1;
            ci -= 1;
        }
        
        ci = 0;
        j = 0;
        Rcpp::checkUserInterrupt();
    }
     
    return out;
}