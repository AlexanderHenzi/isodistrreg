#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pavaDec(List cpY, NumericVector thresholds, NumericVector w) {
    /* Compute one-dimensional isotonic distributional regression (IDR) using
     * the Pool-Adjacent Violators Algorithm (PAVA): For two vectors x and y,
     * solve the problem
     *          minimize sum( (1{y_i <= y_k} - p_i)^2 )
     *          subject to p_i >= p_j if x_i <= x_j
     * for each k, where 1{A} = 1 if A is TRUE and 0 otherwise. If the values of
     * x are not distinct, then there are
     *      xx_1 < xx_2 < ... < xx_m
     * such that all entries x_i are contained in {xx_1, ..., xx_m}. Denote by
     * cpY[[i]] the y values belongig to xx[i]. Then the above minimizaion
     * problem is equivalent to
     *          minimize sum(w * (mean(cpY[[j]] <= y_k) - q_j)^2 )
     *          subject to q_1 >= q_2 >= ... >= q_m,
     * where w = (length(crpY[[1]]), ..., length(cpY[[m]])).
     * This can be solved by applying antitonic PAVA to the vectors
     *          (mean(cpY[[1]] <= y_k), ...., mean(cpY[[m]] <= y_k))
     * with weigths w.
     */
    int n = cpY.length();
    int m = thresholds.length() - 1;
    NumericMatrix y(n, m), out(n, m);
    NumericVector weight(n);
    IntegerVector index(n);
    int ci = 0, j = 0, i;
    double nw;
    
    // compute input vectors for PAVA
    for (int ii = 0; ii < n; ii++) {
        NumericVector tmp = cpY[ii];
        double ny = tmp.length();
        for (int jj = 0; jj < m; jj++) {
            y(ii, jj) = sum(tmp <= thresholds[jj]) / ny;
        }
    }
    
    // apply PAVA
    for (int k = 0; k < m; k ++) {
        index[ci] = 0;
        weight[ci] = w[0];
        out(ci, k) = y(0, k);
        
        while (j < n - 1) {
            j += 1;
            ci += 1;
            index[ci] = j;
            weight[ci] = w[j];
            out(ci, k) = y(j, k);
            while (ci >= 1 && out(ci, k) >= out(ci - 1, k)) {
                nw = weight[ci - 1] + weight[ci];
                out(ci - 1, k) = out(ci - 1, k) +
                    (weight[ci] / nw) * (out(ci, k) - out(ci - 1, k));
                weight[ci - 1] = nw;
                ci -= 1;
            }
            
        }
        
        while (j >= 0) {
            for (i = index[ci]; i <= j; i++) {
                out(i, k) = out(ci, k);
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