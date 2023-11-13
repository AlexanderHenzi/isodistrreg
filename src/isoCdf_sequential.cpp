#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List isoCdf_sequential(NumericVector w, NumericVector W, NumericVector Y, 
 IntegerVector posY, NumericVector y) {
  // Input
  int m = w.size();
  int mY = y.size();
  double yMax = y[mY - 1];
  
  
  // yInd is used to check when to store the CDF at y[yInd]
  int yInd = 0;
  y.push_back(R_PosInf);
  while (y[yInd] < Y[0]) yInd++;

  // Prepare object containers
  NumericVector z0 (m);
  IntegerVector PP = rep(-1, m + 1);
  NumericVector WW (m + 1);
  NumericVector MM (m + 1);
  NumericVector CDF = rep(0.0, m * mY);
  NumericVector tmpCDF;
  int lower = m * yInd;
  int upper;
  int len;

  // First iteration
  int CA = m;
  int CM = m - 1;
  int d = 1;
  int j0 = posY[0] - 1;
  z0[j0] = W[0] / w[j0];
  PP[1] = j0;
  WW[1] = sum(w[Range(0, j0)]);
  MM[0] = R_PosInf;
  MM[1] = W[0] / WW[1];
  if (j0 < m - 1) {
    d++;
    PP[2] = m - 1;
    WW[2] = sum(w[Range(j0, m - 1)]);
    CM--;
  }

  // Store results if needed; keep track of how many results have been stored
  if (Y[0] < Y[1]) {
    if (y[yInd] < Y[1]) {
      for (int l = 1; l <= d; l++) {
        len = PP[l] - PP[l - 1];
        upper = lower + len - 1;
        CDF[Range(lower, upper)] = rep(MM[l], len);
        lower = upper + 1;
      }
      yInd++;
      tmpCDF = CDF[Range(lower - m, lower - 1)];
      while (y[yInd] < Y[1]) {
        upper = lower + m - 1;
        CDF[Range(lower, upper)] = tmpCDF;
        yInd++;
        lower = upper + 1;
      }
    }
  }

  // Prepare objects for loop
  int k = 1;
  int b0;
  int a0;
  int s0;
  int dz;
  IntegerVector remPP;
  NumericVector remWW;
  NumericVector remMM;

  while (Y[k] <= yMax) {
    // Update z0
    j0 = posY[k] - 1;
    z0[j0] = z0[j0] + W[k] / w[j0];

    // Find partition to which j0 belongs
    s0 = 0;
    while (j0 > PP[s0]) s0++;
    a0 = PP[s0 - 1] + 1;
    b0 = PP[s0];
    dz = d;
    
    // Update number of operations
    CA = CA + d - s0 + b0 - a0;
    CM = CM + d; 

    // Copy tail of vector
    if (b0 < m - 1) {
      remPP = PP[Range(s0 + 1, dz)];
      remWW = WW[Range(s0 + 1, dz)];
      remMM = MM[Range(s0 + 1, dz)];
    }

    // Update value on new partition
    d = s0;
    PP[s0] = j0;
    WW[s0] = sum(w[Range(a0, j0)]);
    MM[s0] = sum(w[Range(a0, j0)] * z0[Range(a0, j0)]) / WW[s0];

    // Pooling
    while (MM[d-1] <= MM[d]) {
      d--;
      MM[d] = WW[d] * MM[d] + WW[d + 1] * MM[d + 1];
      WW[d] = WW[d] + WW[d + 1];
      MM[d] = MM[d] / WW[d];
      PP[d] = PP[d + 1];
    }

    // Add new partitions, pool
    if (j0 < b0) {
      for (int i = j0 + 1; i <= b0; i++) {
        d++;
        PP[d] = i;
        WW[d] = w[i];
        MM[d] = z0[i];
        while (MM[d - 1] <= MM[d]) {
          d--;
          MM[d] = WW[d] * MM[d] + WW[d + 1] * MM[d + 1];
          WW[d] = WW[d] + WW[d + 1];
          MM[d] = MM[d]/WW[d];
          PP[d] = PP[d + 1];
        }
      }
    }

    // Copy (if necessary)
    if (b0 < m - 1) {
      int l0 = dz - s0;
      IntegerVector ss = Range(d + 1, d + l0);
      PP[ss] = remPP;
      MM[ss] = remMM;
      WW[ss] = remWW;
      d = d + l0;
    }
    
    // Update CM, increase k
    CM = CM + b0 - j0 + 1 - d;
    k++;

    // Store in matrix (if necessary)
    if (Y[k - 1] < Y[k]) {
      while (y[yInd] < Y[k - 1]) yInd++;
      if (y[yInd] < Y[k]) {
        for (int l = 1; l <= d; l++) {
          len = PP[l] - PP[l - 1];
          upper = lower + len - 1;
          CDF[Range(lower, upper)] = rep(MM[l], len);
          lower = upper + 1;
        }
        yInd++;
        tmpCDF = CDF[Range(lower - m, lower - 1)];
        while (y[yInd] < Y[k]) {
          upper = lower + m - 1;
          CDF[Range(lower, upper)] = tmpCDF;
          yInd++;
          lower = upper + 1;
        }
      }
    }
    
    // Check for user interruption in R
    Rcpp::checkUserInterrupt();
  }
  
  // Transform vector 'CDF' to matrix
  if (y[mY - 1] >= Y[Y.size() - 1]) {
    len = m * mY - lower;
    CDF[Range(lower, m * mY - 1)] = rep(1.0, len);
  }
  CDF.attr("dim") = Dimension(m, mY);

  return List::create(_["CDF"] = CDF, _["CA"] = CA, _["CM"] = CM);
}