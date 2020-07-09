#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix idrHazardCpp(NumericVector w, NumericVector W, NumericVector Y, 
 IntegerVector posY, NumericVector y) {
  int m = w.size();
  int ny = y.size();
  NumericMatrix out(m, ny);
  
  NumericVector denom = clone(w); // denominator in hazard ratio fraction
  NumericVector enume = clone(w); // enumerator in hazard ratio fraction
  NumericVector ee (m + 1); // enumerator in partition
  NumericVector dd (m + 1); // denominator in partition
  NumericVector ff (m + 1); // hazard ratio in partition
  IntegerVector pp (m + 1); // last index of partition (0-based)
  ff[0] = -1; // (no verification of cc > 0 needed).
  pp[0] = -1;
  int l = 0; // run through all Y-values to update enum
  int count = 0; // to copy result into columns of 'out'
  int cc = 1; // counter in loop
  
  double maxy = y[y.size() - 1];
  
  while (Y[l] < maxy) {
    int j = posY[l] - 1;
    enume[j] = enume[j] - W[l]; // if !(Y[l] > y[h]), drop the weight of Y[l]
    l++;
    if (Y[l - 1] == Y[l]) continue; // continue until Y[l - 1] < Y[l]
    
    // PAVA
    ee[1] = enume[0];
    dd[1] = denom[0];
    if (dd[1] == 0) {
      ff[1] = 0; // if 0/0, set to 0 ("smaller")
    } else {
      ff[1] = ee[1] / dd[1];
    }
    pp[1] = 0;

    for (int i = 1; i < m; i++) {
      cc++;
      pp[cc] = i;
      ee[cc] = enume[i];
      dd[cc] = denom[i];
      if (dd[cc] == 0) {
        ff[cc] = 0;
      } else {
        ff[cc] = ee[cc] / dd[cc];
      }
      while (ff[cc] <= ff[cc - 1]) {
        cc--;
        ee[cc] = ee[cc] + ee[cc + 1];
        dd[cc] = dd[cc] + dd[cc + 1];
        if (dd[cc] == 0) {
          ff[cc] = 0;
        } else {
          ff[cc] = ee[cc] / dd[cc];
        }
        pp[cc] = pp[cc + 1];
      }
    }

    // copy to martix 'out'
    for (int r = 1; r <= cc; r++) {
      double v = ff[r];
      for (int s = pp[r - 1] + 1; s <= pp[r]; s++) {
        out(s, count) = v;
      }
    }
    count++;
    denom = clone(enume);
    cc = 1;
  }
  
  // compute survival functions from hazard ratios
  for (int nr = 0; nr < m; nr++) {
    for (int nc = 1; nc < ny; nc++) {
      out(nr, nc) = out(nr, nc - 1) * out(nr, nc);
    }
    // int nc = 1;
    // while (out(nr, nc) > 0) {
    //   out(nr, nc) = out(nr, nc - 1) * out(nr, nc);
    //   nc++;
    // }
  }

  // transform survival functions to cdfs
  out = 1 - out;
  
  return out;
}