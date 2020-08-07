#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix idrHazardCpp(NumericVector w, NumericVector W, NumericVector Y,
  IntegerVector posY, NumericVector y) {
  int m = w.size();
  int k = y.size();
  NumericMatrix out(m, k);
  NumericVector ww (m + 1);
  NumericVector mm (m + 1);
  NumericVector pp (m + 1);
  NumericVector v = rep(1.0, m);
  NumericVector neww = clone(w);
  mm[0] = R_NegInf;
  pp[0] = -1;
  int l = 0;
  int j;
  int d = 1;
  double fill;
  int count = 0;
  double maxy = y[y.size() - 1];

  while (Y[l] == y[0]) {
    j = posY[l] - 1;
    v[j] = v[j] - W[l] / w[j];
    l++;
  }
  
  pp[1] = 0;
  ww[1] = w[0];
  mm[1] = v[0];
  for (int i = 1; i < m; i++) {
    d++;
    pp[d] = i;
    mm[d] = v[i];
    ww[d] = w[i];
    while (mm[d] <= mm[d - 1]) {
      d--;
      mm[d] = ww[d] * mm[d] + ww[d + 1] * mm[d + 1];
      ww[d] = ww[d] + ww[d + 1];
      mm[d] = mm[d] / ww[d];
      pp[d] = pp[d + 1];
    }
  }
  
  // no regression needed for parts which are already zero...
  int zerCount = 1;
  while (mm[zerCount] == 0) zerCount++;
  int zerInd = pp[zerCount - 1] + 1;

  for (int r = zerCount; r <= d; r++) {
    fill = mm[r];
    for (int s = pp[r - 1] + 1; s <= pp[r]; s++) {
      out(s, count) = fill;
      neww[s] = neww[s] * fill;
    }
  }
  
  while (Y[l] < maxy) {
    j = posY[l] - 1;
    v[j] = v[j] - W[l] / w[j];
    l++;
    if (Y[l - 1] == Y[l]) continue;
    
    if (zerInd > 0) {
      pp[1] = zerInd - 1;
      ww[1] = 0;
      mm[1] = 0;
      d = 1;
    } else {
      d = 0;  
    }

    for (int i = zerInd; i < m; i++) {
      // if (neww[i] > 0) {
        d++;
        pp[d] = i;
        mm[d] = v[i] / out(i, count);
        ww[d] = neww[i];
      // } else {
        // continue;
      // }
      while (mm[d] <= mm[d - 1]) {
        d--;
        mm[d] = ww[d] * mm[d] + ww[d + 1] * mm[d + 1];
        ww[d] = ww[d] + ww[d + 1];
        mm[d] = mm[d] / ww[d];
        pp[d] = pp[d + 1];
      }
    }
    zerCount = 1;
    while (mm[zerCount] == 0) zerCount++;
    zerInd = pp[zerCount - 1] + 1;
    
    count++;
    for (int r = zerCount; r <= d; r++) {
      fill = mm[r];
      for (int s = pp[r - 1] + 1; s <= pp[r]; s++) {
        out(s, count) = out(s, count - 1) * fill;
        neww[s] = neww[s] * fill;
      }
    }
    Rcpp::checkUserInterrupt();
  }
  
  return 1 - out;
}