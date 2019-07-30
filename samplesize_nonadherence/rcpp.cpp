#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector getoutcome(NumericVector v1, NumericVector v2, NumericVector v3){
  int len = v1.length();
  NumericVector v4(len);
  
  int i;
  for(i=0;i<len;i++){
    if(v3[i]==1) {
      v4[i] = v1[i];
    } else {
      v4[i]= v2[i];
    }
  }
  return v4;
}
