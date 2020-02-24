#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector getoutcome(NumericVector v1, NumericVector v2, NumericVector v3, NumericVector v4){
  int len = v1.length();
  NumericVector v5(len);
  
  int i;
  for(i=0;i<len;i++){
    if(v4[i]==0) {
      v5[i] = v1[i];
    } else if (v4[i]==1){
      v5[i] = v2[i];
    } else {
      v5[i]= v3[i];
    }
  }
  return v5;
}
