#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector tpr_vec(List splitted_classes, NumericVector fpr, NumericVector tpr, NumericVector interpoint, NumericVector ncontrols, double N){   
  int m = fpr.size();
  int nout = interpoint.size();
  int n = m;
  NumericVector fpr_weighted(m);
  NumericVector fpr_weight(m);
  NumericVector tpr_weight(m);
  NumericVector tpr_interpolated(nout);
  NumericVector sum_tpr_fpr(m);
  
  
  for (int k = 0; k < m; k++){
    fpr_weighted[k] = fpr[k] / ncontrols[0];
    tpr_weight[k] = tpr[m-k-1]; 
    fpr_weight[k] = fpr[m-k-1];
  }
  
  for (int r = 0; r < nout; r++){
    double v = interpoint[r];
    int i1 = 0;
    int j1 = n-1;
    int ij; 
    while(i1 < j1-1){
      ij = (i1+j1)/2;
      if (v < fpr_weighted[ij]) {
        j1 = ij;
      } else {
        i1 = ij;
      }
    }
    
    if(v == fpr_weighted[j1]){ 
      tpr_interpolated[r] = tpr[j1]*ncontrols[0];
      }
      
    else if(v == fpr_weighted[i1]){ 
    tpr_interpolated[r] = tpr[i1]*ncontrols[0];
    }
    else {
      tpr_interpolated[r] = (tpr[i1] + (tpr[j1] - tpr[i1]) * ((v - fpr_weighted[i1])/(fpr_weighted[j1] - fpr_weighted[i1])))*ncontrols[0];
    }
    
  }
  
  for(int k = 0; k < m; k++){
    sum_tpr_fpr[k] = tpr_weight[k] + fpr_weight[k];
  }
  
  
  
  for (int i = 1; i < (N-1); i++){
    NumericVector sum_indicator(m);
    SEXP list_element_0 = splitted_classes[i];
    Rcpp::NumericVector list_element(list_element_0);
    std::sort( list_element.begin(), list_element.end() );
    int L = list_element.size();
    int rep_val = list_element[0]; 
    for(int l = 0; l < rep_val; l++){
        sum_indicator[l] = L;
      }
    int indx = rep_val;
    for (int j = 1; j < L; j++){
      int rep_val = list_element[j] - list_element[j-1];
      for(int l = 0; l < rep_val; l++){
        sum_indicator[indx + l] = L-j;
      }
      indx = indx + rep_val;
    }
    
    
    for (int k = 0; k < indx; k++){
      tpr_weight[k] = tpr_weight[k] - sum_indicator[k];
      fpr[k] = (sum_tpr_fpr[k] - tpr_weight[k])/ncontrols[i];
    }
    for (int k = indx; k < m; k++) {
      fpr[k] = (sum_tpr_fpr[k] - tpr_weight[k])/ncontrols[i];
    }
    
    
    
    NumericVector rev_fpr(m);
    NumericVector rev_tpr_weight(m); 
    
    for(int k = 0; k < m; k++){
      rev_fpr[k] = fpr[m-k-1];
      rev_tpr_weight[k] = tpr_weight[m-k-1];
    }
    
    for (int r = 0; r < nout; r++){
    double v = interpoint[r];
    int i1 = 0;
    int j1 = n-1;
    int ij; 
    while(i1 < j1-1){
      ij = (i1+j1)/2;
      if (v < rev_fpr[ij]) {
        j1 = ij;
      } else {
        i1 = ij;
      }
    }
    
    if(v == rev_fpr[j1]){ 
      tpr_interpolated[r] = rev_tpr_weight[j1]*ncontrols[i] + tpr_interpolated[r];
      }
      
    else if(v == rev_fpr[i1]){ 
    tpr_interpolated[r] = rev_tpr_weight[i1]*ncontrols[i] + tpr_interpolated[r];
    }
    else {
      tpr_interpolated[r] = (rev_tpr_weight[i1] + (rev_tpr_weight[j1] - rev_tpr_weight[i1]) * ((v - rev_fpr[i1])/(rev_fpr[j1] - rev_fpr[i1])))*ncontrols[i] + tpr_interpolated[r];
    }
    
  }
// Check for user interruption in R
    Rcpp::checkUserInterrupt();

  }
  
  return tpr_interpolated;
}