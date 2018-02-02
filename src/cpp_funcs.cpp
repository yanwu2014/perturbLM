// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <stdexcept>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector omitNaCpp(NumericVector x){
	std::vector<double> r(x.size());
	int k = 0;
  	for (int i = 0; i < x.size(); ++i) {
    	if (x[i]==x[i]) {
    		r[k] = x[i];
    		k++;
   		}
  	}
 	r.resize(k);
 	return(Rcpp::wrap(r));
}


// [[Rcpp::export]]
NumericVector sortCpp(NumericVector v) {
	std::sort(v.begin(), v.end());
	return(v);
}


// [[Rcpp::export]]
double calcPvalLessCpp(NumericVector v, double x) {
	if (Rcpp::NumericVector::is_na(x)) {
		return NA_REAL;
	}

	if (v.size() == 0) {
		return NA_REAL;
	}

	int num_vals_less = std::lower_bound(v.begin(), v.end(), x) - v.begin() + 1;
	int l = v.size() + 1;
	double p_val = double(num_vals_less)/l;
	return(p_val);
}


// [[Rcpp::export]]
double calcPvalGreaterCpp(NumericVector v, double x) {
	if (Rcpp::NumericVector::is_na(x)) {
		return NA_REAL;
	}

	if (v.size() == 0) {
		return NA_REAL;
	}

	int num_vals_greater = v.end() - std::lower_bound(v.begin(), v.end(), x) + 1;
	//int num_vals_greater = v.end() - std::upper_bound(v.begin(), v.end(), x) + 1;
	int l = v.size() + 1;
	double p_val = double(num_vals_greater)/l;
	return(p_val);
}


// [[Rcpp::export]]
SEXP winsorizeMatrix(SEXP Mat, SEXP Trim){
  arma::mat m=Rcpp::as<arma::mat>(Mat);
  int n=m.n_cols; int k=m.n_rows;
  int ntr=round(n * Rcpp::as<double>(Trim)); // number of positions to trim (from each side)
  if(ntr==0) { return wrap(m); } // nothing needs to be done
  for(int i=0;i<k;i++) { // for every row of the matrix
    arma::rowvec z= m.row(i);
    // determine outliers
    // arma::urowvec o=sort_index(abs(z-median(z)),1);
    arma::ucolvec o=sort_index(z); // ascending order
    // determine range
    double minv=z(o(ntr));
    double maxv=z(o(n-ntr-1));
    for(int j=0;j<ntr;j++) {
      z(o(j))=minv;
      //R_CheckUserInterrupt();
    }
    for(int j=n-ntr;j<n;j++) {
      z(o(j))=maxv;
      //R_CheckUserInterrupt();
    }
    m.row(i)=z;
    R_CheckUserInterrupt();
  }
  return wrap(m);
}


