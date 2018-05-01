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


