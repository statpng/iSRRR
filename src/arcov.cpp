#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// #ifdef _OPENMP
// #include <omp.h>
// #endif

// for calculating time interval
#include <Rcpp/Benchmark/Timer.h>

// // [[Rcpp::plugins(openmp)]]



#include <Rcpp.h>
#include <math.h>

using namespace arma;
using namespace std;
using namespace Rcpp;



//[[Rcpp::export()]]
NumericMatrix arcov(int p, double rho, int threads=1, bool display_progress=true){
// #ifdef _OPENMP
//   omp_set_num_threads(threads);
// #endif
  NumericMatrix out(p);

// # pragma omp parallel for schedule(dynamic)
  for(int i=0; i<p; ++i){

    for(int j=i; j<p; ++j){
      if(j>i){
        double tmp = pow(rho, j-i);
        if( tmp > 0.001 ){
          out(i,j) = out(j,i) = tmp;
          //out(j,i) = tmp;
        }
      } else if (i==j){
        out(i,i) = 1;
      }
    }
  }
  return out;
}
