#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// for calculating time interval
#include <Rcpp/Benchmark/Timer.h>

#include <Rcpp.h>
#include <math.h>

using namespace arma;
using namespace std;
using namespace Rcpp;




// // [[Rcpp::export]]
// int tmp(arma::mat Y, arma::mat X, arma::vec pvec, arma::vec lam, arma::mat A, int nrank, double conv, int maxiter) {
// 
//   int p=X.n_cols;  //
//   int n=Y.n_rows, q=Y.n_cols;
// 
//   Rcpp::Rcout << "A:" << A << "\n";
// 
//   printf("n=%d, p=%d, q=%d \n", n, p, q);
// 
//   printf("the 1st element of pvec = %f \n", pvec(0));
//   printf("the 1st element of A = %f \n", A(0,0));
// 
// 
//   return pvec.size();
// }



// // [[Rcpp::export]]
// Rcpp::List which2(arma::vec x, int value) {
// 
//   Rcpp::List out;
// 
//   int nx=x.size();
//   std::vector<int> y;
// 
//   for(int i=0; i<nx; i++){
//     if( x(i) == value ) y.push_back(i+1);
//   }
// 
// 
//   Rcpp::Rcout << y.size() <<"\n";
// 
//   out["which"] = y;
//   return out;
// }



// // [[Rcpp::export]]
// int aa(arma::vec x, int y) {
// 
//   if( x.size() != y ) printf("abc \n");
// 
//   Rcpp::Rcout << arma::eye(3,3) << "\n";
// 
//   return 1;
// }



// // [[Rcpp::export]]
// arma::mat aa(arma::mat x, arma::uvec y) {
//
//
//   return x.submat(y, y) / 3;
// }







// // [[Rcpp::export]]
// arma::mat MatMul(const arma::mat& X, const arma::mat& Y) {
//   return X * Y;
// }


// arma::mat CrossProd(const arma::mat& X) {
//   return X.t() * X;
// }





// // [[Rcpp::export]]
// Rcpp::List UpdateB_SGM(arma::mat X, arma::mat Y, arma::mat A, arma::mat B,
//                    arma::vec pvec, int nrank, int d,
//                    arma::mat lam, int maxit, double eps){
// 
//   int i, k, it=0;
//   int n = Y.n_rows;
//   double lambda_ik;
//   arma::mat Xi, Xmi, Bmi, Rmik, Sik, SXiXi, SXiXi_inv, Bold;
//   arma::uvec kuvec;
//   arma::vec diff(maxit+1);
//   diff.fill(2*eps);
//   arma::uvec pi, pmi;
//   Rcpp::List out;
//   arma::vec range_i = linspace<arma::vec>(1, d, d);
//   arma::vec range_k = linspace<arma::vec>(1, nrank, nrank);
//   Timer timer;
// 
// 
//   // arma::mat zeros(10, 1);
//   // zeros.fill(0);
// 
// 
//   while ((it < maxit) & (diff(it) > eps)) {
//     it += 1;
//     Bold = B;
// 
//     for (i = 0; i < d; i++) {
// 
//       pi = find(pvec == i+1);
//       pmi = find(pvec != i+1);
// 
//       arma::mat I_pi = arma::eye(pi.size(), pi.size());
// 
//       Xi = X.cols(pi);
//       Xmi = X.cols(pmi);
//       Bmi = B.rows(pmi);
// 
//       for (k = 0; k < nrank; k++) {
//         lambda_ik = lam(i,k); // * pow(pi.size(),0.5) ;
// 
//         if(it==0 & i==0 & k==0) timer.step("start");
//         if( pmi.size() != 0 ){
//           Rmik = (Y * A.col(k) - Xmi * Bmi.col(k));
//           Sik = Xi.t() * Rmik / lambda_ik;
//         } else {
//           Rmik = Y * A.col(k);
//           Sik = Xi.t() * Rmik / lambda_ik;
//         }
//         if(it==0 & i==0 & k==0) timer.step("t1");      // record the first step
// 
//         kuvec = find(range_k == k+1);
// 
//         arma::mat zeros(pi.size(), 1);
//         zeros.fill(0);
// 
//         if( pow(accu(square(Sik)), 0.5) < 1 ){
// 
//           B.submat(pi, kuvec) = zeros;
// 
//         } else {
// 
//           if( pow(accu(square(B.submat(pi, kuvec))), 0.5) < 1e-10 ){
// 
//             B.submat(pi, kuvec) = zeros;
// 
//           } else {
// 
//             SXiXi = Xi.t() * Xi + I_pi * lambda_ik / pow(accu(square(B.submat(pi, kuvec))), 0.5);
//             SXiXi_inv = arma::pinv(SXiXi);
// 
//             if( it == 0 & (i+k)==0 ){
//               // Rcpp::Rcout << SXiXi_inv.submat(1,1,1,1) << "\n";
//               // Rcpp::Rcout << SXiXi_inv << "\n";
//             }
// 
//             B.submat(pi, kuvec) = SXiXi_inv * Xi.t() * Rmik;
// 
//           }
//         }
// 
//         if(it==0 & i==0 & k==0) timer.step("t2");  // record the second step
// 
//         diff(it) = pow(accu(square(B - Bold)), 0.5);
// 
//       }
// 
//     }
// 
// 
//   }
// 
// 
//   NumericVector time(timer);   //
//   for (int i=0; i<time.size(); i++) {
//     time[i] = time[i] / n;
//   }
// 
// 
//   out["B"] = B;
//   out["time"] = time;
//   out["iter"] = it;
//   out["diff"] = diff;
//   out["eps"] = eps;
//   // out["pinv"] = aa;
// 
// 
//   return out;
// }










// // [[Rcpp::export]]
// Rcpp::List TestEigen(arma::mat X){
//   Rcpp::List out;
//   arma::mat u, v;
//   arma::vec s;
// 
//   svd(u, s, v, X);
// 
//   out["u"] = u;
//   out["s"] = s.t()*s;
//   out["s2"] = pow(s, 2);
//   out["v"] = v;
//   out["smax"] = max(s.t()*s);
// 
// 
//   return out;
// }






// [[Rcpp::export]]
Rcpp::List UpdateB_BMD(arma::mat X, arma::mat Y, arma::mat A, arma::mat B,
                   arma::vec pvec, int nrank, int d,
                   arma::mat lam, int maxit, double eps, int threads=1){
  
  
// #ifdef _OPENMP
//   omp_set_num_threads( threads );
//   // REprintf("Number of threads=%i\n", omp_get_max_threads());
// #endif
  
  
  
  int i, k, it;
  int n = Y.n_rows;
  double lambda_ik;
  Rcpp::List diff_list(nrank);

  // B.fill(0);

  arma::uvec kuvec;
  arma::vec diff(maxit+1);
  arma::uvec pi, pmi;
  
  Rcpp::List out;
  arma::vec range_i = linspace<arma::vec>(1, d, d);
  arma::vec range_k = linspace<arma::vec>(1, nrank, nrank);
  Timer timer;


  // arma::mat zeros(10, 1);
  // zeros.fill(0);

  arma::mat Xi, Bi, Bk, Rk, Ak, Ri, SXiXi, SXiXi_inv, Bold, Ui;
  arma::mat u, v;
  arma::vec eta_i(d), s, s2, iter(nrank);
  double gamma_i, norm_Ri, ms2;


  for (int j=0; j<d; j++){
    pi = find(pvec == j+1);
    Xi = X.cols(pi);
    svd( u, s, v, Xi );
    s2 = pow(s, 2);
    ms2 = max( s2 );
    eta_i(j) = ms2;
  }



// #pragma omp parallel for// schedule(static)
  for (k = 0; k < nrank; k++) {
    Bk = B.col(k);
    Ak = A.col(k);
    
    it=0;
    diff.fill(2*eps);
    while ((it < maxit) & (diff(it) > eps)) {
      it += 1;
      Bold = Bk;

      for (i = 0; i < d; i++) {
        lambda_ik = lam(i,k);

        pi = find(pvec == i+1);
        pmi = find(pvec != i+1);

        arma::mat I_pi = arma::eye(pi.size(), pi.size());

        Xi = X.cols(pi);
        Bi = Bk.rows(pi);

        Rk = (Y * Ak - X * Bk);

        Ui = Xi.t() * Rk / n;
        gamma_i = (1 + 1e-6)*eta_i(i);
        Ri = ( Ui + gamma_i * Bi );

        norm_Ri = pow(accu(square(Ri)), 0.5);
        Bk.rows(pi) = (1/gamma_i) * Ri * max(0.0, 1 - lambda_ik/norm_Ri );


        if(it==0 & i==0 & k==0) timer.step("t2");  // record the second step

      }

      iter(k) = it;
      B.col(k) = Bk;
      diff(it) = pow( accu(square(Bk - Bold)),  0.5 );
    }

    diff_list(k) = diff;
  }


  NumericVector time(timer);   //
  for (int i=0; i<time.size(); i++) {
    time[i] = time[i] / n;
  }


  out["B"] = B;
  out["time"] = time;
  out["diff"] = diff_list;
  out["iter"] = iter;
  // out["pinv"] = aa;


  return out;
}





















// [[Rcpp::export]]
Rcpp::List png_quartimax(arma::mat X, int maxit = 1000, double eps = 1e-5, double al = 1.0){
  
  int r = X.n_cols;
  
  mat Q, Qold, mQ;
  Q.eye(r, r); 
  
  // arma::mat W = arma::mat(4, 4, arma::fill::ones);
  // arma::mat D = arma::mat(4, 4, arma::fill::zeros);
  // B.fill(0);
  
  mat Z = X * Q;
  mat dQ = -pow(Z, 3);
  mat G = X.t() * dQ ;
  
  arma::mat u, v;
  arma::vec s;
  
  Rcpp::List out;
  
  int it = 0;
  double diff = 10;
  
  
  while( (it < maxit) & (diff > eps) ){
    
    Qold = Q;
    
    mQ = Q - al*G;
    svd( u, s, v, mQ );
    Q = u * v.t();
    
    Z = X * Q;
    dQ = -pow(Z, 3);
    G = X.t() * dQ ;
    
    it = it + 1;
    
    diff = pow(accu(square(Q - Qold)), 0.5);
  }
  
  // Rcpp::Rcout << "it:" << it << "\n";
  
  out["loadings"] = X * Q;
  out["Th"] = Q;
  
  return out;
}




// [[Rcpp::export]]
Rcpp::List png_varimax(arma::mat X, int maxit = 1000, double eps = 1e-5, double al = 1.0){
  
  int p = X.n_rows, r = X.n_cols;
  
  mat Q, Qold, mQ, D;
  Q.eye(r, r); 
  
  arma::mat op = arma::vec(p, arma::fill::ones);
  
  mat Z = X * Q;
  mat dQ = (pow(Z,3) - ( Z * diagmat( sum(pow(Z,2), 0) ) ) / p ) / p;
  mat G = X.t() * dQ ;
  
  arma::mat u, v;
  arma::vec s;
  
  Rcpp::List out;
  
  int it = 0;
  double diff = 10;
  
  
  while( (it < maxit) & (diff > eps) ){
    
    Qold = Q;
    
    mQ = Q + al*G;
    svd( u, s, v, mQ );
    Q = u * v.t();
    
    Z = X * Q;
    mat dQ = (pow(Z,3) - ( Z * diagmat( sum(pow(Z,2), 0) ) ) / p ) / p;
    G = X.t() * dQ ;
    
    it = it + 1;
    
    diff = pow(accu(square(Q - Qold)), 0.5);
  }
  
  // Rcpp::Rcout << "it:" << it << "\n";
  
  out["loadings"] = X * Q;
  out["Th"] = Q;
  
  return out;
}










