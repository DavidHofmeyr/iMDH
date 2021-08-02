#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

List smdh_subspace_cpp(arma::mat X, arma::mat V, arma::vec b, arma::vec ss, double hmult, int t_max, arma::vec clusters, arma::vec scl, arma::vec mu0, int scale, int t_init0, double C){
  int n = arma::size(X)[0];
  int d = arma::size(X)[1];
  arma::vec x(d);
  int dp = arma::size(V)[1];
  arma::mat I(dp,dp);
  //for(int j=0; j<dp; j++) I(j,j) = 1.0;
  arma::vec p(1);
  double h;
  int t;
  double div;
  double nv;
  double dstmin, dst;
  int cix;
  for(int tt = 0; tt < n; tt++){
    if(scale==1){
      t = tt + t_init0;
      if(t > t_max) t = t_max;
      mu0 = t/(t+1.0)*mu0 + X.row(tt).t()/(t+1.0);
      scl = t/(t+1.0)*scl + pow(X.row(tt).t(), 2)/(t+1.0);
      for(int j = 0; j < d; j++){
        div = scl[j] - mu0[j]*mu0[j];
        if(div > 0) x[j] = (X(tt,j) - mu0[j])/pow(div, 0.5);
        else x[j] = (X(tt,j) - mu0[j]);
      }
    }
    else{
      x = X.row(tt).t();
    }
    dstmin = pow(10.0, 20);
    for(int j = 0; j < dp; j++){
      t = tt + t_init0;
      if(t > t_max) t = t_max;
      nv = norm(V.col(j), 2);
      V.col(j) /= nv;
      p = x.t()*V.col(j);
      dst = pow(b[j]-p[0], 2);
      if(dst < dstmin){
        cix = j;
        dstmin = dst;
      }
      ss[j] = ss[j]*t/(t+1.0) + p[0]*p[0]/(t+1.0);
      h = pow(ss[j], .5)*hmult/pow(t+1.0, .2);
      V.col(j) += (b[j]-p[0])*exp(-pow(b[j]-p[0], 2)/2.0/h/h)*x/(t+1.0)/h/h/h;
      b[j] -= (b[j]-p[0])*exp(-pow(b[j]-p[0], 2)/2.0/h/h)/(t+1.0)/h/h/h;
    }
    clusters[tt] = cix;
    I = V.t()*V;
    for(int j = 0; j < dp; j++){
      I(j,j) -= 1.0;
      V.col(j) -= C*V*I.col(j)/(1+1.0);
    }
    //V += C*V*(V.t()*V - I)/(t+1.0);
  }
  List ret(8);
  ret[0] = V;
  ret[1] = b;
  ret[2] = ss;
  ret[3] = clusters;
  ret[4] = scl;
  ret[5] = mu0;
  ret[6] = n + t_init0;
  ret[7] = V.t()*V - I;
  return(ret);
}

// [[Rcpp::export]]

arma::mat srenyi_cpp(arma::mat X, arma::mat V, double hmult, int t_max, int t_init0, double C){
  int n = arma::size(X)[0];
  int d = arma::size(X)[1];
  arma::rowvec x(d);
  int dp = arma::size(V)[1];
  arma::vec p(1);
  arma::mat I(dp,dp);
  double h;
  int t;
  double nv;
  for(int tt = 0; tt < n-1; tt++){
    x = X.row(tt) - X.row(tt+1);
    for(int j = 0; j < dp; j++){
      t = tt + t_init0;
      if(t > t_max) t = t_max;
      nv = norm(V.col(j), 2);
      V.col(j) /= nv;
      p = x*V.col(j);
      h = hmult/pow(t+1.0, .2);
      V.col(j) -= p[0]*exp(-pow(p[0], 2)/2.0/h/h)*x.t()/(t+1.0)/h/h/h;
    }
    I = V.t()*V;
    for(int j = 0; j < dp; j++){
      I(j,j) -= 1.0;
      V.col(j) -= C*V*I.col(j)/(1+1.0);
    }
  }
  return(V);
}


// [[Rcpp::export]]

arma::mat sica_cpp(arma::mat X, arma::mat V, int t_max, int t_init0, double C){
  int n = arma::size(X)[0];
  int d = arma::size(X)[1];
  arma::rowvec x(d);
  int dp = arma::size(V)[1];
  arma::vec c(dp);
  arma::vec p(1);
  arma::mat I(dp,dp);
  double h;
  int t;
  double nv;
  //double Gu = sqrt(2.0/3.141593);
  for(int tt = 0; tt < n; tt++){
    x = X.row(tt);
    for(int j = 0; j < dp; j++){
      t = tt + t_init0;
      if(t > t_max) t = t_max;
      nv = norm(V.col(j), 2);
      V.col(j) /= nv;
      p = x*V.col(j);
      c[j] = c[j]*t/(t+1.0) + (1-2*p[0]*p[0])*exp(-p[0]*p[0]/2)/(t+1.0);
      if(c[j]>0) V.col(j) -= p[0]*exp(-p[0]*p[0]/2)*x.t()/(t+1.0);
      else V.col(j) += p[0]*exp(-p[0]*p[0]/2)*x.t()/(t+1.0);
      //if(p[0] >= 0) V.col(j) += 2*(p[0]-Gu)*x.t()/(t+1.0);
      //else V.col(j) += 2*(p[0]+Gu)*x.t()/(t+1.0);
      //if(p[0] >= 0) V.col(j) += x.t()/(t+1.0);
      //else V.col(j) -= x.t()/(t+1.0);
    }
    I = V.t()*V;
    for(int j = 0; j < dp; j++){
      I(j,j) -= 1.0;
      V.col(j) -= C*V*I.col(j)/(1+1.0);
    }
  }
  return(V);
}


