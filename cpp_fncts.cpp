#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// good help:
// https://adv-r.hadley.nz/rcpp.html
// https://www.cplusplus.com/reference/cmath/abs/


// [[Rcpp::export]]
double cppdeg2rad(double deg) {
  double rad = (deg*3.141593)/180;
  return rad;
}

// [[Rcpp::export]]
double cpprad2deg(double rad) {
  double deg = (rad * 180)/3.141593;
  return deg;
}

// [[Rcpp::export]]
int cppGetSingleIndex(int x, int y, int nr) {
  int i = nr*(x-1)+y;
  return i;
}
// [[Rcpp::export]]
double cppMODULO(int i, int nr) {
  int y = i%nr;
  return y;
}




// NumericVector cppGetCircle(int r, int x, int y, int nr) {
//   
//   for (int i = 0; i < 5; i++) {
//     cout << i << "\n";
//   }
//   
// }
// 
// lapply(c(-r:r), function(X){
//   y=sqrt(r^2-X^2)
//   bottom <- round(VP[["y"]]-y)
//   top <- round(VP[["y"]]+y)
//   lapply(c(bottom:top), cppGetSingleIndex, x=X+round(VP[["x"]]), nr=nr_orig) %>%
//     return()
// })



// [[Rcpp::export]]
NumericVector cppGetXYIndex(int i, int nr) {
  int y = i%nr;
  if(y == 0){
    int yret = nr;
    int x = (i-y)/nr;
    NumericVector v =
      NumericVector::create(Named("x",x) , Named("y")=yret);
    return v;
  }
  else {
    int yret = y;
    int x = (i-y)/nr+1;
    NumericVector v =
      NumericVector::create(Named("x",x), Named("y")=yret);
    return v;
  }
}

// [[Rcpp::export]]
#include <cmath>
#include <cstdlib>
NumericVector cppElongateLine(double ha, double va, int xs, int ys, int zs, double l) {
  int x = round(xs+(l*cos(cppdeg2rad(ha))*cos(cppdeg2rad(va))));
  int y = round(ys+(l*sin(cppdeg2rad(ha))*cos(cppdeg2rad(va))));
  int z = round(zs+(l*sin(cppdeg2rad(va))));
  NumericVector v =
    NumericVector::create(Named("x",x), Named("y")=y , _["z"]=z);
  return v;
}



// [[Rcpp::export]]
#include <math.h>
double cppDistPts(double xs, double ys, double zs, double xe, double ye, double ze) {
  double xd = xs-xe;
  double yd = ys-ye;
  double zd = zs-ze;
  double d = sqrt(pow(xd, 2)+pow(yd, 2)+pow(zd, 2));
  return d;
}




// [[Rcpp::export]]
double cppSelectIntensity(int xs, int ys, int zs, List IMG) {
  if(zs > IMG.length()) {
    double r = R_NaN;
    return r;
  }
  else if(zs < 1){
    double r = R_NaN;
    return r;
  }
  else {
    int z = zs-1;
    NumericMatrix layer = IMG[z];
    if(xs > layer.ncol()) {
      double r = R_NaN;
      return r;
    }
    else if(xs < 1){
      double r = R_NaN;
      return r;
    }
    else if(ys > layer.nrow()) {
      double r = R_NaN;
      return r;
    }
    else if(ys < 1) {
      double r = R_NaN;
      return r;
    }
    else {
      int x = xs-1;
      int y = ys-1;
      double r = layer(y,x); 
      return r;
    }
  }
}



// [[Rcpp::export]]
#include <cmath>
#include <cstdlib>
double cppGetVA(double d, double ze, double zs) {
  // ze = destination, zs=start
  double a = cpprad2deg(asin((abs(zs-ze))/d));
  if(zs>ze) {
    double f = 0-a;
    return f;
  }
  else {
    double f = 0+a;
    return f;
  }
}

// [[Rcpp::export]]
#include <cmath>
#include <cstdlib>
double cppGetHA(double xe, double ye, double xs, double ys) {
  // z1 = destination, z2=start
  double a = cpprad2deg(atan((abs(xe-xs))/abs(ye-ys)));
    if(xs>=xe&ys>=ye){
      double f = 270-a;
      return f;
    } 
    else if(xs>=xe&ys<=ye){
      double f = 90+a;
      return f;
    } 
    else if(xs<=xe&ys>=ye){
      double f = 270+a;
      return f;
    } 
    else if(xs<=xe&ys<=ye){
      double f = 90-a;
      return f;
    }
}
// get_ha <- function(x2,y2,x1,y1){
// ## x2 = destination... x1 = begin
//   angle <- cpprad2deg(atan((abs(x1-x2))/abs(y1-y2)))
//   if(x1>=x2&y1>=y2){
//     fin <- 270-angle
//   } else if(x1>=x2&y1<=y2){
//     fin <- 90+angle
//   } else if(x1<=x2&y1>=y2){
//     fin <- 270+angle
//   } else if(x1<=x2&y1<=y2){
//     fin <- 90-angle
//   } 
//   return(fin)
// }











