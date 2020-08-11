#ifndef _hsltm_RCPP_HELLO_WORLD_H
#define _hsltm_RCPP_HELLO_WORLD_H

#include <RcppArmadillo.h>

RcppExport SEXP bias_p_xy1(SEXP Ix, SEXP Iy, SEXP Ihfun, SEXP Ib, SEXP Ipcu, SEXP IPi, SEXP Idelta, SEXP Iymax, 
                           SEXP Idy, SEXP Itheta_f, SEXP Itheta_b, SEXP Ially, SEXP Icdf);


RcppExport SEXP bias_gety_obs(SEXP Ix, SEXP Inull_yobs, SEXP Iyobs, SEXP Itheta_f,
                              SEXP Itheta_b, SEXP Iymax, SEXP Idy);

RcppExport SEXP hsltm_get_tfm(SEXP Istype, SEXP Iv);

RcppExport SEXP hsltm_get_invtfm(SEXP Istype, SEXP Iv);

#endif
