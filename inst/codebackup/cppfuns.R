#                      -------------------
#--------------------- Start C++ functions --------------------
#' @title Return all (discrete) forward distances in veiw.
#'
#' @description
#'  C++ function to return all (discrete) forward distances in veiw.
#'  
#' @param Ix Perpendicular distance.
#' @param Inull_yobs REDUNDANT; must = NULL.
#' @param Iyobs Forward distance.
#' @param Itheta_f REDUNDANT; must = 0
#' @param Itheta_b REDUNDANT; must = 90
#' @param Iymax Maximum forward distance in view.
#' @param Idy Forward distance step increment.
#' 
gety.obs <- function(Ix, Inull_yobs, Iyobs, Itheta_f, Itheta_b, Iymax, Idy) {
  .Call( "bias_gety_obs", Ix, Inull_yobs, Iyobs, Itheta_f, Itheta_b, Iymax, Idy, PACKAGE = "hmltm" )
}

#' @title Calculate probability of detection at (x,y).
#'
#' @description
#'  C++ function to calculate probability of detection at (x,y), or if cdf=TRUE, BY (x,y).
#'  
#' @param Ix Perpendicular distance.
#' @param Iy Forward distance.
#' @param Ihfun Hazard function name.
#' @param Ib Hazard function parameter vector.
#' @param Ipcu Bernoulli state-dependent probability parameters.
#' @param IPi Markov model transition probability matrix.
#' @param Idelta Markov model stationary distribution.
#' @param Iymax Maximum forward distance in view.
#' @param Idy Forward distance step increment.
#' @param Itheta_f REDUNDANT; must = 0
#' @param Itheta_b REDUNDANT; must = 90
#' @param Ially Flag for whether or not to return probabilities at all forward distances in view.
#' @param Icdf Flag for whether or not to return cumulative distribution function in forward dimension.
#' This differes from specifying Ially=TRUE in that Ially=TRUE calculates the cdf from ymax 
#' to y=0, whereas Icdf=TRUE calculates the cdf from ymax to y.
p.xy1 <- function(Ix, Iy, Ihfun, Ib, Ipcu, IPi, Idelta, Iymax, Idy, Itheta_f, Itheta_b, Ially, Icdf){
  .Call( "bias_p_xy1", Ix, Iy, Ihfun, Ib, Ipcu, IPi, Idelta, Iymax, Idy, Itheta_f, Itheta_b, Ially, Icdf, PACKAGE = "hmltm" )
}

#' @title Parameter transformation function for all models.
#'
#' @description
#'  Parameter transformation for all models
#'  
#' @param b parameters on original scale.
#' @param fun detection hazard function name (character) - see details below.
#' 
#' @details
#' Valid detection hazard function names for detection hazards with certain detection at radial
#' distance zero are "h.EP1.0", "h.EP1x.0", "h.EP2.0", "h.EP2x.0", "h.IP.0".
#' #' Transformations are:
#' \describe{
#' \item{h.EP1.0}{log transform}
#' \item{h.EP1x.0}{log transform}
#' \item{h.EP2.0}{log transform}
#' \item{h.EP2x.0}{log transform}
#' \item{h.IP.0}{log transform}
#' }
tfm <- function(b, fun) {
  .Call("hmltm_get_tfm", b, fun, PACKAGE = "hmltm" )
}

#' @title Parameter inverse transformation function for all models.
#'
#' @description
#'  Parameter inverse transformation for all models.
#'  
#' @param b parameters on transformed scale
#' @param fun detection hazard function name (character) - see details below.
#'  
#' @details
#' Valid detection hazard function names for detection hazards with certain detection at radial
#' distance zero are "h.EP1.0", "h.EP1x.0", "h.EP2.0", "h.EP2x.0", "h.IP.0".
#' Inverse transformations are:
#' \describe{
#' \item{h.EP1.0}{exponential}
#' \item{h.EP1x.0}{exponential}
#' \item{h.EP2.0}{exponential}
#' \item{h.EP2x.0}{exponential}
#' \item{h.IP.0}{exponential}
#' }
invtfm <- function(b, fun) {
  .Call("hmltm_get_invtfm", b, fun, PACKAGE = "hmltm" )
}



#---------------------  End C++ functions --------------------
#                      -------------------