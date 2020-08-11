
#                          -----------------------------------------------
#------------------------- Start Analytic Variance Estimation functions -----------------------

#' @title Calculates asymptotic variance using Hessian.
#'
#' @description
#' Calculates asymptotic variance of a statistic derived from a hsltm fit, using the inverse of the 
#' Hessian matrix and estimated derivative of the statististic with respect to the parameter vector.
#' The derivative is approximated using function \code{\link{splinefun}}.
# 
#' @param stat name of statistic to calculate. Valid statistics are "p0" for estimated probability 
#' at perpendicular distance zero, "p" for mean estimated detection probability over all perpendicular
#' distances, "invp" for the inverse of mean estimated detection probability over all perpendicular
#' distances, "esw" for estimated effective strip width, and "invesw" for estimated inverse of effective 
#' strip width.
#' @param hmmlt object output by \code{\link{fit.hmltm}}.
#' @param nx4spline is number of points at which to evaluate function for spline interpolation.
#' @param doplot if TRUE, plots the spline approximation to the statistic as a function of each 
#' parameter.
#' 
#' @seealso \code{\link{hmltm.stat}} for calculation of derived statistics, and 
#' \code{\link{splinefun}} for approximation of function derivatives.
#' 
invHessvar=function(stat,hmmlt,nx4spline=50,doplot=FALSE){
  # extract things we need:
  models=hmmlt$models
  if(!is.nullmodel(models)) stop("This function only coded for models without covariates.")
  cov=hmmlt$xy[1,] # extract only 1 obs 'cause assuming no covars, so stats same for all obs
  hfun=hmmlt$h.fun
  b=hmmlt$fit$b
  hmm.pars=hmmlt$fitpars$hmm.pars
  survey.pars=hmmlt$fitpars$survey.pars
  vcv=solve(hmmlt$fit$hessian)
  vlen=dim(vcv)[1]
  # calculate derivative vector
  dstat.db=rep(0,vlen)
  no.se=FALSE
  for(i in 1:vlen){
    se=sqrt(vcv[i,i])
    if(is.nan(se)) {
      no.se=TRUE
      x=seq(b[i]-b[i]/10,b[i]+b[i]/10,length=nx4spline)
    } else {x=seq(b[i]-se/2,b[i]+se/2,length=nx4spline)}
    y=x*0
    for(j in 1:nx4spline) {
      bi=b
      bi[i]=x[j]
      y[j]=hmltm.stat(stat,bi,hfun,models,cov,survey.pars,hmm.pars)
    }
    sf=splinefun(x,y)
    if(doplot) {
      plot(x,y,type="l")
      lines(x,sf(x),col="red",lty=2)
    }
    dstat.db[i]=sf(b[i],deriv=1)
  }
  # use inverse Hessian to calculate variance of statistic:
  if(no.se) return(list(invHess=vcv,dstat.db=dstat.db,stat=stat,se=-999,cv=-999))
  else {
    VAR=t(dstat.db%*%vcv%*%dstat.db)
    SE=sqrt(VAR)
    stat=hmltm.stat(stat,b,hfun,models,cov,survey.pars,hmm.pars) # called with original b parameters
    CV=SE/stat
    return(list(invHess=vcv,dstat.db=dstat.db,stat=stat,se=SE,cv=CV))
  }
}
