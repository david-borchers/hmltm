#                       -----------------------------------
#---------------------- Start estimation functions for x and y ----------------------


#' @title Maximises hidden Markov line transect model likelihood.
#'
#' @description
#' \code{fit.xy} A wrapper function for \code{\link{optim}} to get maximum likelihood estimates of
#' detection hazard function paramters using forward and perpendicular distance data 
#' (and any associated covariates), using estimated Markov model or hidden Markov model availability 
#' prameters. 
#'
#' @param pars starting parameter values.
#' @param xy data frame with one line per detection containing perpendicular distance ($x) and
#' forward distance ($y) and any covariates in the model.
#' @param FUN detection hazard functional form name (character).
#' @param models list of characters with elements \code{$y} and \code{$x} specifying y- and 
#' x-covariate models. Either \code{NULL} or regression model format (without response on left).
#' @param pm is a vector of state-dependent Bernoulli distribution parameters (interpretation: the
#' probability of being available, given state).
#' @param Pi is a Markov model transition matrix governing the transition between states. (Square matrix
#' with each dimension equal to that of \code{pm}.)
#' @param delta is a vector specifying the stationary distribution of the Markov model governed by Pi.
#' @param theta.f REDUNDANT parameter determining the max forward angle in view (must=0).
#' @param theta.b REDUNDANT parameter determining the max forward angle in view (must=90).
#' @param W perpendicular truncation distance for fitting. Must be greater than \code{max(xy$x)}.
#' @param ymax maximum forward distance that things could be detected. Must be greater than \code{max(xy$y)}.
#' @param dy resolution of forward distances for Markov model (typically observer speed times time step size).
#' @param nx NOT SURE - I think it is redundant - CHECK
#' @param hessian Logical; if TRUE Hessian matrix is returned
#' @param control list controlling \code{\link{optim}} optimisation. 
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
#'
#' @return A list comprising \code{\link{optim}} ouptut plus element $par containing the estimated parameters 
#'  on the same scale as input parameters pars (which are transformed before calling \code{\link{optim}}).
fit.xy=function(pars,xy,FUN,models=list(y=NULL,x=NULL),pm,Pi,delta=delta,
                theta.f=0,theta.b=90,W,ymax,dy,nx=50,hessian=FALSE,
                control=list(trace=5,reltol=1e-6,maxit=200),groupfromy=NULL)
{
  b=tfm(pars,FUN)
  fit=optim(par=b,fn=negllik.xy,control=control,hessian=hessian,xy=xy,
            FUN=FUN,models=models,pm=pm,Pi=Pi,delta=delta,
            theta.f=theta.f,theta.b=theta.b,W=W,ymax=ymax,dy=dy,nx=nx,
            groupfromy=groupfromy)
  fit$par=invtfm(fit$par,FUN)
  return(fit)
}


#' @title Fit detection hazard with a Markov or hidden Markov availability model.
#'
#' @description
#' \code{fit.hmltm} Estimates detection probability from (1) line transect data that includes forward 
#' detection distances and (2) estimated Markov model or hidden Markov model availability prameters. 
#'
#' @param xy data frame with one line per detection containing perpendicular distance ($x) and
#' forward distance ($y) and any covariates in the model.
#' @param pars starting parameter values.
#' @param FUN detection hazard functional form name (character)
#' @param models list of characters with elements \code{$y} and \code{$x} specifying y- and 
#' x-covariate models. Either \code{NULL} or regression model format (without response on left).
#' @param survey.pars list.
#' @param hmm.pars HMM parameters of animals' availability data.
#' @param control.fit list with elements:
#'              \code{$hessian (logical)} - if TRUE Hessian is estimated and returned, else not;
#'              \code{$nx (scalar)} - determines the number of intervals to use with Simpson's rule 
#'              integration over y. \code{nx=64} seems safe; smaller number makes computing faster.
#' @param control.optim as required by \code{\link{optim}}.
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
#'
#' @return A list comprising:
#' \itemize{
#'  \item{xy} {dat used in fitting (input reflection).}
#'  \item{phats} {estimated detection probabilities of all detections.}
#'  \item{phat} {1/mean(1/phats).}
#'  \item{pzero} {estimated detection probabilities at perpendicular distance zero for 
#'  detected population.}
#'  \item{h.fun} {=FUN (input reflection).}
#'  \item{models} {=models (input reflection).}
#'  \item{fit} {output from \code{\link{fit.xy}}.}
#'  \item{Loglik} {log-likelihood function at MLE.}
#'  \item{AIC} {AIC.}
#'  \item{x} {vector of x-values for plotting perpendicular distance fit.}
#'  \item{p} {vector of detection function values for plotting perpendicular distance fit.}
#'  \item{fitpars} {a list containing all the given parameters controlling the fit (survey.pars,hmm.pars,
#'  control.fit,control.optim).}
#' }
#'  
#' @references Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 
#' 2013. Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
fit.hmltm=function(xy,pars,FUN,models=list(y=NULL,x=NULL),survey.pars,hmm.pars,
                   control.fit,control.optim,groupfromy=NULL)
{
  # unpack things:
  hessian=control.fit$hessian
  nx=control.fit$nx
  theta.f=survey.pars$theta.f
  theta.b=survey.pars$theta.b
  W=survey.pars$W
  ymax=survey.pars$ymax
  dy=survey.pars$dy
  pm=hmm.pars$pm
  Pi=hmm.pars$Pi
  delta=hmm.pars$delta
  
  # fit model
  est=fit.xy(pars=pars,xy=xy,FUN=FUN,models=models,pm=pm,Pi=Pi,delta=delta,theta.f=theta.f,theta.b=theta.b,W=W,ymax=ymax,dy=dy,nx=nx,hessian=hessian,
             control=control.optim,groupfromy=groupfromy)
  
  # store transformed parameters
  b=tfm(est$par,FUN)
  est$b=b
  
  n=length(xy$x)
  xs=seq(0,W,length=nx)
  if(is.nullmodel(models)) {
    px=rep(0,nx)
    px=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=rep(b,nx),pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    phat=sintegral(px,xs)/W
  } else {
    px=matrix(rep(0,nx*n),nrow=n)
    phat=parea=rep(0,n)
    covb=make.covb(b,FUN,models,xy) # put covariates into paramerters
    nb=length(covb)/n
    for(i in 1:n) {
      start=(i-1)*nb+1
      bi=c(rep(covb[start:(start+nb-1)],nx)) # nx replicates of covb for ith detection
      px[i,]=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    }
    phat=apply(px,1,sintegral,xs)/W
  }
  
  # log-likelihood and AIC
  llik=-est$value
  npar=length(est$par)
  aic=-2*llik+2*npar
  
  # return it all:
  if(is.nullmodel(models)) pzero=px[1]
  else pzero=mean(px[,1])
  fitpars=list(survey.pars=survey.pars,hmm.pars=hmm.pars,control.fit=control.fit,control.optim=control.optim)
  hmmlt=list(xy=xy,phats=phat,phat=1/mean(1/phat),pzero=pzero,h.fun=FUN,models=models,fit=est,Loglik=llik,AIC=aic,x=xs,p=px,fitpars=fitpars)
  class(hmmlt)="hmmlt"
  return(hmmlt)
  
}


#' @title HMM line transect model likelihood with perpendicular and forward distances.
#'
#' @description
#' Evaluates HMM line transect model likelihood with both perp. and forward distance data.
#'  
#' @param b likelihood function parameters.
#' @param xy detections data, including x and y values for each detection, plus any covariates 
#' in \code{model}).
#' @param FUN detection hazard function name.
#' @param models model specification, as in \code{\link{est.hmltm}}.
#' @param pm state-dependent Bernoulli parameters.
#' @param Pi Markov model transition probability matrix.
#' @param delta Markov model stationary distribution.
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param W perpendicular truncation distance for fitting.
#' @param ymax forward distance by which detection hazard has fallen to zero.
#' @param dy Markov model distance step size.
#' @param nx number of x-values to use in evaluating detection function.
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
negllik.xandy=function(b,xy,FUN,models=list(y=NULL,x=NULL),pm,Pi,delta,theta.f,theta.b,W,ymax,dy,nx=100,groupfromy=NULL)
{
  # If asked to group close y's, split data accordingly:
  if(!is.null(groupfromy)){
    if(groupfromy<0) stop("groupfromy must be greater than zero.")
    if(groupfromy>max(xy$y)) stop(paste("With groupfromy=",groupfromy,"you are grouping all forward distances: can't do this."))
    grouped=xy$y<=groupfromy
    xgrp=xy$x[grouped];ygrp=xy$y[grouped]
    x=xy$x[!grouped];y=xy$y[!grouped]
    covb=make.covb(b,FUN,models,xy[!grouped,]) # put covariates into paramerters; ungrouped data
    covbgrp=make.covb(b,FUN,models,xy[grouped,]) # put covariates into paramerters; ungrouped data
  } else {
    x=xy$x;y=xy$y
    covb=make.covb(b,FUN,models,xy) # put covariates into paramerters
  }
  
  # Deal with the ungrouped data:
  # ----------------------------
  n=length(x)
  if(length(y)!=n) stop("Length of y: ",length(y)," not equal to length of x: ",n,"\n")
  llik=0
  # calculate p(see first @ y|x) at each (x,y):
  li=p.xy(x,y,hfun=FUN,b=covb,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b)
  if(any(li<.Machine$double.xmin)) return(.Machine$double.xmax)
  xs=seq(0,W,length=nx)
  nb=length(covb)/n # number of parameters for each observation
  if(is.nullmodel(models)) {
    bi=c(rep(covb[1:nb],nx)) # nx replicates of covb for ith detection
    ps=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    p=sintegral(ps,xs)/W
  } else {
    ps=matrix(rep(0,n*nx),nrow=n)
    for(i in 1:n) {
      start=(i-1)*nb+1
      bi=c(rep(covb[start:(start+nb-1)],nx)) # nx replicates of covb for ith detection
      ps[i,]=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    }
    p=apply(ps,1,sintegral,xs)/W
  }
  
  llik=sum(log(li)-log(p))
  
  
  # Deal with the grouped data:
  # ----------------------------
  if(!is.null(groupfromy)){
    ngrp=length(xgrp)
    if(ngrp>0){
      if(length(ygrp)!=ngrp) stop("Length of ygrp: ",length(ygrp)," not equal to length of xgrp: ",ngrp,"\n")
      llik.grp=0
      xs=seq(0,W,length=nx)
      nb=length(covbgrp)/ngrp # number of parameters for each observation
      if(is.nullmodel(models)) {
        bi=c(rep(covbgrp[1:nb],nx)) # nx replicates of covbgrp for ith detection
        ps=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
        p=sintegral(ps,xs)/W
      } else {
        ps=matrix(rep(0,ngrp*nx),nrow=ngrp)
        for(i in 1:ngrp) {
          start=(i-1)*nb+1
          bi=c(rep(covbgrp[start:(start+nb-1)],nx)) # nx replicates of covbgrp for ith detection
          ps[i,]=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
        }
        p=apply(ps,1,sintegral,xs)/W
      }
      #  # calculate p(see first at or after groupfromy|x) at each x:
      pfar=p.xy(x=xgrp,y=rep(groupfromy,ngrp),hfun=FUN,b=covbgrp,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,cdf=TRUE)
      p0=p.xy(x=xgrp,y=rep(0,ngrp),hfun=FUN,b=covbgrp,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,cdf=TRUE)
      pnear=p0-pfar
      llik=llik+sum(log(pnear)-log(p))
    }
  }
  
  return(-llik)
}


#' @title HMM line transect model likelihood with only perpendicular distances.
#'
#' @description
#' Evaluates HMM line transect model likelihood with both perpendicular distance but no forward distance data.
#'  
#' @param b likelihood function parameters.
#' @param xy detections data, including x values for each detection, plus any covariates 
#' in \code{model}).
#' @param FUN detection hazard function name.
#' @param models model specification, as in \code{\link{est.hmltm}}.
#' @param pm state-dependent Bernoulli parameters.
#' @param Pi Markov model transition probability matrix.
#' @param delta Markov model stationary distribution.
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param W perpendicular truncation distance for fitting.
#' @param ymax forward distance by which detection hazard has fallen to zero.
#' @param dy Markov model distance step size.
#' @param nx number of x-values to use in evaluating detection function.
#' 
negllik.x=function(b,xy,FUN,models,pm,Pi,delta,theta.f,theta.b,W,ymax,dy,nx=100)
{
  covb=make.covb(b,FUN,models,xy) # put covariates into paramerters
  x=xy$x;y=xy$y
  n=length(x)
  if(length(y)!=n) stop("Length of y: ",length(y)," not equal to length of x: ",n,"\n")
  llik=0
  # calculate p(see BY y=0 |x) at each x:
  li=p.xy(x,y=rep(0,n),hfun=FUN,b=covb,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
  # calculate p(see by min y in window|x) over all x:
  #*# xs=seq(0,W,length=nx)
  #*# px=p.xy(x=xs,y=xy$y,hfun=FUN,b=b,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
  #*# wt=c(0.5,rep(1,(nx-2)),0.5)*W/nx
  #*# p=sum(px*wt) # effective strip width
  # cat("pars=",invtfm(b,FUN),"\n")
  xs=seq(0,W,length=nx)
  nb=length(covb)/n # number of parameters for each observation
  if(is.nullmodel(models)) {
    bi=c(rep(covb[1:nb],nx)) # nx replicates of covb for ith detection
    ps=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    p=sintegral(ps,xs)/W
  } else {
    ps=matrix(rep(0,n*nx),nrow=n)
    for(i in 1:n) {
      start=(i-1)*nb+1
      bi=c(rep(covb[start:(start+nb-1)],nx)) # nx replicates of covb for ith detection
      ps[i,]=p.xy(x=xs,y=rep(0,nx),hfun=FUN,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=TRUE)
    }
    p=apply(ps,1,sintegral,xs)/W
  }
  
  if(any(li<.Machine$double.xmin)) return(.Machine$double.xmax)
  else{
    llik=sum(log(li)-log(p))
    return(-llik)
  }
  return(-llik)
}


#' @title HMM line transect model likelihood with at least some forward distances.
#'
#' @description
#' Evaluates HMM line transect model likelihood with perpendicular and forward distance data.
#'  
#' @param b likelihood function parameters.
#' @param xy detections data, including x values for each detection and y values for some, 
#' plus any covariates in \code{model}).
#' @param FUN detection hazard function name.
#' @param models model specification, as in \code{\link{est.hmltm}}.
#' @param pm state-dependent Bernoulli parameters.
#' @param Pi Markov model transition probability matrix.
#' @param delta Markov model stationary distribution.
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param W perpendicular truncation distance for fitting.
#' @param ymax forward distance by which detection hazard has fallen to zero.
#' @param dy Markov model distance step size.
#' @param nx number of x-values to use in evaluating detection function.
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
#' 
negllik.xy=function(b,xy,FUN,models=list(y=NULL,x=NULL),pm,Pi,delta,theta.f,theta.b,W,ymax,dy,nx=100,groupfromy=NULL)
{
  xydat=xy[!is.na(xy$y),] # detections with x (perp) and y (forward) data
  xdat=xy[is.na(xy$y),]   # detections with x (perp) data only (no y data)
  xy.negllik=negllik.xandy(b,xydat,FUN,models,pm,Pi,delta,theta.f,theta.b,W,ymax,dy,nx,groupfromy)
  x.negllik=0
  if(length(xdat$x)>0) negllik.x(b,xdat,FUN,models,pm,Pi,delta,theta.f,theta.b,W,ymax,dy,nx)
  negllik=xy.negllik+x.negllik
  return(negllik)
}



#' @title Calculates probability of first observing animal at (x,y).
#'
#' @description
#' Just calls \code{\link{p.xy1}} >= once, cycling through 3rd index of Pi and 2nd indices of pm and delta.
# 
#' @param x perpendicular distance.
#' @param y forward distance.
#' @param hfun detection hazard function name.
#' @param b detection hazard function parameters.
#' @param pm state-dependent Bernoulli parameters.
#' @param Pi Markov model transition probability matrix.
#' @param delta Markov model stationary distribution.
#' @param ymax forward distance at or beyond which things are undetectable.
#' @param dy distance step in forward distance direction corresponding to the time step size of the 
#' Markov model governed by \code{Pi}.
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param ally If TRUE calculates probability of detection at all (i.e., at some y. (Argument
#' \code{cdf} is redundant in this case).
#' @param cdf If TRUE returns 1 -(the cdf at y), which is equal to the probability of detection BY y.
#' This differes from specifying ally=TRUE in that ally=TRUE calculates the cdf from ymax 
#' to y=0, whereas cdf=TRUE calculates the cdf from ymax to y.
#' 
#' @seealso \code{\link{p.xy1}}
#' 
p.xy=function(x,y,hfun,b,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=FALSE)
{
  if(is.vector(pm)&!is.matrix(Pi) | !is.vector(pm)&is.matrix(Pi)) stop("Single animal: pm is not a vector or Pi is not a matrix")
  if(is.vector(pm)) { # convert to matrix and array so can use loop below
    #    Pi=array(Pi,dim=c(2,2,1))
    dimPi=dim(Pi)
    Pi=array(Pi,dim=c(dimPi[1:2],1))
    pm=matrix(pm,ncol=1)
    delta=matrix(delta,ncol=1)
  }
  nw=dim(pm)[2]
  p=rep(0,length(x))
  if(dim(Pi)[3]!=nw) stop("Last dimensions of pm and Pi must be the same.")
  for(w in 1:nw) {
    Pi.w=Pi[,,w]
    pm.w=pm[,w]
    delta.w=delta[,w]
    p.w=p.xy1(x,y,hfun,b,pm.w,Pi.w,delta.w,ymax,dy,theta.f,theta.b,ally,cdf)
    ##print(c('p.w',p.w,' p',p))
    p=p+p.w
  }
  ##print(p/nw);
  ##print("step");
  return(p/nw)
}



#---------------------- End estimation functions for x and y ----------------------
#                       -----------------------------------
