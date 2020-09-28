#                          -------------------------------
#------------------------- Start Goodness of Fit functions -----------------------


#' @title Goodness-of-fit in forward dimension.
#'
#' @description
#' Calculates goodness-of-fit in forward dimension, plots fit, and returns p-value and other stuff.
#' Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs), and Cramer-von Mises
#' p-value (which is based on all cdf values).
# 
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param ks.plot If TRUE, does CDF-EDF plot (similar to a Q-Q plot, but using cumulative distribution function 
#' values rather than quantiles). Point corresponding to largest difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
#' @param seplots if TRUE does additional diagnostic plots
#' @param smult multiplier to size circles in third plot.
#' @param ymax forward distance at which detection probability is assumed to be zero. 
#' @param breaks breaks for Chi-squared goodness of fit test.
hmmlt.gof.y=function(hmltm,ks.plot=TRUE,seplots=FALSE,smult=5,ymax=hmmlt$fitpars$survey.pars$ymax,breaks=NULL)
{
  hmmlt=hmltm$hmltm.fit
  dat=hmmlt$xy[!is.na(hmmlt$xy$y),]
  hfun=hmmlt$h.fun
  b=hmmlt$fit$b
  theta.b=hmmlt$fitpars$survey.pars$theta.b
  theta.f=hmmlt$fitpars$survey.pars$theta.f
  
  if(is.null(dat$y)) stop("Can't do goodness-of-fit in y-dimension when don't have y observations!")
  models=hmmlt$models
  ###  ymax=hmmlt$fitpars$survey.pars$ymax
  dy=hmmlt$fitpars$survey.pars$dy
  Pi=hmmlt$fitpars$hmm.pars$Pi
  pm=hmmlt$fitpars$hmm.pars$pm
  delta=hmmlt$fitpars$hmm.pars$delta
  n=length(dat$x)
  covb=make.covb(b,hfun,models,dat) # put covariates into paramerters
  Fy=p.xy(dat$x,dat$y,hfun,b=covb,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=TRUE)
  F0=p.xy(dat$x,dat$y,hfun,b=covb,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=TRUE,cdf=FALSE)
  #  Fy=Fx.cox.Q(dat$x,dat$y,mu,ystart,Q,b) # area up to y
  #  F0=Fx.cox.Q(dat$x,rep(0,n),mu,ystart,Q,b) # area up to y=0
  Fy0=Fy/F0
  Fy0.order=order(Fy0)
  yy=dat$y[Fy0.order]
  xx=dat$x[Fy0.order]
  cdf=Fy0[Fy0.order]
  e.cdf=order(cdf)/n
  # K-S statistic
  dF=cdf-e.cdf
  worst=which(abs(dF)==max(abs(dF))) # mark point on which Kolmogarov test hinges
  Dn=max(abs(dF))*sqrt(n)
  p=p.kolomogarov(Dn)
  p.cvm=cvm.test(Fy0)$p.value
  # plots
  if(ks.plot) {
    plot(1-e.cdf,cdf,xlab="Empirical Distribution Function",ylab="Cumulative Distribution Function",main="Forward Dist. CDF-EDF Plot",xlim=c(0,1),ylim=c(0,1),pch="+")
    lines(c(0,1),c(1,0))
    points(1-e.cdf[worst],cdf[worst],col="red") # mark point on which Kolmogarov test hinges
    if(seplots){
      plot(xx,dF,xlab="Perpendicular distance",ylab="CDF-Empirical CDF")
      lines(c(0,max(xx)),c(0,0),lty=2)
      size=smult*dF/max(dF)
      dFcol=c("red","black","black")
      plot(dat$x[Fy0.order],dat$y[Fy0.order],xlab="Perpendicular distance",ylab="Forward distance",cex=abs(size),col=dFcol[sign(dF)+2],main="CDF-Empirical CDF (red=negative)")
    }
  }
  
  # Now calulate Chi-squared gof
  if(is.null(breaks)) breaks =seq(0,ymax,length=11)
  chisq.y = chisq.gof.y(hmltm,breaks,nys=250)
  
  return(list(p.ks=1-p,p.cvm=p.cvm,qq.x=e.cdf,qq.y=cdf,y=yy,p.chisq=chisq.y$p.chisq,chisq.gof=chisq.y))
}

#' @title Chi-squared goodness-of-fit in forward dimension.
#'
#' @description
#' Calculates Chi-squared goodness-of-fit in forward dimension, plots fit, and returns p-value and other stuff.
#' Returns a data frame  (\code{goftable}) with columns for start and end of bins, observed and expected
#' frequencies in the bins, and the Chi-squared statistic for the bin. Also 
#' returns the Chi-squared statistic(\code{X2stat}), the degrees of freedom of the 
#' test (\code{df}), and the associated p-value (\code{p.Chisq}).
# 
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param breaks cutpoints on forward-distance axis defining bins
#' @param nys number of y-values at which to calculate function before using 
#' \code{approxfun} approximate it at arbitrary resolution.
#' 
#' @details 
#' Uses \code{approxfun} to approximate the forward distance pdf and then uses
#' \code{integrate} to integrate it within each bin.
#'   
chisq.gof.y = function(fit,breaks,nys=250) {
  fplot = fyfit.plot(fit,breaks=breaks,allx=FALSE,nys=nys,doplot=FALSE)
  model.df = length(fit$hmltm.fit$fit$par)
  n = length(!is.na(fplot$xy$object))
  xcol = which(names(fplot)=="x")
  ycol = which(names(fplot)=="y")
  if(length(xcol)==0 & length(ycol)==0) stop("Must have column named x or named y")
  if(length(xcol)!=0 & length(ycol)!=0) stop("Can't have column named x AND named y")
  if(length(xcol)!=0) {
    x = fplot[[xcol]]
  } else {
    x = fplot[[ycol]]
  }
  f = approxfun(x,fplot$f)
  nbins = length(fplot$breaks)-1
  o = e = rep(NA,nbins)
  for(i in 1:nbins) {
    o[i] = fplot$barheight[i]* (fplot$breaks[i+1]-fplot$breaks[i]) 
    e[i] = integrate(f,fplot$breaks[i],fplot$breaks[i+1])$value
  }
  o = n * o/sum(o)
  p = e/sum(e)
  e = n * p
  d2 = (o-e)^2
  X2 = d2/e
  X2stat = sum(X2)
  nobs = length(o)
  df = nobs-model.df
  pval = 1 - pchisq(X2stat,df)
  goftable = data.frame(From=breaks[-nobs],To=breaks[-1],Observed=o, Expected=e,X2=X2)
  if(any(e<5)) warning("Some expected values < 5; p-value may not be reliable.")
  return(list(goftable,X2stat=X2stat,df=df,p.chisq=pval))
}




#' @title Kolmogarov-Smirnov goodness-of-fit p-value.
#'
#' @description
#' Kolmogarov-Smirnov goodness-of-fit p-value calculation.
# 
#' @param x value of Kolmogarov-distributed random variable at which to evaluate.
#' @param inf approximation to infinity (a large number).
#' @param dp approximation convergence criterion.
#' 
#' @details
#' Calculates p-value for Kolmogarov distribution at x, approximating infite sum
#' \code{sqrt(2*pi)/x*sum{i=1}^infty exp(-(2i-1)^2*pi^2/(8x^2)))}
#' by a finite sum to inf (default 1000) if sum to inf and inf+1 differ by less
#' than dp (default 1e-4), else sum until difference is less than dp.
p.kolomogarov=function(x,inf=1000,dp=1e-4)
{
  infsum=rep(0,inf)
  i=1:inf
  K=sqrt(2*pi)/x
  p=p1=K*sum(exp(-(2*i-1)^2*pi^2/(8*x^2)))
  dp=1
  while(dp>1e-4) {
    inf=inf+1
    p=p1+K*exp(-(2*inf-1)^2*pi^2/(8*x^2))
    dp=abs(p-p1)
  }
  return(p=p)
}


#' @title Goodness-of-fit in perpendicular dimension.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular dimension, plots fit, and returns p-value and 
#' other stuff. Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs), and Cramer-von Mises
#' p-value (which is based on all cdf values).
# 
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param ks.plot If TRUE, does CDF-EDF plot (similar to a Q-Q plot, but using cumulative distribution function 
#' values rather than quantiles). Point corresponding to largest difference between#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
hmmlt.gof.x=function(hmltm,ks.plot=TRUE){
  hmmlt=hmltm$hmltm.fit
  n=length(hmmlt$xy$x)
  edf=(1:n)/n
  cdf=fitted.esw(hmmlt,to.x=TRUE,all=TRUE)/fitted.esw(hmmlt,all=TRUE)
  cdf.order=order(cdf)
  cdf=cdf[cdf.order]
  e.cdf=cdf.order/n
  # K-S statistic
  dF=cdf-edf
  worst=which(abs(dF)==max(abs(dF)))
  Dn=max(abs(dF))*sqrt(n)
  p.ks=p.kolomogarov(Dn)
  p.cvm=cvm.test(cdf)$p.value # Under model, cdf values are from uniform; default for cvm.test is "punif"
  if(ks.plot) {
    plot(edf,cdf,pch="+",xlim=c(0,1),ylim=c(0,1),xlab="Empirical Distribution Function",ylab="Cumulative Distribution Function",main="Perp. Dist. CDF-EDF Plot")
    lines(c(0,1),c(0,1))
    points(edf[worst],cdf[worst],col="red")
  }
  
  return(list(p.ks=1-p.ks,p.cvm=p.cvm,qq.x=edf,qq.y=cdf,x=hmmlt$xy$x[cdf.order]))
}




#-------------------------- End Goodness of Fit functions ------------------------
#                          ------------------------------

