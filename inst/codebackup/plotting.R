#                        ------------------------
#----------------------- Start Plotting functions --------------------------------------

#' @title Detection hazard plotting.
#'
#' @description
#'  Plots detection hazard contours.
#'  
#' @param hfun detection hazard function name
#' @param pars detection hazard function parameter vector
#' @param dat data frame.
#' @param models model list, as for \code{\link{est.hmltm}} for example.
#' @param xrange range of x-axis.
#' @param yrange range of y-axis.
#' @param nx number of points on x-axis at which to evaluate detection hazard function.
#' @param ny number of points on y-axis at which to evaluate detection hazard function.
#' @param type "contour", "persp", "image" or "both" (for image and contour).
#' @param nlevels number of contour levels.
#' @param add if TRUE adds to existing plot, else creates new plot.
#' @param col colour of plot.
#' @param logscale If TRUE, plots hazard values on log scale.
#' @param xlab x label.
#' @param ylab y label.
#' @param theta argument for \code{\link{image}}.
#' @param phi argument for \code{\link{image}}.
#' @param ... other arguments to image, contour or persp.
h.plot=function(hfun,pars,dat=NULL, models=NULL, xrange=c(0,50),yrange=xrange,nx=50,ny=nx,
                type="contour",nlevels=20,add=FALSE,col="black",logscale=FALSE,
                xlab="Perpendicular distance",ylab="Forward distance",theta=90,phi=35,...)
{
  if(type!="persp" & type!="both" & type!="contour" & type!="image") stop("Agrument `type' must be `persp', `contour', `image'or `both' (for contour on image).")
  b=tfm(pars,hfun)
  nb=length(b)
  n=dim(dat)[1]
  if(n>1) {
    warning("Only allowed to plot at one level of dat. Only using first row of dat")
    dat=dat[1,]
    n=1
  }
  h=match.fun(hfun)
  if(!is.null(dat)) {
    covb=make.covb(b,hfun,models,dat)
    nb=length(covb)/n
  }
  else covb=b
  x=seq(xrange[1],xrange[2],length=nx)
  y=seq(yrange[1],yrange[2],length=ny)
  for(i in 1:n){
    start=(i-1)*nb+1
    bi=covb[start:(start+nb-1)]
    h.xy=outer(x,y,FUN=hfun,b=bi)
    if(logscale) h.xy=log(h.xy)
    if(i>1) add=TRUE
    if(type=="image" | type=="both") image(x,y,h.xy,add=add,xlab=xlab,ylab=ylab)
    if(type=="both") add=TRUE
    if(type=="contour" | type=="both") contour(x,y,h.xy,add=add,nlevels=nlevels,xlab=xlab,ylab=ylab,col=col)
    if(type=="persp") persp(x,y,h.xy,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab="Detection hazard")
  }
}

#' @title detection function plotting.
#'
#' @description
#' Plots the detection function (1- the survivor function) of a fitted model in two dimensions.
#'  
#' @param hmltm output from \code{\link{est.hmltm}}
#' @param obs row (observation) in data to use in plotting
#' @param new.ymax a forward distance beyond which things can't be detected (used to override the 
#' ymax already in \code{hmltm}).
#' @param new.pars model parameter vector  (used to override the pars already in \code{hmltm}).
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param xrange range of x-axis.
#' @param yrange range of y-axis.
#' @param nx number of points on x-axis at which to evaluate detection function.
#' @param ny number of points on y-axis at which to evaluate detection function.
#' @param xlab x label.
#' @param ylab y label.
#' @param type "contour", "persp", "image" or "both" (for image and contour).
#' @param nlevels number of contour levels.
#' @param add If TRUE, adds to current plot, else makes new plot.
#' @param col colour of plot.
#' @param logscale If TRUE, plots hazard values on log scale.
#' @param theta argument for \code{\link{image}}.
#' @param phi argument for \code{\link{image}}.
#' @param zlab label in z-dimension for persp plots
#' @param values If TRUE, returns hazard function values on nx by ny grid.
#' @param ... other arguments to \code{\link{image}}, \code{\link{contour}} or \code{\link{persp}}.
f.plot=function(hmltm,obs=1:length(hmltm$hmltm.fit$xy$x),new.ymax=NULL,new.pars=NULL,
                theta.f=0,theta.b=90,
                xrange=c(0,max(hmltm$hmltm.fit$xy$x)),yrange=c(0,1.5*max(na.omit(hmltm$hmltm.fit$xy$y))),
                nx=50,ny=nx,
                xlab="Perpendicular distance",ylab="Forward distance",type="contour",
                nlevels=20,add=FALSE,col="black",logscale=FALSE,theta=90,phi=35,
                zlab="pdf of first detections",values=FALSE,...)
{
  hfun=hmltm$hmltm.fit$h.fun
  pars=hmltm$hmltm.fit$fit$par
  dat=hmltm$hmltm.fit$xy
  models=hmltm$hmltm.fit$models
  Pi=hmltm$hmltm.fit$fitpars$hmm.pars$Pi
  pm=hmltm$hmltm.fit$fitpars$hmm.pars$pm
  delta=hmltm$hmltm.fit$fitpars$hmm.pars$delta
  ymax=hmltm$hmltm.fit$fitpars$survey.pars$ymax
  dy=hmltm$hmltm.fit$fitpars$survey.pars$dy
  xmax=hmltm$hmltm.fit$fitpars$survey.pars$W
  if(!is.null(new.pars)) pars=new.pars
  if(!is.null(new.ymax)) ymax=new.ymax
  
  if(type!="persp" & type!="both" & type!="contour" & type!="image") stop("Agrument `type' must be `persp', `contour', `image'or `both' (for contour on image).")
  #  if(max(xrange)>ymax) {
  #    ymax=max(xrange)+dy
  #    warning("ymax<max(xrange) so ymax set equal to max(xrange)+dy")
  #  }
  if(max(yrange)>ymax) {
    ymax=max(yrange)+dy
    warning("ymax<max(yrange) so ymax set equal to max(yrange)+dy")
  }
  #  h=match.fun(hfun)
  b=tfm(pars,hfun)
  nb=length(b)
  n=dim(dat[obs,])[1]
  if(n>1) {
    warning("Only allowed to plot at one level of covariates. Only using first row of dat[obs,]")
    dat=dat[1,]
    n=1
  }
  if(!is.null(dat)) {
    covb=make.covb(b,hfun,models,dat)
    nb=length(covb)/n
  }
  else covb=b
  x=seq(xrange[1],xrange[2],length=nx)
  y=seq(yrange[1],yrange[2],length=ny)
  f.xy=matrix(rep(0,nx*ny),nrow=nx)
  for(i in 1:nx){
    f.xy[i,]=p.xy(rep(x[i],ny),y,hfun=hfun,b=rep(covb,ny),pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b)
  }
  if(logscale) f.xy=log(f.xy)
  if(type=="image" | type=="both") image(x,y,f.xy,add=add,xlab=xlab,ylab=ylab,...)
  if(type=="both") add=TRUE
  if(type=="contour" | type=="both") contour(x,y,f.xy,add=add,nlevels=nlevels,xlab=xlab,ylab=ylab,col=col,...)
  if(type=="persp") persp(x,y,f.xy,theta=theta,phi=phi,xlab=xlab,ylab=ylab,zlab=zlab,...)
  if(values) return(list(x=x,y=y,z=f.xy))
}


#' @title Plotting pdf of forward detection distances.
#'
#' @description
#' Plots pdf of forward detection distances (f(y)) averaged over detectoins on top of forward detection 
#' distance distribution. 
#'  
#' @param hmltm object produced by \code{\link{est.hmltm}}.
#' @param values If TRUE, returns function values used in plot.
#' @param breaks histogram break points.
#' @param allx If TRUE, returns matrix with p(see at y, given x) for each x in data in hmmlt (dashed gray lines).
#' @param nys number of y-values to use in plotting, spaced equally from 0 to 
#' \code{hmmlt$fitsurvey.pars$ymax}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main title.
#' @param ylim y-value limits.
#' @param text.cex relative text size.
#' @param theta.b REDUNDANT must = 90.
#' @param doplot If FALSE, does not produce a plot, else does.
#' @param trackprogress If TRUE, prints line for each of the x-values for which
#' the function is computed, as it is computed.
#' 
fyfit.plot=function(hmltm,values=TRUE,breaks=NULL,allx=FALSE,nys=250,
                    xlab="Forward distance (y)",ylab="pdf(y)",main="",ylim=NULL,theta.b=90,
                    text.cex=0.66,doplot=TRUE,trackprogress=FALSE)
{
  hmmlt=hmltm$hmltm.fit
  if(theta.b !=90) stop("This function only allows theta.b=90 degrees at present.")
  # unpack things from hmmlt object:
  pm=hmmlt$fitpars$hmm.pars$pm
  Pi=hmmlt$fitpars$hmm.pars$Pi
  delta=hmmlt$fitpars$hmm.pars$delta
  theta.f=hmmlt$fitpars$survey.pars$theta.f
  theta.b=hmmlt$fitpars$survey.pars$theta.b
  W=hmmlt$fitpars$survey.pars$W
  ymax=hmmlt$fitpars$survey.pars$ymax
  dy=hmmlt$fitpars$survey.pars$dy
  hfun=hmmlt$h.fun
  b=tfm(hmmlt$fit$par,hfun)
  models=hmmlt$models
  xy=hmmlt$xy[!is.na(hmmlt$xy$y),]
  if(is.null(breaks)) {
    breaks=seq(0,max(xy$y),nclass=10)
  } 
  else if(breaks[length(breaks)]<max(xy$y)){
    warning("Found max(breaks) less than max y observed.")
  }
  breakmax=max(breaks)
  
  #  x=sort(xy$x)
  x=xy$x
  nx=length(x)
  covb=make.covb(b,hfun,models,xy) # put covariates into parameters
  nb=length(covb)/nx # number of parameters for each observation
  
  maxy=0;maxi=0
  for(i in 1:nx) { # find maximum y
    maxyi=max(gety.obs(min(x),TRUE,0,theta.f,theta.b,ymax,dy))
    #CJ#    maxyi=max(gety.obs(min(x),0,theta.f,theta.b,ymax,dy))
    if(maxy<maxyi) {
      maxy=maxyi
      maxi=i
    }
  }
  ys=gety.obs(x[maxi],TRUE,0,theta.f,theta.b,ymax,dy) # get y's spanning maximum forward distance
  #CJ#  ys=gety.obs(x[maxi],0,theta.f,theta.b,ymax,dy) # get y's spanning maximum forward distance
  fyx=matrix(rep(NA,nx*length(ys)),ncol=nx) # set matrix to largest y-dimension
  
  for(i in nx:1){
    yi=gety.obs(x[i],TRUE,0,theta.f,theta.b,ymax,dy) # TRUE,0 used to be NULL. get all y's in view at this x
    #CJ#    yi=gety.obs(x[i],0,theta.f,theta.b,ymax,dy) # TRUE,0 used to be NULL. get all y's in view at this x
    #    ndys=length(yi)
    #    if(!is.null(nys)&length(yi)>nys) {
    #      ndys=ceiling(length(yi)/nys)
    #      if(ndys>1) yi=seq(min(yi),ndys*nys*dy,length=nys)
    #    } else nys=length(yi)
    nyi=length(yi)
    start=(i-1)*nb+1
    bi=c(rep(covb[start:(start+nb-1)],nyi)) # nx replicates of covb for ith detection
    fyx[1:nyi,i]=p.xy(x=rep(x[i],nyi),y=yi,hfun=hfun,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b)
    if(trackprogress) cat("Done ",nx-i+1," of ",nx,"\n")
  }
  f=apply(na.omit(fyx),1,mean)
  #  wt=c(0.5,rep(1,(dim(fyx)[1]-2)),0.5)
  #  farea=sum(f*wt*dy)
  farea=sintegral(f,ys)
  #  fy=approxfun(ys,f)
  #  farea=integrate(fy,lower=min(y),upper=max(y))$value
  #  if(is.null(breaks)) {hst=hist(xy$y,plot=FALSE);breaks=hst$breaks} 
  #  else {
  #    if(max(breaks)<max(y)) {
  #      breaks=c(breaks,max(y))
  #      warning("Set max(breaks)=maxy(y)=",max(y)," for plotting")
  #    }
  hst=hist(xy$y[xy$y<=breakmax],breaks=breaks,plot=FALSE)
  #  }
  barheight=hst$density*farea
  if(is.null(ylim)){
    if(allx) ylim=c(0,max(fyx,barheight))
    else ylim=c(0,max(f,barheight))  
  }
  if(doplot) {
    histline(barheight,breaks=breaks,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
    if(allx) for(i in 1:dim(fyx)[2]) lines(ys,fyx[,i],col="gray",lty=2)
    lines(ys,f,lwd=2)
    nb=length(hst$breaks)
    text(hst$breaks[1:(nb-1)]+diff(hst$breaks)/2,0,hst$count,cex=text.cex)    
  }
  
  if(values) return(list(xy=xy,y=ys,fyx=fyx,f=f,breaks=breaks,barheight=barheight))
}


#' @title Detection function fit plotting.
#'
#' @description
#' Plots detection function on top of forward detection distance distribution. Returns a bunch of stuff.
#'  
#' @param hmltm object produced by \code{\link{est.hmltm}}.
#' @param allx If TRUE, returns matrix with p(see at x, given covars) for each x in data in hmmlt (dashed gray 
#' lines).
#' @param values If TRUE, returns function values used in plot.
#' @param breaks histogram break points.
#' @param ylim y-value limits.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param type "prob" for probability function, "density" for pdf.
#' @param text.cex relative text size.
#' 
fxfit.plot=function(hmltm,allx=FALSE,values=TRUE,breaks=NULL,ylim=NULL,xlab=NULL,ylab=NULL,
                    type="prob",text.cex=0.66)
{
  hmmlt=hmltm$hmltm.fit
  if(type!="density" & type!="prob") stop("Invalid type: must be `density' or `prob'.")
  nx=length(hmmlt$x)
  if(is.null(breaks)) h=hist(hmmlt$xy$x,plot=FALSE)
  else h=hist(hmmlt$xy$x,plot=FALSE,breaks=breaks)
  h.area=sum(diff(h$breaks)*h$density)
  #  esw=((apply(hmmlt$p,1,sum)-((hmmlt$p[,1]+hmmlt$p[,nx])/2))/nx)*max(hmmlt$x)
  if(is.nullmodel(hmmlt$models)) {
    esw=sintegral(hmmlt$p,hmmlt$x)
    mean.esw=esw
    mean.p=hmmlt$p
    scaled.density=h$density*mean.esw/h.area
    E.density=hmmlt$p*h.area/mean.esw
    mean.E.density=E.density
  } else {
    n=dim(hmmlt$p)[1]
    mean.p=rep(0,dim(hmmlt$p)[2])
    esw=rep(0,n)
    E.density=matrix(rep(esw,nx),nrow=n)
    for(i in 1:n) esw[i]=sintegral(hmmlt$p[i,],hmmlt$x)
    wt=(1/esw)/sum(1/esw)
    mean.esw=mean(esw)
    for(i in 1:n) mean.p=mean.p+hmmlt$p[i,]*wt[i]
    #    apply(hmmlt$p,2,mean)
    scaled.density=h$density*sintegral(mean.p,hmmlt$x)/h.area
    for(i in 1:n) E.density[i,]=hmmlt$p[i,]*h.area/esw[i]
    mean.E.density=apply(E.density,2,mean)
  }
  
  if(is.null(xlab)) xlab="Perp. dist. (x)"
  nb=length(h$breaks)
  if(type=="prob") {
    # Detection function scale
    if(is.null(ylim) & allx) ylim=range(0,scaled.density,hmmlt$p)
    if(is.null(ylim) & !allx) ylim=range(0,scaled.density,mean.p)
    if(is.null(ylab)) ylab="p(x)"
    histline(scaled.density,h$breaks,ylim=ylim,xlab=xlab,ylab=ylab)
    if(allx & !is.nullmodel(hmmlt$models)) for(i in 1:n) lines(hmmlt$x,hmmlt$p[i,],col="gray")
    lines(hmmlt$x,mean.p,lwd=2)
    text(h$breaks[1:(nb-1)]+diff(h$breaks)/2,0,h$count,cex=text.cex)    
  } else if(type=="density") {
    # Count scale
    if(is.null(ylim) & allx) ylim=range(h$density,E.density)
    if(is.null(ylim) & !allx) ylim=range(h$density,mean.E.density)
    if(is.null(ylab)) ylab="Density"
    histline(h$density,h$breaks,ylim=ylim,xlab="Perp. dist. (x)",ylab="Density")
    if(allx & !is.nullmodel(hmmlt$models)) for(i in 1:n) lines(hmmlt$x,E.density[i,],col="gray")
    lines(hmmlt$x,mean.E.density,lwd=2)
    text(h$breaks[1:(nb-1)]+diff(h$breaks)/2,0,h$count,cex=text.cex)    
  } else stop(paste("Invalid type:",type))
  
  if(values) return(list(x=hmmlt$x,breaks=breaks,scaled.density=scaled.density,p=hmmlt$p,mean.p=mean.p,E.density=E.density,mean.E.density=mean.E.density))
}  


#----------------------- End Plotting functions --------------------------------------
#                        ----------------------
