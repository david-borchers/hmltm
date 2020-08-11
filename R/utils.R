#                      -----------------------------
#--------------------- Start DLB's utility functions --------------------
#' @title Clone of \code{sample} without weird behaviour for length(x)==1.
#'
#' @usage bsample(x, size, replace = FALSE, prob = NULL)
#' 
#' @description
#'  Identical to \code{sample} but when x is an integer it just returns x rather than an
#'  integer in 1:x (which is what \code{sample} does).
#'  
#' @param x Either a vector of one or more elements from which to choose, or a positive integer. 
#' See \link{sample} for details.
#' @param size a positive number, the number of items to choose from.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the elements of the vector being sampled.
#' 
bsample=function(x,size,replace=FALSE,prob=NULL) {
  if(length(x)==1) return(x)
  else return(sample(x,size,replace,prob))
}


#' @title Caclulates coefficient of variation.
#'
#' @description
#'  Utility function
#'  
#' @param x Random variable.
cv=function(x) sd(x)/mean(x) # calculates coefficient of variation

#' @title Converts nm to metres.
#'
#' @description
#'  Utility function
#'  
#' @param x Distance in nautical miles.
nm2m=function(x) return(x*1852) # converts nautical miles to metres

#mod=function(x,y) return(x-floor(x/y)*y) # returns the integer part of x/y


#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
histline=function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                  xlab="x",ylab="y",...)
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}


#' @title Integration using Simpson's method.
#'
#' @description
#'  Modified version of function sintegral from library Bolstad.
#'  
#' @param fx Values at \code{x} (see below) of the function to be integrated.
#' @param x Values at which function fx has been evaluated.
#' @param n.pts The number of points to be used in integration. If \code{x} contains more than
#'\code{n.pts} then \code{n.pts} will be set to \code{length(x)}.
#' @param type if equal to `\code{cdf}' a list comprising x-values and cdf values at each x-value 
#' is returned, else the integral over the range of \code{x} is returned.
sintegral=function (fx, x, n.pts=16, type="int") 
{
  #  if (class(fx) == "function") 
  #    fx = fx(x)
  #  n.x = length(x)
  #  if (n.x != length(fx)) 
  #    stop("Unequal input vector lengths")
  #  if (n.pts < 64) 
  #    n.pts = 64
  ap = approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h*(ap$y[2*(1:n.pts)-1]+4*ap$y[2*(1:n.pts)]+ap$y[2*(1:n.pts)+1])/3
  if(type!="cdf") return(sum(integral)) 
  else return(list(x=ap$x[2*(1:n.pts)],y=cumsum(integral)))
}



#' @title Logit function.
#'
#' @description
#'  Logit function
#'  
#' @param p probability (scalar or vector).
logit=function(p) return(log(p/(1-p)))  # returns logit of p

#' @title Inverse logit function.
#'
#' @description
#'  Inverse logit function
#'  
#' @param x scalar or .
inv.logit=function(x) return(1/(1+exp(-x))) # returns p from x=(logit of p)

#' @title Lognormal confidence interval with lower bound
#'
#' @usage lnci.nmin(n,Nhat,cv)
#' 
#' @description
#'  Calculates 95\% confidence interval for abundance estimate \code{Nhat}, based on the assumption that 
#'  \code{Nhat} is lognormally distributed and that \code{n} population members were observed (so that 
#'  lower bound of confidence interval can't be below \code{n}.
#'  
#' @param n observed number of population members.
#' @param Nhat point estimate of abundance.
#' @param cv coefficient of variation of \code{Nhat}.
#' 
#' @export
lnci.nmin=function(n,Nhat,cv){
  varNhat=(Nhat*cv)^2
  cfactor=exp(1.96*sqrt(log(1+varNhat/(Nhat-n)^2)))
  lower=n+(Nhat-n)/cfactor
  upper=n+(Nhat-n)*cfactor
  return(list(lower=lower,upper=upper))
}

#' @title Data truncation.
#'
#' @description
#' Left- and/or right-truncates any variable in data frame dat, inserting effort-only record if 
#' truncation removes all detections on a transect.
#' Subtracts left-trunction point off all variable values - after all truncation.

#' 
#' @param dat distance data frame. Must have columns for stratum stratum.area, transect, transect.length, 
#' <any-observer-level_variable>, and if \code{twosit}==TRUE then object as well. These can 
#' have any names, but the names must be specified, via argument \code{colnames}.
#' @param minvalue left-truncation point.
#' @param maxvalue right-truncation point.
#' @param twosit If TRUE, assumes this is an mrds-type dataset with two lines per detection.
#' @param colnames name of columns containing stratum, stratum area, transect, transect length, 
#' <any-observer-level_variable>, and if \code{twosit}==TRUE then object as well, IN THIS ORDER, 
#' in a character vector. The default value is 
#' colnames=c("stratum","area","transect","L","x","obs").
truncdat=function(dat,minval=0,maxval=NULL,twosit=FALSE,colnames=c("stratum","area","transect","L","x","obs")){
  tdat=dat
  keepcols=rep(NA,4)
  for(i in 1:4) {
    keepcols[i]=which(names(dat)==colnames[i])
    if(is.null(keepcols[i])) stop("No column in dat called ",colnames[i])
  }
  NAs=dat[1,,drop=FALSE]
  NAs[,-keepcols]=NA
  xcol=which(names(dat)==colnames[5])
  if(is.null(xcol)) stop("No column in dat called ",colnames[5])
  if(is.null(maxval)) maxval=max(na.omit(dat[,xcol]))
  if(twosit){
    if(length(colnames)<6) stop("With double-observer data, need 6th colname for observer; only 5 colnames:",colnames)
    obscol=which(names(dat)==colnames[6])
    if(is.null(obscol)) stop("No column in dat called ",colnames[6])
    out1=which(dat[,obscol]==1 & (dat[,xcol]<minval | dat[,xcol]>maxval))
    out2=which(dat[,obscol]==2 & (dat[,xcol]<minval | dat[,xcol]>maxval))
    if(length(out1) != length(out2)) stop("Different number of obs1 and obs2 detections for left-truncation")
    if(unique(out2-out1) != 1) stop("Looks like non-consecutive obs1 and obs2 detections for left-truncation")
    nout=length(out1)
    
    if(nout>0) {
      svalues=tvalues=rep(NA,nout)
      if(is.factor(dat[,keepcols[1]])) {
        svalues=rep(levels(dat[,keepcols[1]])[1],nout)
        levels(svalues)=levels(dat[,keepcols[1]])
        NAs[keepcols[1]]=svalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
        levels(NAs)=levels(svalues)
      }
      if(is.factor(dat[,keepcols[3]])) {
        tvalues=rep(levels(dat[,keepcols[3]])[3],nout)
        levels(tvalues)=levels(dat[,keepcols[3]])
        NAs[keepcols[3]]=tvalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
        levels(NAs)=levels(tvalues)
      }
      outeff=data.frame(stratum=svalues,area=rep(NA,nout),transect=tvalues,L=rep(NA,nout))
      if(is.factor(outeff$stratum)) levels(outeff$stratum)=levels(svalues)
      if(is.factor(outeff$transect)) levels(outeff$stratum)=levels(tvalues)
      for(i in 1:nout) outeff[i,]=dat[out1[i],keepcols] # store effort info for truncated sightings
      tdat=dat[-c(out1,out2),] # remove all truncated sightings
      for(i in 1:nout) {
        got=which(tdat[,keepcols[1]]==outeff$stratum[i] & tdat[,keepcols[2]]==outeff$area[i] & 
                    tdat[,keepcols[3]]==outeff$transect[i] & tdat[,keepcols[4]]==outeff$L[i])
        if(length(got)==0) { # transect no longer in truncated dataset
          tdat=rbind(tdat,NAs) # add row of NAs
          tdat[dim(tdat)[1],keepcols]=outeff[i,] # add in missing transect info
        }
      }  
    }
  } else {
    out=which(dat[,xcol]<minval | dat[,xcol]>maxval)
    nout=length(out)
    
    if(nout>0) {
      svalues=tvalues=rep(NA,nout)
      if(is.factor(dat[,keepcols[1]])) {
        svalues=rep(levels(dat[,keepcols[1]])[1],nout)
        levels(svalues)=levels(dat[,keepcols[1]])
        NAs[keepcols[1]]=svalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
        levels(NAs)=levels(svalues)
      }
      if(is.factor(dat[,keepcols[3]])) {
        tvalues=rep(levels(dat[,keepcols[3]])[3],nout)
        levels(tvalues)=levels(dat[,keepcols[3]])
        NAs[keepcols[3]]=tvalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
        levels(NAs)=levels(tvalues)
      }
      outeff=data.frame(stratum=svalues,area=rep(NA,nout),transect=tvalues,L=rep(NA,nout))
      if(is.factor(outeff$stratum)) levels(outeff$stratum)=levels(svalues)
      if(is.factor(outeff$transect)) levels(outeff$stratum)=levels(tvalues)
      #    for(i in 1:nout) outeff[i,]=dat[out1[i],keepcols] # store effort info for truncated sightings
      for(i in 1:nout) outeff[i,]=dat[out[i],keepcols] # store effort info for truncated sightings
      tdat=dat[-out,] # remove all truncated sightings
      for(i in 1:nout) {
        got=which(tdat[,keepcols[1]]==outeff$stratum[i] & tdat[,keepcols[2]]==outeff$area[i] & 
                    tdat[,keepcols[3]]==outeff$transect[i] & tdat[,keepcols[4]]==outeff$L[i])
        if(length(got)==0) { # transect no longer in truncated dataset
          tdat=rbind(tdat,NAs) # add row of NAs
          tdat[dim(tdat)[1],keepcols]=outeff[i,] # add in missing transect info
        }
      }
    }
  }
  
  if(nout>0) {
    tdat=tdat[order(tdat[keepcols[1]],tdat[keepcols[3]]),]
    tdat[!is.na(tdat[,xcol]),xcol]=tdat[!is.na(tdat[,xcol]),xcol]-minval # shift all left
  }
  return(tdat)
}


#' @title Reduces MRDS data frame to CDS data frame.
#'
#' @param dat distance data frame in mrds form. Must have columns 'object' (unique detection 
#' identifier), 'seen' (binary indicating detecteb by observer or not) and 'y' (forward) 
#' detection distance.
#' @param prefer which of the two observers' data to prefer when forward distances are 
#'   missing/equal must be 1 or 2.
#'   
#' @description
#' Reduces mark-recapture Distance sampling (MRDS) data frame dat, with two lines per detection, to a 
#' conventional distance sampling (CDS) data frame with a single line per detection. In the case of 
#' duplicates ("recaptures"), takes the information from the detection made farthest ahead (i.e. that
#' with larger \code{y}) With duplicates that have the same \code{y} or neither of which have a 
#' \code{y}, it chooses according to the parameter \code{prefer}. If only one of the duplicates 
#' has a \code{y}, it chooses that one.
#' 
#' @param dat MRDS data frame.
make.onesit=function(dat,prefer=1) {
  if(prefer!=1 & prefer!=2) stop("Argument 'prefer' must be 1 or 2.")
  n=dim(dat)[1]
  out=rep(FALSE,n)
  i=1
  while(i<=n){
    if(!is.na(dat$object[i])){
      if(dat$seen[i]==0) out[i]=TRUE
      if(dat$seen[i+1]==0) out[i+1]=TRUE
      if(dat$seen[i]==1 & dat$seen[i+1]==1) {
        if(is.na(dat$y[i]) & is.na(dat$y[i+1])) { # no y's; remove not preferred detection
          out[i+(3-prefer)-1]=TRUE
        } else if(is.na(dat$y[i]) & !is.na(dat$y[i+1])) { # keep only non-NA y
          out[i]=TRUE
        } else if(!is.na(dat$y[i]) & is.na(dat$y[i+1])){ # keep only non-NA y
          out[i+1]=TRUE
        } else if(dat$y[i]==dat$y[i+1]){# remove not preferred detection
          out[i+(3-prefer)-1]=TRUE
        } else if(dat$y[i]<dat$y[i+1]){ # remove later (closer) detection
          out[i]=TRUE
        } else {
          out[i+1]=TRUE
        }
      }
      i=i+2
    } else i=i+1
  }
  dat=dat[!out,] # remove worse observer's lines
  return(dat)
}



#' @title Makes survey parameter list (\code{survey.pars}) suitable for passing to 
#' \code{\link{est.hmltm}}.
#'
#' @description
#' \code{make.survey.pars} just puts a bunch of stuff in a list in a format that \code{\link{est.hmltm}}
#'  expects.
#'
#' @param spd is observer speed.
#' @param W is perpendicular right-truncation distance for estimation.
#' @param ymax is maximum forward distance for estimation - it should be far enough ahead that it is
#' not possible to detect anyting beyond \code{ymax}.
#' @param Wl is perpendicular left-truncation distance for estimation. (After truncating \code{Wl} is 
#' subtracted from all perpendicular distances.)
#' @param dT is the time step on which the availability hidden Markov model operates.
#' @param theta.f is REDUNDANT and should be removed. It must equal 0. 
#' @param theta.b is REDUNDANT and should be removed. It must equal 90. 
#'
#' @details Packs the above in a list suitable for passing as \code{survey.pars} to 
#' \code{\link{est.hmltm}}.
make.survey.pars=function(spd,W,ymax,Wl=0,dT=1,theta.f=0,theta.b=90){
  return(list(spd=spd,W=W,ymax=ymax,Wl=Wl,dT=dT,dy=spd*dT,theta.f=theta.f,theta.b=theta.b))
}

#' @title Decides if model is a null model.
#'
#' @description
#'  Logical function: true if \code{models} includes no covariates
#'  
#' @param models list of characters with elements \code{$y} and \code{$x} specifying y- and 
#' x-covariate models. Either \code{NULL} or regression model format (without response on left).
is.nullmodel=function(models){
  null=TRUE
  for(i in 1:length(models)) null=null & is.null(models[[i]])
  return(null)
}


#' @title Constructs linear predictor.
#'
#' @description
#' Returns parameter vector on linear predictor scale for hazard function \code{FUN}, using 
#' \code{model} and data frame \code{dat}. Makes a matrix with each row the parameter values for an 
#' observation, then concatenates these into a single vector (so can pass as vector to C++ code).
#'  
#' @param b parameter vector.
#' @param FUN detection hazard function name (character).
#' @param models list with two components (\code{$y} and \code{$x}) specifying models for detection 
#' hazard function scale parameters in forward (y) and perpendicular (x) dimensions. If detection 
#' hazard function form does not have separate scale parameters in y- and x- dimensions, the model
#' given in \code{$y} is take as the scale parameter model. \code{$y} and \code{$x} must be either 
#' \code{NULL} or "~<regspec>", where <regspec> is a regression model spefication (e.g. 
#' \code{height+weight} or \code{height:weight} or \code{height*weight}, etc.).
#' @param dat data frame, which must have columns corresponding to variable names in \code{$y} and
#' \code{$x}.
make.covb=function(b,FUN,models,dat)
{
  nfixed=switch(FUN,
                h.IP.0 = 2,
                h.EP1.0 = 2,
                h.EP2.0 = 3,
                h.EP1x.0 = 3,
                h.EP2x.0 = 4,
                0
  )
  if(nfixed==0) stop("Hazard model ",FUN," is not progremmed (yet).")
  if(is.nullmodel(models)) {
    if(length(b) != nfixed) stop("length of b inconsistent with model")
    n=dim(dat)[1]
    #    covb=matrix(c(rep(b,rep(n,nfixed))),ncol=nfixed)
    covb=rep(b,n)
  } else {
    if(FUN=="h.EP2.0" | FUN=="h.EP1.0" | FUN=="h.IP.0"){
      X=model.matrix(as.formula(models$y),data=dat)
      n=dim(X)[1]
      nb=dim(X)[2]
      nfixed=nfixed-1
      if(length(b) != (nfixed+nb)) stop("length of b inconsistent with model")
      covb=matrix(c(rep(b[1:nfixed],rep(n,nfixed)),X%*%b[nfixed+(1:nb)]),ncol=(nfixed+1))
      covb=as.vector(t(covb))
    }
    else if(FUN=="h.EP2x.0" | FUN=="h.EP1x.0"){
      if(!is.null(models$y)){
        X=model.matrix(as.formula(models$y),data=dat)
        n=dim(X)[1]
        nb=dim(X)[2]        
        if(!is.null(models$x)){ # here if have covariates for x and y
          X.x=model.matrix(as.formula(models$x),data=dat)
          nb.x=dim(X.x)[2]
          nfixed=nfixed-2
          if(length(b) != (nfixed+nb+nb.x)) stop("length of b inconsistent with model and model.x")
          covb=matrix(c(rep(b[1:nfixed],rep(n,nfixed)),X%*%b[nfixed+(1:nb)],X.x%*%b[nfixed+nb+(1:nb.x)]),ncol=(nfixed+2))  
          covb=as.vector(t(covb))
        } else { # here if have covariates for y only
          nfixed=nfixed-1
          if(length(b) != (nfixed+nb)) stop("length of b inconsistent with model and model.x")
          #      covb=matrix(c(rep(b[1:nfixed],rep(n,nfixed)),X%*%b[nfixed+(1:nb)]),ncol=(nfixed+1))   #### bitchange
          covb=matrix(c(rep(b[1:(nfixed-1)],rep(n,(nfixed-1))),X%*%b[(nfixed-1)+(1:nb)],rep(b[length(b)],n)),ncol=(nfixed+1))  #### bitchange 
          covb=as.vector(t(covb))
        }
      } else { # here if have covariates for x only
        X.x=model.matrix(as.formula(models$x),data=dat)
        n.x=dim(X.x)[1]
        nb.x=dim(X.x)[2]
        nfixed=nfixed-1
        if(length(b) != (nfixed+nb.x)) stop("length of b inconsistent with model and model.x")
        covb=matrix(c(rep(b[1:nfixed],rep(n.x,nfixed)),X.x%*%b[nfixed+(1:nb.x)]),ncol=(nfixed+1))  
        covb=as.vector(t(covb))        
      }
    } else {
      stop("Invalid FUN.")
    }
  }
  return(covb)
}
