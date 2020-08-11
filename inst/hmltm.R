# Changed all pcu to pm 6th May 2013. Not yet tested with bootstrap.
#
#
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


# ============== Other availability correction methods and realted ==========================

#' @title Laake's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the probability that an animal is available at least once while within detectable
#' region, using the method of Laake et al. (1997), Eqn (4). Does this for each of a set of m 
#' 2-state Markov model availability parameters passed to it and returns probabilities and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param ymax is maximum forward distance to consider (`w' in Laake et al. (1997), Eqn (4)).
#' @param spd speed of observer. Needed to convert time (units of \code{hmm.pars}) to distance (units
#'  of \code{ymax}).
#'
#' @details See Laake et al. (1997), Eqn (4), or Borchers et al. (2013) for details of Laake's method.
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' Laake, J., Calambokidis, J., Osmek, S., and Rugh, D. 1997. Probability of detecting harbor 
#' porpoise from aerial surveys: estimating g(0). Journal of Wildlife Management 61, 63-75.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{mclaren.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' laake.a(hmm.pars,ymax=200,spd=4) # because hmm.pars is in units of time, not distance, you need to specify spd.
laake.a=function(hmm.pars,ymax,spd=NULL){
  if(length(dim(hmm.pars$Pi))==2) hmm.pars$Pi=array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  nav=dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a=rep(NA,nav)
  for(i in 1:nav){
    a[i]=jeffa(hmm.pars$Pi[,,i],ymax,spd)
  }
  return(list(mean=mean(a),a=a))
}


#' @title Laake's availability correction factor calculation for a single 2-state Markov model.
#'
#' @description
#' Calculates the probability that an animal is availabel at least once whil within detectable
#' region, using the method of Laake et al. (1997), Eqn (4).
#'
#' @param Pi is Markov model transition matrix (state 1 is UNavailable).
#' @param w is maximum forward distance things can be detected.
#' @param spd is observer speed. If NULL w is assumed to be in the units of the Markov chain.
#' @param E is mean times unavailable (E[1]) and available (E[2]) (in unist of Markov chain). 
#' If E is not NULL, it is used in preference to Pi for calculations, else it is ignored.
#'
#' @details See Laake et al. (1997), Eqn (4), or Borchers et al. (2013) for details of the method.
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' Laake, J., Calambokidis, J., Osmek, S., and Rugh, D. 1997. Probability of detecting harbor 
#' porpoise from aerial surveys: estimating g(0). Journal of Wildlife Management 61, 63-75.
jeffa=function(Pi,w,spd=NULL,E=NULL)
{
  if(!is.null(spd)) w=w/spd # convert distance to time if spd given
  if(!is.null(E)) {
    warning("Used E, not Pi for calculations")
  } else {
    E=makeE(Pi) # expected time up, expected time down
  }
  a=(E[2] + E[1]*(1-exp(-w/E[1])))/sum(E)    
  return(a)
}


#' @title Instantaneous availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the probability that an animal is available at an instant, for each of a set of m 
#' 2-state Markov model availability parameters passed to it and returns probabilities and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param Et is a 2xm matrix with expect times Unavailable (row 1) and Available (row 2).
#'
#' @details If Et is given, it is used and hmm.pars is ignored.
#' 
#' @seealso \code{\link{mclaren.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3);pm=NULL
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' instant.a(hmm.pars)
#' instant.a(NULL,Et=matrix(c(Eu,Ea),ncol=2,byrow=TRUE))
#' instant.a(NULL,c(Eu[1],Ea[1]))
instant.a=function(hmm.pars,Et=NULL){
  if(!is.null(Et)) {
    if(!is.null(hmm.pars)) warning("Used Et, not Pi for calculations")
    if(is.vector(Et)) Et=matrix(Et,ncol=1)
    if(dim(Et)[1]!=2) stop("1st dimension of Et must be 2.")
    nav=dim(Et)[2]
  }else {
    if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) stop("1st two dimensions of hmm.pars$Pi must be 2.")
    if(length(dim(hmm.pars$Pi))==2) hmm.pars$Pi=array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
    nav=dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  }
  a=rep(NA,nav)
  for(i in 1:nav){
    if(is.null(Et)) a[i]=simplea(hmm.pars$Pi[,,i],Et)
    else a[i]=simplea(NULL,Et[,i])
  }
  return(list(mean=mean(a),a=a))
}


#' @title Simple availability correction factor calculation.
#'
#' @description
#' Calculates proportion of time an animal with a Markov availability process is available.
#'
#' @param Pi is Markov model transition matrix (state 1 is UNavailable).
#' @param E is expected time in each state (state 1 in UNavailable).
#'
#' @details If \code{E} is NULL, uses Pi to calculate proportion of time available, else uses \code{E}.
simplea=function(Pi,E=NULL)
{
  if(is.null(E)) E=makeE(Pi) # expected time up, expected time down
  a=E[2]/sum(E)    
  return(a)
}

#' @title McLaren's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates McLaren's availability correction factor, for each of a set of m 2-state Markov model 
#' availability parameters passed to it and returns these and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element \code{$Pi} (where m is number of availability parameter sets).
#' @param w is max forward distance things can be seen at (or max forward time). Must be scalar.
#' @param spd is observer speed; omit if w is max forward TIME.
#'
#' @references
#' McLaren, I.A. 1961. Methods of determining the numbers and availability of ringed seals in the 
#' eastern Canadian Arctic. Arctic 14:162 – 175.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' mclaren.a(hmm.pars,w=10,spd=4)
#' mclaren.a(hmm.pars,w=100,spd=4) # can be greater than 1 (!)
mclaren.a=function(hmm.pars,w,spd=1){
  if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) stop("1st two dimensions of hmm.pars$Pi must be 2.")
  if(length(dim(hmm.pars$Pi))==2) hmm.pars$Pi=array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  nav=dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a=rep(NA,nav)
  for(i in 1:nav){
    Eti=makeE(hmm.pars$Pi[,,i])*spd # Pi is TIME; if w is DISTANCE need to multiply time by speed
    a[i]=(Eti[2]+w)/sum(Eti)
  }
  return(list(mean=mean(a),a=a))
}

#' @title Richard's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the availability correction factor "C_{ca}" of Richard et al. (2010). Does this for each 
#' of a set of m 2-state Markov model availability parameters passed to it and returns probabilities 
#' and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param w vector of forward distances.
#' @param spd observer speed: must be entered if y is not time, since hmm.pars always time.
#'
#' @details See Richard et al. (2010), equation on botto mof page 91 for details of method.
#' 
#' @references 
#' Richard, P.R., Laake, J.L., Hobbs, R.C., Heide-Jorgensen, M.P., Asselin, N.C. and Cleator, H. 2010.
#' Baffin Bay narwhal population distribution and numbers: aerial surveys in the Canadian High Arctic,
#' 2002-04. Arctic 63: 85-99.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{mclaren.a}}, \code{\link{laake.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3) # mean avail and unavail times (& var)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt) # make hmm.pars object
#' richard.a(hmm.pars,w=10,spd=4)
#' richard.a(hmm.pars,w=100,spd=4) # can be greater than 1 (!)
#' richard.a(hmm.pars,w=rexp(20,1/100),spd=4)
richard.a=function(hmm.pars,w,spd=1){
  y=na.omit(w)
  n=length(y)
  if(length(dim(hmm.pars$Pi))==2) hmm.pars$Pi=array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  nav=dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  ina=matrix(rep(instant.a(hmm.pars)$a,n),nrow=n,byrow=TRUE)
  mca=matrix(rep(NA,n*nav),nrow=n)
  for(i in 1:n) mca[i,]=mclaren.a(hmm.pars,y[i],spd)$a
  sumfb=apply(ina/mca,2,sum)
  a=1/((1/ina[1,])*(sumfb/n))
  return(list(mean=mean(a),a=a))
}


#' @title Makes Markov transition matrices from mean times available and unavailable.
#'
#' @description
#' Makes Markov transition matrices from mean times available (Ea) and unavailable (Eu). 
#' If Ea and Eu are vectors of length m (they must be the same length), returns a 2x2xm array 
#' in which element [,,i] is the ith Markov transition matrix; else returns a single 2x2 matrix.
#'
#' @param Eu is the mean time UNavailable in one available-unavailable cycle.
#' @param Ea is the mean time available in one available-unavailable cycle.
#'
#' @examples
#' Ea=c(10,12);Eu=c(20,22)
#' makePi(Eu,Ea)
#' makePi(Eu[1],Ea[1])
makePi=function(Eu,Ea)
{
  nav=length(Eu)
  if(length(Ea)!=nav) stop("Lengths of Eu and Ea must be the same")
  #  if(nav==1) {
  #    Pi=matrix(rep(0,4),nrow=2,
  #              dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available")))
  #    Pi[1,2]=1/Eu
  #    Pi[2,1]=1/Ea
  #    Pi[1,1]=1-Pi[1,2]
  #    Pi[2,2]=1-Pi[2,1]
  #  } else {
  Pi=array(rep(0,2*2*nav),dim=c(2,2,nav),
           dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                         Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Pi[1,2,i]=1/Eu[i]
    Pi[2,1,i]=1/Ea[i]
    Pi[1,1,i]=1-Pi[1,2,i]
    Pi[2,2,i]=1-Pi[2,1,i]    
  }
  #  }
  return(Pi)
}

#' @title Returns expected time in state for 2-state Markov transition matrices.
#'
#' @description
#' Returns expected time in state for the m 2-state Markov transition matrices Markov transition 
#' matrices contained in the matrix or array Pi. If m>1 returns a matrix with column i being the
#' expected time in state 1 ("unavailable") as its first element and expected time in state 2 
#' ("available") as its second; else returns a vector.
#'
#' @param Pi either a 2x2 Markov transition probability matrix or a 2x2xm array with element [,,i]
#' being a 2x2 Markov transition probability matrix.
#'
#' @examples
#' Ea=c(10,12);Eu=c(20,22)
#' Pi2=makePi(Eu,Ea) # make a set of transition matrices
#' makeE(Pi2) # recover c(Eu,Ea)
#' Pi1=makePi(Eu[1],Ea[1]) # make a single transition matrix
#' makeE(Pi1) # recover c(Eu,Ea)
makeE=function(Pi){
  #----------------------------------------------------------
  # Returns expected time in states 1 and 2 for the 2x2 
  # probability transition matrix Pi for a 2-state
  # Markov process.
  #----------------------------------------------------------
  if(length(dim(Pi))==2){
    E=c(1/Pi[1,2],1/Pi[2,1])
    names(E)=c("Unavailable","Available")
  } else {
    nav=dim(Pi)[3]
    E=matrix(rep(NA,2*nav),nrow=2,
             dimnames=list(State=c("Unavailable","Available"),Animal=as.character(1:nav)))
    for(i in 1:nav){
      E[,i]=c(1/Pi[1,2,i],1/Pi[2,1,i])      
    }
  }
  return(E)
}
  

#' @title Makes 2-state hidden Markov model list suitable for passing to 
#' \code{\link{est.hmltm}} via argument \code{hmm.pars}.
#'
#' @description
#' \code{make.hmm.pars.from.Et} creates 2-state hidden Markov availability model specification from 
#' mean times unavailable (Eu) and available (Ea).
#'
#' @param Ea vector of length m>=1 specifying mean distance (=mean time * observer speed) animals are in 
#' the more available state in one cycle (e.g. dive cycle: surface-dive).
#' @param Eu vector of length m>=1 specifying mean distance (=mean time * observer speed) animals are in 
#' the more UNavailable state in one cycle.
#' @param seEa standard error of Ea.
#' @param seEu standard error of Eu.
#' @param covEt vector of length m>=1 containing covariance of each pair of Ea and Eu.
#' @param pm is 2xm matrix containing state-dependent Bernoulli distribution parameters for the m
#' pairs of Ea and Eu, with first being probability of being available when in state i (i=1,2), 
#' where i=1 is the more UNavailable state and i=2 is the more available state.
#'
#' @details Calculates 2-state hidden Markov model parameters such that the Markov process is in states 1 
#' (more unavailable) and 2 (more available) for mean times Ea and Eu, with Bernoulli state-dependent response
#' probability pm[1,] and pm[2,], respectively. Also constructs a covariance matrix for Ea and Eu. If 
#' pm[1,i]=0 and pm[2,i]=1 for animal i, the availability process is a Markov process and the states 
#' are actual unavailbile (state 1) and availabile (state 2).
#' 
#'  @examples
#'  # Some arbitrary numbers for illustration:
#'  Ea=c(10,12);Eu=c(100,120);seEa=c(2,3);seEu=c(20,30);covEt=c(10,12)
#'  make.hmm.pars.from.Et(Ea[1],Eu[1],seEa[1],seEu[1],covEt[1]) # single animal
#'  make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt) # two animals
#'  
#'  # Here's how the data porpoise.hmm.pars was created from numbers in Westgate et al. (1995):
#'  ppn=c(40,52,36,34,49,60,33)/100 # proportion of time available
#'  ET=c(76,44,52,64,70,46,103) # mean dive cycle duration
#'  seET=c(48,37,52,65,59,32,67) # SE of mean dive cycle duration
#'  cvET=seET/ET # CV of mean dive cycle duration
#'  Ea=ET*ppn  # mean time available
#'  Eu=ET*(1-ppn)  # mean time UNavailable
#'  # For lack of better info, assume independence of Ea and Eu, and that cv(Ea)=cv(Eu)=cv, 
#'  # which means that
#'  cv=sqrt((cvET*ET)^2/(Ea^2+Eu^2))
#'  # and hence:
#'  seEa=Ea*cv
#'  seEu=Eu*cv
#'  covEt=seET*0 # assume independence
#'  porpoise.hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#'  
#'  # Here's how the dataset beaked.hmm.pars were created:
#'  Ea=121.824
#'  Eu=1580.256
#'  seEa=9.618659
#'  seEu=134.9212
#'  beaked.hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu)
#'  
#' @references 
#' Westgate, A. J., Read, A. J., Berggren, P., Koopman, H. N., and Gaskin, D. E. 1995. Diving behaviour 
#' of harbour porpoises, Phocoena phocoena. Canadian Journal of Fisheries and Aquatic Sciences 52, 
#' 1064-1073.
#' 
make.hmm.pars.from.Et=function(Ea,Eu,seEa,seEu,covEt=0,pm=NULL) {
  nav=length(Ea)
  if(length(Eu)!=nav |length(seEa)!=nav |length(seEu)!=nav |length(covEt)!=nav) stop("Lengths of Ea, Eu, seEa, seEu, covEt must all be the same.")
  if(is.null(pm)) pm=matrix(c(rep(0,nav),rep(1,nav)),nrow=2,byrow=TRUE)
  if(is.vector(pm)) {
    if(length(pm)!=2) stop("pm must either be a vector of length 2 or a matrix of dimension length(Ea)x2.")
    pm=matrix(c(pm[1],pm[2]),ncol=2)
  }
  if(dim(pm)[2]!=nav) stop("Inconsistent dimensions of Ea and pm.")
  Pi=Sigma.Et=array(rep(NA,2*2*nav),dim=c(2,2,nav),
                    dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                                  Animal=as.character(1:nav)))
  Et=delta=newpm=matrix(rep(NA,2*nav),ncol=nav,dimnames=list(State=c("Unavailable","Available"),
                                                          Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Et[,i]=c(Eu[i],Ea[i])
    Sigma.Et[,,i]=diag(c(seEu[i],seEa[i])^2)
#    cvEt=c(seEu[i]/Et[1,i],seEa[i]/Et[2,i])
#    Sigma.Et[,,i]=diag((cvEt*Et)^2)
    Sigma.Et[1,2,i]=Sigma.Et[2,1,i]=covEt[i]
    Pi[,,i]=makePi(Et[1,i],Et[2,i])
    delta[,i]=compdelta(Pi[,,i])
    newpm[,i]=pm[,i]
  }
  hmm.pars=list(pm=newpm,Pi=Pi,delta=delta,Et=Et,Sigma.Et=Sigma.Et)
  return(hmm.pars)  
}


#' @title Calculate lengths of all runs in availability.
#'
#' @description
#' Calculates lengths of the runs of 0s and 1s in a vector.
#'
#' @param avail is a vecto of binary data (0s and 1s).
#' 
#' @details Returns a list with components \code{run1} being a vector with 
#' the lengths of all runs of 1s, in order, and  \code{run0} being a vector with 
#' the lengths of all runs of 0s, in order.
#' 
getruns = function(avail) {
  n = length(avail)
  run1 = run0 = NULL
  state = avail[1]
  runlen = 1
  for(i in 2:n) {
    if(avail[i]==state) runlen = runlen+1
    else {
      if(state==1) run1 = c(run1,runlen)
      else run0 = c(run0,runlen)
      state = avail[i]
      runlen = 1
    }
  }
  return(list(run1=run1,run0=run0))
}

  
# ==============- utility functions specific to avail estimation ------------------

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



#' @title Constructs a HMM equivalent to a Poisson process.
#'
#' @description
#'  Returns list with $Pi, $pm and $delta of Poisson with same event rate as in the input HMM object 
#'  availhmm.
#'  
#' @param availhmm availability model list of the sort passed to \code{\link{est.hmltm}}.
#' @param zero REDUNDANT (I think - CHECK)
poiss.equiv=function(availhmm,zero=0){
  Pois.availhmm=availhmm
  Pi=availhmm$Pi
  pcu=FALSE
  if(is.element("pcu",names(availhmm))) {
    names(availhmm)[which(names(availhmm)=="pcu")]="pm"
    pcu=TRUE
  }
  pm=availhmm$pm
  delta=availhmm$delta
  if(is.vector(pm)&!is.matrix(Pi) | !is.vector(pm)&is.matrix(Pi)) stop("Single animal: pcu/pm is not a vector or Pi is not a matrix")
  if(is.vector(pm)) { # convert to matrix and array so can use loop below
    Pi=array(Pi,dim=c(2,2,1))
    pm=matrix(pm,ncol=1)
    delta=matrix(delta,ncol=1)
  }
  nw=dim(pm)[2]
  
  PiPoiss=matrix(c(zero,1-zero,zero,1-zero),byrow=TRUE,nrow=2)
  for(i in 1:nw) {
    Pi[,,i]=PiPoiss
    delta[,i]=compdelta(PiPoiss)
  }
  E=Estate=matrix(rep(0,2*nw),nrow=2)
  if(nw>1) {
    for(w in 1:nw) {
      Estate[,w]=c(1/availhmm$Pi[1,2,w],1/availhmm$Pi[2,1,w])
      events=apply(Estate*pm,2,sum)
      duration=apply(Estate,2,sum)
      eventrate=rep(0,length(duration))
      eventrate[duration>0]=events[duration>0]/duration[duration>0]
      Pois.availhmm$Pi=Pi
      if(pcu) Pois.availhmm$pcu=matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      else Pois.availhmm$pm=matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      Pois.availhmm$delta=delta
    }
  }else {
    Estate=c(1/availhmm$Pi[1,2],1/availhmm$Pi[2,1])
    events=sum(Estate*pm)
    duration=sum(Estate)
    eventrate=events/duration
    Pois.availhmm$Pi=Pi[,,1]
    if(pcu) Pois.availhmm$pcu=c(0,eventrate)
    else Pois.availhmm$pm=c(0,eventrate)
    Pois.availhmm$delta=as.vector(delta)
  }
  
  return(Pois.availhmm)
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

#--------------------- End DLB's utility functions --------------------
#                      ---------------------------


#                        --------------------------------
#----------------------- Start detection hazard functions --------------------------------------
#                      (C++ code makes this all redundant)

#h.EP1=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.EP1(b)
#  return(par[1]*exp(-(abs(x)^par[2] + y^par[2])/(par[3]^par[2])))
#}
#tfm.EP1=function(par) return(c(logit(par[1]),log(par[2:3]))) # Transforms EP1 parameters from natural to link scale
#invtfm.EP1=function(b) return(c(inv.logit(b[1]),exp(b[2:3]))) # Transforms EP1 parameters from link to natural scale

#h.EP1.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #----------------------------------------------------------
#{
#  par=invtfm.EP1.0(b)
#  return(exp(-(abs(x)^par[1] + abs(y)^par[1])/(par[2]^par[1])))
#}
#tfm.EP1.0=function(par) return(log(par))
#invtfm.EP1.0=function(b) return(exp(b))

#h.EP1x.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #
  # Modified to have separate x- and y- scale functions
  #----------------------------------------------------------
#{
#  par=invtfm.EP1x.0(b)
#  return(exp(-((abs(x)/par[3])^par[1] + (abs(y)/par[2])^par[1])))
#}
#tfm.EP1x.0=function(par) return(log(par))
#invtfm.EP1x.0=function(b) return(exp(b))

#h.EP2=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.EP2(b)
#  return(par[1]*exp(-((abs(x)/par[4])^par[2] + (abs(y)/par[4])^par[3])))
#}
#tfm.EP2=function(par) return(c(logit(par[1]),log(par[2:4])))
#invtfm.EP2=function(b) return(c(inv.logit(b[1]),exp(b[2:4])))

#h.EP2.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #----------------------------------------------------------
#{
#  par=invtfm.EP2.0(b)
#  return(exp(-((abs(x)/par[3])^par[1] + (abs(y)/par[3])^par[2])))
#}
#tfm.EP2.0=function(par) return(log(par))
#invtfm.EP2.0=function(b) return(exp(b))

#h.EP2x.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #
  # Modified to have separate x- and y- scale functions
  #----------------------------------------------------------
#{
#  par=invtfm.EP2x.0(b)
#  return(exp(-((abs(x)/par[4])^par[1] + (abs(y)/par[3])^par[2])))
#}
#tfm.EP2x.0=function(par) return(log(par))
#invtfm.EP2x.0=function(b) return(exp(b))

#h.IP=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.IP(b)
#  return(par[1]*exp(par[2]*log(par[3])-(par[2]/2)*log(par[3]^2+x^2+y^2)))
#}
#tfm.IP=function(par) return(c(logit(par[1]),log(par[2:3])))
#invtfm.IP=function(b) return(c(inv.logit(b[1]),exp(b[2:3])))
                                                            
#h.IP.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.IP.0(b)
#  return(exp(par[1]*log(par[2])-(par[1]/2)*log(par[2]^2+x^2+y^2)))
#}
#tfm.IP.0=function(par) return(log(par))
#invtfm.IP.0=function(b) return(exp(b))

# Derivative stuff below from earlier version: now redundant
#dinvt.dpar=function(b,fun)
#{
#  switch(fun,
#         h.EP1=d.dp.EP1(b),
#         h.EP1.0=d.dp.EP1.0(b),
#         h.EP1x.0=d.dp.EP1x.0(b),
#         h.EP2=d.dp.EP2(b),
#         h.EP2.0=d.dp.EP2.0(b),
#         h.EP2x.0=d.dp.EP2x.0(b),
#         h.IP=d.dp.IP(b),
#         h.IP.0=d.dp.IP.0(b),
#  )
#}
#d.dp.EP1=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:3]))) # Returns derivative w.r.t p of link functions at parameter estimates (parameters on link scale)
#d.dp.EP1.0=function(b) return(exp(b))
#d.dp.EP1x.0=function(b) return(exp(b))
#d.dp.EP2=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:4])))
#d.dp.EP2.0=function(b) return(exp(b))
#d.dp.EP2x.0=function(b) return(exp(b))
#d.dp.IP=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:3])))
#d.dp.IP.0=function(b) return(exp(b))

#d.dx.inv.logit=function(x) return(exp(-x)/(1+exp(-x))^2) # dervivative of inverse logit w.r.t. x=logit(p)

#----------------------- End detection hazard functions --------------------------------------
#                        ------------------------------


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


#                          -----------------------------------------
#------------------------- Start Derived Stats Calculation functions -----------------------

#' @title Calculates perp dist detection probability p(x).
#'
#' @description
#' Calculates perp dist detection probability p(x).
#'
#' @param x perpendicular distances at which to evaluate function.
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param cov covariate matrix with each row corresponding to an observation. Must contain columns 
#' with variables appearing in \code{models}, and named accordingly, as well as column of perpendicular
#' distance, named "x". (Perpendicular distances only used if \code{to.x} is TRUE.)
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars availability hmm parameters, as for \code{\link{est.hmltm}}. Must have elements
#' \code{$Et} and \code{$Sigma.Et}
#' @param type if "link", parameter vector \code{pars} is assumed to be on the link scale, else on 
#' the natural scale
#' @param ally If TRUE calculates detection probability at all forward distances, else at zero.
#' 
hmltm.px=function(x,pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,
                  type="response",ally=TRUE){
  theta.f=survey.pars$theta.f
  theta.b=survey.pars$theta.b
  ymax=survey.pars$ymax
  dy=survey.pars$dy
  pm=hmm.pars$pm
  Pi=hmm.pars$Pi
  delta=hmm.pars$delta
  
  b=pars
  if(type!="link") b=tfm(pars,hfun)
  
  nx=length(x)
  n=1
  covb=b
  if(!is.null(cov) & !(is.nullmodel(models))) { # only use covars if have them and model uses them
    n=dim(cov)[1]
    covb=make.covb(b,hfun,models,cov) # put covariates into paramerters
  } 
  nb=length(covb)/n
  px=matrix(rep(0,nx*n),nrow=n)
  for(i in 1:n) {
    start=(i-1)*nb+1
    bi=c(rep(covb[start:(start+nb-1)],nx)) # nx replicates of covb for ith detection
    px[i,]=p.xy(x=x,y=rep(0,nx),hfun=hfun,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,ally=ally)
  }
  return(px)
}


#' @title Calculates a bunch of derived statistics from model.
#'
#' @description
#' Calculates a one of a variety of derived statistics from model (see below).
# 
#' @param stat the statistic name (character variable).
#' @param hmmlt the model (as ouput by \code{\link{est.hmltm}}.
#' @param obs indices of rows of \code{hmmlt$xy} (i.e. which observations) to use in calculating 
#' the statistics.
#' 
#' @details
#' The following are the options for argument \code{stat}:
#' \describe{
#' \item{esw:}{effective strip width estimate.}
#' \item{invesw:}{inverse effective strip width estimate.}
#' \item{p0:}{estimated probability of detection at perpendicular distance zero.}
#' \item{p:}{estimated mean probability of detection.}
#' \item{invp:}{estimated inverse mean probability of detection.}
#' }
calc.derived=function(stat,hmmlt,obs=1:dim(hmmlt$xy)[1]){
  if(stat=="esw") {return(fitted.esw(hmmlt,obs))}
  else if(stat=="invesw") {return(fitted.invesw(hmmlt,obs))}
  else if(stat=="p0") {return(fitted.px(hmmlt,obs,at.x=0))}
  else if(stat=="p") {return(fitted.p(hmmlt,obs))}
  else if(stat=="invp") {return(fitted.invp(hmmlt,obs))}
  else stop(paste(stat," is an invalid stat type"))
}


#' @title Calculates fitted values for p(x) from model.
#'
#' @description
#' Calculates fitted values, p(x) for given observations, from model (optionally at given x-value).
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param at.x values at which to evaluate p(x) - see below.
#' 
#' @details 
#' If \code{hmmlt$models} is NULL and
#' \describe{
#'   \item{at.x is NULL}{returns vector of values of p(x=hmmlt$xy$x[obs]) if obs given, or vector 
#'   of values of p(x=hmmlt$xy$x) if obs not given;}
#'   \item{at.x is specified}{returns vector of values of p(x=at.x).}
#' }
#' If \code{hmmlt$models} is not NULL and
#' \describe{
#'   \item{at.x is NULL}{returns vector of values of p_obs[i](x=hmmlt$xy$x[j]) if obs given, or 
#'   vector of values of p_i(x=hmmlt$xy$x[j]) for all i if obs not given;}
#'   \item{at.x is specified}{as when at.x is NULL, but with x=at.x instead of hmmlt$xy$x.}
#'}
fitted.px=function(hmmlt,obs=1:dim(hmmlt$xy)[1],at.x=NULL){
  cov=hmmlt$xy
  if(max(obs)>dim(hmmlt$xy)[1]) stop("obs greater than number observations in hmmlt$xy")
  if(min(obs)<1) stop("obs < 1")
  cov=cov[obs,]
  if(!is.null(at.x)) {
    if(length(at.x)!=1 & length(at.x)!=dim(cov)[1]) stop("Length of at.x inconsistent with covariate data frame.")
    if(length(at.x)==1) {cov$x=rep(at.x,length(cov$x));at.x=cov$x}
    else {cov$x=at.x}
  } else {
    at.x=cov$x
  }
  if(is.nullmodel(hmmlt$models)){
    cov=cov[1,]
    #    at.x=at.x[1]
  }
  pars=hmmlt$fit$par
  hfun=hmmlt$h.fun
  models=hmmlt$models
  survey.pars=hmmlt$fitpars$survey.pars
  hmm.pars=hmmlt$fitpars$hmm.pars
  px=hmltm.px(at.x,pars,hfun,models,cov,survey.pars,hmm.pars,ally=TRUE)
  if(!is.nullmodel(hmmlt$models)){
    #    rownames(px)=paste("obs",obs,sep="")
    #    colnames(px)=paste("x=",at.x,sep="")
    px=diag(px)
  }
  return(px)
}

#' @title Calculates E[p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of fitted values, p(x) for given observations, 
#' from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
fitted.p=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100){
  W=hmmlt$fitpars$survey.pars$W
  esw=fitted.esw(hmmlt,obs,nx)
  return(esw/W)
}

#' @title Calculates E[1/p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of inverse of fitted values, 1/p(x) for given 
#' observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
fitted.invp=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100){
  W=hmmlt$fitpars$survey.pars$W
  esw=fitted.esw(hmmlt,obs,nx)
  return(W/esw)
}

#' @title Calculates 1/esw from model.
#'
#' @description
#' Calculates inverse of effective strip half-width, 1/(W*E[p(x)]) (where W is actual half-width) for 
#' given observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
fitted.invesw=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100){
  return(1/fitted.esw(hmmlt,obs,nx))
}

#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip half-width, W*E[p(x)], where W is actual half-width, from model.
# 
#'  @param hmmlt output from \code{\link{fit.hmltm}}
#'  @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#'  @param nx number of x values to use to implement Simpson's rule in perp dist dimension;
#'  @param to.x If TRUE integrates only to observed x, else integrates to W
#'  @param all If TRUE then returns esw for every observation, else returns only that for first obs 
#'  if there are no covariates; always returns esw for every observation if there are covariates.
#' 
#' @details 
#' Calls \code{\link{hmltm.esw}} to calclate effective stript width (esw) for fitted object \code{hmmlt}.
#
fitted.esw=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100,to.x=FALSE,all=FALSE){
  if(!is.null(obs)){
    if(max(obs)>dim(hmmlt$xy)[1]) stop("obs greater than number observations in hmmlt$xy")
    if(min(obs)<1) stop("obs < 1")
    cov=hmmlt$xy[obs,]
  } 
  if(is.nullmodel(hmmlt$models) & !all){
    cov=hmmlt$xy[1,]
  } else {
    cov=hmmlt$xy[obs,]
  }
  pars=hmmlt$fit$par
  hfun=hmmlt$h.fun
  models=hmmlt$models
  survey.pars=hmmlt$fitpars$survey.pars
  hmm.pars=hmmlt$fitpars$hmm.pars
  esw=hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type="response",to.x)
  return(esw)
}


#' @title Calculates various statistics from model.
#'
#' @description
#' Calculates various statistics from model.
# 
#' @param stat name of statistic to calculate. Valid statistics are "p0" for estimated probability 
#' at perpendicular distance zero, "p" for mean estimated detection probability over all perpendicular
#' distances, "invp" for the inverse of mean estimated detection probability over all perpendicular
#' distances, "esw" for estimated effective strip width, and "invesw" for estimated inverse of effective 
#' strip width.
#' @param b detection hazard parameter vector.
#' @param hfun detection hazard name (character).
#' @param models covariate models (see \code{\link{est.hmltm}} for details).
#' @param cov covariate values.
#' @param survey.pars survey parameter specification  (see \code{\link{est.hmltm}} for details).
#' @param hmm.pars hidden Markov model parameter specification  (see \code{\link{est.hmltm}} for 
#' details).
#' @param nx number of points at which to evaluate detection function in perpendicular distance
#' dimension.
#' @param type if "link", assumes that parameter vector \code{b} is on link scale, else assumes 
#' it is on natural scale.
hmltm.stat=function(stat,b,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                    type="link"){
  if(type!="link") b=tfm(b,hfun)
  if(stat=="p0") {return(hmltm.px(x=0,b,hfun,models,cov,survey.pars,hmm.pars,type))}  
  else if(stat=="p") {return(hmltm.p(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))}  
  else if(stat=="invp") {return(1/hmltm.p(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))}
  else if(stat=="esw") {return(hmltm.esw(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))}
  else if(stat=="invesw") {return(1/hmltm.esw(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))}
  else stop("Invalid stat")
}


#' @title Calculates E[p(x)] from model.
#'
#' @description
#' Calculates expected p(x) from model.
# 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' 
hmltm.p=function(pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                 type="response"){
  esw=hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type=type)
  return(esw/survey.pars$W)
}


#' @title Calculates 1/esw from model.
#'
#' @description
#' Calculates inverse of effective strip half-width from model.
# 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' 
hmltm.invp=function(pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                    type="response"){
  esw=hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type=type)
  return(survey.pars$W/esw)
}

#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip half-width from model.
#'
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param cov covariate matrix with each row corresponding to an observation. Must contain columns 
#' with variables appearing in \code{models}, and named accordingly, as well as column of perpendicular
#' distance, named "x". (Perpendicular distances only used if \code{to.x} is TRUE.)
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars availability hmm parameters, as for \code{\link{est.hmltm}}. Must have elements
#' \code{$Et} and \code{$Sigma.Et}
#' @param nx number of x-values (perpendicular distances) to use in evaluating esw.
#' @param type if "link", parameter vector \code{pars} is assumed to be on the link scale, else on 
#' the natural scale
#' @param to.x if TRUE integrates only out to \code{cov$x[i]} for observation i (else integrates to 
#' \code{survey.pars$W}).
#' 
#' @details
#' Returns effective strip half-width (esw) for fitted object hmmlt, integrating 
#' using Simpson's rule.
#'
hmltm.esw=function(pars,hfun,models,cov,survey.pars,hmm.pars,nx=100,type="response",to.x=FALSE){
  n=dim(cov)[1]
  if(to.x) {maxx=cov$x}
  else maxx=rep(survey.pars$W,n)
  p=rep(0,n)
  for(i in 1:n) {
    if(maxx[i]>0){
      xs=seq(0,maxx[i],length=nx) # set of poiints on which to evaluate p(see|x)
      px=hmltm.px(xs,pars,hfun,models,cov[i,],survey.pars,hmm.pars,type)
      p[i]=sintegral(px,xs)
    }
  }
  return(p)
}


#                          ---------------------------------------
#------------------------- End Derived Stats Calculation functions -----------------------


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
#' @param ks.plot If TRUE, does Q-Q plot. Point corresponding to largest difference between
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
    plot(1-e.cdf,cdf,xlab="Empirical Distribution Function",ylab="Cumulative Distribution Function",main="Forward Dist. Q-Q Plot",xlim=c(0,1),ylim=c(0,1),pch="+")
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
  n = sum(fplot$xy$seen)
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
#' @param ks.plot If TRUE, does Q-Q plot. Point corresponding to largest difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
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
    plot(edf,cdf,pch="+",xlim=c(0,1),ylim=c(0,1),xlab="Empirical Distribution Function",ylab="Cumulative Distribution Function",main="Perp. Dist. Q-Q Plot")
    lines(c(0,1),c(0,1))
    points(edf[worst],cdf[worst],col="red")
  }
  
  return(list(p.ks=1-p.ks,p.cvm=p.cvm,qq.x=edf,qq.y=cdf,x=hmmlt$xy$x[cdf.order]))
}




#-------------------------- End Goodness of Fit functions ------------------------
#                          ------------------------------



#                          -----------------------------------------------
#------------------------- Start Analytic Variance Estimation functions -----------------------



#------------------------- End Analytic Variance Estimation functions -----------------------
#                          ------------------------------------------


#                     ------------------------------------------------
#-------------------- Start Density and Abundance Estimation Functions --------------------

#' @title Calculates E[1/p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of inverse of fitted values, 1/p(x) for given 
#' observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
#' @param W actual half-width over which to integrate (overrides W in \code{hmmlt}).
#' 
#' @details
#' Identical to fitted.invp but returns data frame instead of numerical scalar or vector.
#' This is to allow it to be used in NDest for estimating density and abundance. (Also has extra 
#' parameter: W)
fitted.invp1=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100,W=NULL){
  if(is.null(W)) W=hmmlt$fitpars$survey.pars$W
  esw=fitted.esw1(hmmlt,obs,nx,W=W)
  return(data.frame(stratum=esw$stratum,transect=esw$transect,object=esw$object,invp=W/esw$esw))
}


#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip width from fitted model.
# 
#' @param hmmlt output from fit.hmltm()
#' @param obs observations (row numbers of hmmlt$xy) for which to calculate esw.
#' @param nx number of x values to use to implement Simpson's rule in perp dist dimension.
#' @param to.x If TRUE integrates only to observed x, else integrates to W.
#' @param all If TRUE then returns esw for every observation, else returns only that for
#'        first obs if there are no covariates; always returns esw for every observation
#'        if there are covariates.
#' @param W limit of perp. dist integration. If NULL, uses survey.pars$W.
#' 
#' @details
#' Designed to be called by \code{\link{fitted.invp1}}.
#' Identical to fitted.esw but returns list instead of numerical scalar or vector. This is to allow 
#' it to be used in \code{\link{NDest}} for estimating density and abundance.
#
#' Calls \code{\link{hmltm.esw}} to calclate effective stript width (esw) for 1 observer for fitted 
#' object \code{hmmlt}.
fitted.esw1=function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100,to.x=FALSE,all=FALSE,W=NULL){
  if(!is.null(obs)){
    if(max(obs)>dim(hmmlt$xy)[1]) stop("obs greater than number observations in hmmlt$xy")
    if(min(obs)<1) stop("obs < 1")
    cov=hmmlt$xy[obs,]
  } 
  if(is.nullmodel(hmmlt$models) & !all){
    cov=hmmlt$xy[1,]
  } else {
    cov=hmmlt$xy[obs,]
  }
  pars=hmmlt$fit$par
  hfun=hmmlt$h.fun
  models=hmmlt$models
  survey.pars=hmmlt$fitpars$survey.pars
  hmm.pars=hmmlt$fitpars$hmm.pars
  sID=list(stratum=hmmlt$xy$stratum,transect=hmmlt$xy$transect,object=hmmlt$xy$object) # unique sighting identifier
  esw=hmltm.esw1(pars,hfun,models,cov,survey.pars,hmm.pars,ID=sID,nx,to.x,type="response",W=W)
  return(esw)
}


#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip width from fitted model.
#' 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param ID list with elements \code{$stratum}, \code{$transect}, \code{$object} - to uniquely 
#' identify detections.
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' @param to.x If TRUE integrates only to observed x, else integrates to W
#' @param W limit of perpendicular dist integration. If NULL, uses \code{survey.pars$W}
#' 
#' @details
#' Designed to be called by \code{\link{fitted.esw1}}.
#' Identical to hmltm.esw but returns list instead of numerical scalar or vector, and allows limit 
#' of integration (W) to be specified explicitly, which overrides limit \code{survey.pars$W}.
#' This is to allow it to be used in \code{\link{NDest}} for estimating density and abundance.
#
#' Returns effective stript half-width (esw) for 1 observer for fitted object \code{hmmlt}, 
#' integrating using Simpson's rule.
#' 
hmltm.esw1=function(pars,hfun,models,cov,survey.pars,hmm.pars,ID,nx=100,type="response",to.x=FALSE,
                    W=NULL){
  nmax=dim(cov)[1]
  if(is.null(models$y) & is.null(models$x)) smax=length(ID$object) else smax=nmax # number detectoins
  if(to.x) {maxx=cov$x}
  else {
    maxx=rep(survey.pars$W,nmax)
    if(!is.null(W)) maxx=rep(W,nmax)
  }
  ustrat=utrans=uobject=esw=rep(0,smax) # vectors for unique stratum, transect, sighting numbers and esws
  for(i in 1:nmax) {
    if(maxx[i]>0){
      xs=seq(0,maxx[i],length=nx) # set of poiints on which to evaluate p(see|x)
      px=hmltm.px(xs,pars,hfun,models,cov[i,],survey.pars,hmm.pars,type)
      ustrat[i]=ID$stratum[i];utrans[i]=ID$transect[i];uobject[i]=ID$object[i] # record sighting ID
      esw[i]=sintegral(px,xs)
    }
  }
  if(smax>nmax){ # repeat single esw for all smax detections:
    for(i in 2:smax) {
      ustrat[i]=ID$stratum[i];utrans[i]=ID$transect[i];uobject[i]=ID$object[i] # record sighting ID
      esw[i]=esw[1]      
    }
  }
  return(list(stratum=ustrat,transect=utrans,object=uobject,esw=esw))
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



#' @title Estimates density and abundance.
#'
#' @description
#' Horvitz-Thompson like estimation of density and abundance of groups and of individuals, as well as
#' of group size (estimated as the ratio of individual density and group density estimates). Produces 
#' estimates by stratum and over all strata.
#' 
#' @param dat MRDS data frame. Must have cols "stratum","area","transect","L","size","object","x","y" 
#' (and possibly others).
#' @param hmltm.fit output from \code{\link{fit.hmltm}}.
#' @param W perpendicular truncation distance for estimation.
NDest=function(dat,hmltm.fit,W){
  maxx=max(na.omit(dat$x))
  if(maxx>W) {
    cat("Maximum perp. dist=",maxx,"is greater than W=",W,"\n")
    stop("You need a bigger W.")
  }
  # Add 1/p column
  dat$invp=rep(NA,dim(dat)[1])
  invp=fitted.invp1(hmltm.fit,W=W)
  for(i in 1:length(invp$object)) {
    row=which(dat$stratum==invp$stratum[i] & dat$transect==invp$transect[i] & dat$object==invp$object[i])
    if(length(row)>1) {
      cat("Target stratum:",invp$stratum[i],"\n")
      cat("Target transect:",invp$transect[i],"\n")
      cat("Target sighting:",invp$object[i],"\n")
      cat("Found >1: at rows",row,"\n")
      stop("")
    }
    dat$invp[row]=invp$invp[i]
  }
  # Calculate density and abundance by stratum
  m2km=1/1000
  strat=unique(dat$stratum)
  nstrat=length(strat)
  n=L=a=A=Dg=D=Ng=N=sbar=rep(0,nstrat+1)
  stratname=rep("",nstrat+1)
  for(i in 1:nstrat){
    stratname[i]=as.character(strat[i])
    vdat=dat[dat$stratum==strat[i],]
    trans=unique(vdat$transect)
    L.tr=0
    for(tr in 1:length(trans)) {
      L.tr=L.tr+vdat$L[min(which(vdat$transect==trans[tr]))]
    }
    L[i]=L.tr
    a[i]=L[i]*2*W*m2km
    A[i]=vdat$area[1]
    svdat=vdat[!is.na(vdat$object),]
    n[i]=length(svdat$invp)
    Dg[i]=sum(svdat$invp)/a[i]
    D[i]=sum(svdat$size*svdat$invp)/a[i]
    sbar[i]=D[i]/Dg[i]
    Ng[i]=Dg[i]*A[i]
    N[i]=D[i]*A[i]
  }
  stratname[nstrat+1]="Total"
  Ng[nstrat+1]=sum(Ng[1:nstrat])
  N[nstrat+1]=sum(N[1:nstrat])
  A[nstrat+1]=sum(A[1:nstrat])
  Dg[nstrat+1]=Ng[nstrat+1]/sum(A[1:nstrat])
  D[nstrat+1]=N[nstrat+1]/sum(A[1:nstrat])
  n[nstrat+1]=sum(n[1:nstrat])
  L[nstrat+1]=sum(L[1:nstrat])
  a[nstrat+1]=sum(a[1:nstrat])
  sbar[nstrat+1]=D[nstrat+1]/Dg[nstrat+1]
  # add transect frequency:
  tfreq=apply(table(dat$stratum,dat$transect)>0,1,sum)
  k=c(tfreq,sum(tfreq))
  
  return(list(invp=invp,ests=data.frame(stratum=stratname,n=n,k=k,L=L,covered.area=a,stratum.Area=A,Dgroups=signif(Dg,3),Ngroups=signif(Ng,3),
                                        mean.size=round(sbar,1),D=signif(D,5),N=round(N,1))))
  #                                        mean.size=signif(sbar,3),D=signif(D,3),N=signif(N,3))))
}




#' @title Line transect estimation with a hidden Markov availability model.
#'
#' @description
#' \code{est.hmltm} estimates group and individual density and abundance, together with mean 
#' group size, by stratum, from (1) line transect data that includes forward detection distances and 
#' (2) estimated Markov model or hidden Markov model availability prameters. 
#'
#' @param dat data frame in distance-like format, but including forward distances of detections. The 
#' following are compulsory elements (field name in quotes, contents in brackets): 
#' "stratum" (survey stratum: must be numeric), "area" (stratum area), "transect" (transect number: must
#' be numeric), "L" (transect length), "size" (group size), "object" (unique detection identifier: must
#' be numeric), "x" (perpendicular distance), "y" (forward distance).
#' @param pars starting parameter values.
#' @param FUN detection hazard functional form name (character). Currently implemented forms are 
#' "h.IP.0", "h.EP1.0", "h.EP2.0", "h.EP1x.0", "h.EP2x.0". (See Vignette "Specifying models and parameter 
#' starting values" for details.)
#' @param models list of characters with elements \code{$y} and \code{$x} specifying models for the y- 
#' and x-dimension detection hazard scale parameters. Must be either \code{NULL} or regression model 
#' format (without response on left, e.g. "~size").
#' @param survey.pars a list containing the following elements (in any order):
#'  \itemize{
#'  \item {$spd} {speed of observer,}
#'  \item {$W} {perpendicular distance right-truncation point,}
#'  \item {$Wl} {perpendicular distance left-truncation point,}
#'  \item {$ymax} {forward distance by which detection probability is effectively zero,}
#'  \item {$dT} {availability process (Markov chain) time step size.}
#' }
#' @param hmm.pars a list containing the parameters of animals' availability processs hidden Markov 
#' model (HMM), as follows (in any order):
#' \itemize{
#'  \item {$Pi} {a 2x2xm HMM transition probability matrix, where m is the number of availability HMMs 
#'  being used to model animal availability. If m>1, each set of HMM parameters is treated as a random
#'  sample from the set of HMM parameters in the population.}
#'  \item {$pm} {a2xm matrix of HMM state-dependent Bernoulli distribution parameters (the probabilities
#'  of being available, given the animal's "behavioural" state - i.e. the state of the hidden Markov 
#'  chain)}
#'  \item {$delta} {a 2xm matrix of stationary distribution of a Markov chain, the ith of which has
#'  transition probability matrix Pi[,,i].}
#' }
#' And if the HMM was constructed from mean times animals are available and unavailable (by means of
#' function \code{\link{make.hmm.pars.from.Et}} for example), then also
#' \itemize{
#'  \item {$Et} {a 2xm matrix in which the first element is the mean time animals are UNavailable
#'  in a single available-unavailable cycle, and the second element is the corresponding mean time that
#'  they are available,}
#'  \item {Sigma.Et} {a 2x2xm matrix, in which Sigma.Et[,,i] is the variance-covariance matrix of 
#'  Et[,i.] (i.e. the variance-covariance matrix of Et for the ith availability model).}
#' }
#' @param control.fit list with elements
#' \itemize{
#'  \item{$hessian} {logical) - if TRUE Hessian is estimated and returned, else not,}
#'  \item{$nx} {(scalar) - the number of intervals to use with Simpson's rule integration over y. 
#'    \code{nx=64} seems safe; smaller number makes computing faster.}
#' }
#' @param control.opt as required by \code{\link{optim}} (and hence by \code{\link{fit.hmltm}}).
#' @param twosit TRUE if \code{dat} is in mrds format (with two lines per detection), else assumes
#'              that \code{dat} is in cds format (with one line per detection).
#' @param notrunc if TRUE, does not do any perp dist truncation, else uses \code{survey.pars$W}
#' and \code{$Wl} to do perp dist truncation.
#' @param W.est right truncation perpendicular distance for estimation. Can't be less than maximum 
#' perpendicular distance (x) in the line transect data frame \code{dat}, but can be less than 
#' the max perpendicular distance used for fitting (\code{survey.pars$W}).
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
#'
#' @return A list with four elements: \code{hmltm.fit}, \code{point}, \code{dat}, \code{W.est}. 
#' Their contents are as follows:
#' 
#' \code{hmltm.fit} is the output from \code{fit.hmltm}, i.e. a list containing the following elements:
#' \itemize{
#'  \item{xy} {dat used in fitting (input reflection).}
#'  \item{phats} {estimated detection probabilities of all detections.}
#'  \item{phat} {1/mean(1/phat).}
#'  \item{pzero} {estimated detection probabilities at perpendicular distance.}
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
#' \code{point} is a list containing two elements:
#' \itemize{
#'  \item{invp} {is a data frame containing one row for every observation, with the first three columns
#'   giving the stratum, transect and object identifier for the observation, and the final column 
#'   (invp) giving the estimate of the inverse of the probability of detection for the observation.}
#'  \item{ests} {is a data frame with one row per stratum and a final row for all strata combined,
#'  and columns giving the number of detections in the stratum (n), the line lingth in the stratum 
#'  (L), the covered area in the stratum (covered.area=2WL), the stratum area (stratum.Area), the 
#'  estimated group density in the stratum (Dgroups), the estimated group abunance in the stratum 
#'  (Ngroups), the estimated mean group size in the stratum (mean.size), the individual denstiy
#'  in the stratum (D), and the abundance in the stratum (N).}
#' }
#' \code{dat} is the data frame passed to \code{est.hmltm}.
#' 
#' \code{W.est} is the right perpendicular distance used for estimation (and passed to 
#' \code{est.hmltm}.)
#' 
#' @references Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 
#' 2013. Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
est.hmltm=function(dat,
                   pars,FUN,models=list(y=NULL,x=NULL),
                   survey.pars,hmm.pars,control.fit,control.opt,
                   twosit=FALSE,notrunc=FALSE,W.est=NULL,groupfromy=NULL){
  # convert mrds data to cds data format (combining detections)
  if(twosit) dat1=make.onesit(dat) else dat1=dat
  # data truncation:
  if(!notrunc) if(min(na.omit(dat1$x)) < survey.pars$Wl | survey.pars$W < max(na.omit(dat1$x))) 
    dat1=truncdat(dat1,minval=survey.pars$Wl,maxval=survey.pars$W,twosit=FALSE)
  # extract sightings rows only
  #  sits1=!is.na(dat1$seen)
  #  srows1=sits1 & dat1$seen==1 # mark rows that have sightings
  srows1=!is.na(dat1$object)
  sdat1=dat1[srows1,]
  # fit the model:
  hmltm.fit=fit.hmltm(sdat1,pars,FUN,models,survey.pars,hmm.pars,control.fit,control.opt,groupfromy=groupfromy)
  # estimate density, etc:
  if(is.null(W.est)) W.est=survey.pars$W
  if(!is.null(survey.pars$Wl)) W.est=survey.pars$W-survey.pars$Wl # since in this case data all shifted left by $Wl
  point=NDest(dat1,hmltm.fit,W.est)
  hmltm.obj=list(hmltm.fit=hmltm.fit,point=point,dat=dat1,W.est=W.est)
  class(hmltm.obj)=c(class(hmltm.obj),"hmltm")
  return(hmltm.obj)
}

#----------------------- End Density and Abundance Estimation Functions ----------------------
#                        ----------------------------------------------


#                          ------------------------------
#-------------------------- Start Bootstrap functions ------------------------

#' @title Summarise bootstrap results.
#'
#' @description
#' Produces summarry of bootstrap results
#' 
#' @param bests output from \code{\link{bs.hmltm}}.
#' @param ests hmltm object with points estimates.
#' @param cilevel confidence level for confidence intervals by percentile method. 
#' @param write.csvs if TRUE, writes each output (see Value below) to separate .csv file.
#' @param dir directory to which to write outputs if \code{write.csvs} is TRUE. 
#'
#' @return If \code{ests} is NULL:\cr
#' Returns a list with elements as follows for each statistic in \code{bests}:
#' \itemize{
#'  \item{nbad} {number of bad estimates of the statistic in question that were excluded from the 
#'  summary. (Bad estimates occur, for example, when a bootstrap sample involves no detections and
#'  the estimation function tries to calculate mean group size.)}
#'  \item{means} {Bootstrap means of the statistic in question.}
#'  \item{cv} {Percentage CVs of the statistic in question.}
#'  \item{se} {SEs of the statistic in question.}
#'  \item{lower} {Lower \code{cilevel} percentiles of the statistic in question.}
#'  \item{upper} {Upper \code{cilevel} percentiles of the statistic in question.}
#'} \cr
#' If \code{ests} is an object of class 'hmltm':\cr
#' Returns a data frame with rows only for every stratum with detections (and the total) and columns as follows: \cr
#' \itemize{
#' \item{Stratum} stratum number \cr
#' \item{n} original number of detections \cr
#' \item{n.L} original encounter rate \cr
#' \item{CV.n.L} percentage coefficient of variation of encounter rate\cr
#' \item{N.grp} original estimate of group abundance\cr
#' \item{CV.N.grp} percentage coefficient of variation of group abundance\cr
#' \item{N.grp.lo} lower bound of group abundance confidence interval\cr
#' \item{N.grp.hi} upper bound of group abundance confidence interval\cr
#' \item{Es} original mean group size estimate\cr
#' \item{CV.Es} percentage coefficient of variation of mean group size\cr
#' \item{Es.lo} lower bound of mean group size confidence interval\cr
#' \item{Es.hi} upper bound of mean group size confidence interval\cr
#' \item{N} original individual abundance estimate\cr
#' \item{CV.N} percentage coefficient of individual abundance\cr
#' \item{N.lo} lower bound of individual abundance confidence interval\cr
#' \item{N.hi} upper bound of individual abundance confidence interval\cr
#'}
bootsum=function(bests,ests=NULL,cilevel=0.95,write.csvs=FALSE,dir=getwd()){
  if(cilevel<=0 | cilevel>=1) stop("cilevel must be greater than 0 and less than 1.")
  if(!is.null(ests) & !inherits(ests,"hmltm")) stop("ests must be of class 'hmltm'.")
  ns=apply(bests[,2,],1,sum) # sum sample sizes by stratum
  keepstrat=(ns>0) # only consider strata with some detections
  if(!is.null(ests)) keepcols=c("stratum","n","L","Ngroups","mean.size","N")
  else keepcols=colnames(bests)
  if(!is.null(ests)) bsests=bests[keepstrat,keepcols,]
  else bsests=bests[keepstrat,,]
  # replace L with encounter rate
  Lcol=which(colnames(bests[,,1])=="L")
  bsests[,Lcol,]=bsests[,"n",]/bsests[,"L",] 
  colnames(bsests)[Lcol]="n/L"
  # replace stratum.Area with p
  Acol=which(colnames(bests[,,1])=="stratum.Area")
  bsests[,Acol,]=bsests[,"n",]/(bsests[,"covered.area",]*bsests[,"Dgroups",]) 
  colnames(bsests)[Acol]="p"
  bdim=dim(bsests)
  nstrat=bdim[1]
  nests=bdim[2]
  cv=matrix(rep(NA,nstrat*nests),ncol=bdim[2])
  rownames=dimnames(bsests)[[1]]
  rownames[length(rownames)]="Total"
  colnames=dimnames(bsests)[[2]]
  dimnames(cv)=list(rep("",length(rownames)),colnames)
  nbad=se=means=lower=upper=cv
  B=bdim[3]
  cat("Results from ",B," bootstrap replicates:\n",sep="")
  cat("----------------------------------------\n")
  for(i in 1:nstrat) {
    for(j in 1:nests){
      goodests=na.omit(bsests[i,j,])
      nbad[i,j]=B-length(goodests)
      means[i,j]=mean(goodests)
      se[i,j]=sd(goodests)
      cv[i,j]=se[i,j]/means[i,j]
      perc=quantile(goodests,probs=c((1-cilevel)/2,(1-(1-cilevel)/2)))
      lower[i,j]=perc[1]
      upper[i,j]=perc[2]
    }
  }
  # remove stats of stratum name, and make stratum a character
  nbad=as.data.frame(nbad,row.names=1:dim(nbad)[1])
  means=as.data.frame(means,row.names=1:dim(nbad)[1])
  se=as.data.frame(se,row.names=1:dim(nbad)[1])
  cv=as.data.frame(cv,row.names=1:dim(nbad)[1])
  lower=as.data.frame(lower,row.names=1:dim(nbad)[1])
  upper=as.data.frame(upper,row.names=1:dim(nbad)[1])
  nbad[,1]=rownames
  means[,1]=rownames
  se[,1]=rownames
  cv[,1]=rownames
  lower[,1]=rownames
  upper[,1]=rownames
  
  outsum=list(nbad=nbad,mean=means,se=se,cv=cv,lower=lower,upper=upper)
  if(!is.null(ests)) outsum=strat.estable(ests$point$ests[keepstrat,],outsum$cv)
  
  if(write.csvs){
    if(is.null(ests)) {
      write.csv(nbad,file=paste(dir,"nbad.csv",sep=""),row.names=FALSE)
      write.csv(means,file=paste(dir,"means.csv",sep=""),row.names=FALSE)
      write.csv(se,file=paste(dir,"se.csv",sep=""),row.names=FALSE)
      write.csv(cv,file=paste(dir,"cv.csv",sep=""),row.names=FALSE)
      write.csv(lower,file=paste(dir,"lower.csv",sep=""),row.names=FALSE)
      write.csv(upper,file=paste(dir,"upper.csv",sep=""),row.names=FALSE)
    } else write.csv(outsum,file=paste(dir,"bs-summary.csv",sep=""),row.names=FALSE)
  }
  
  return(outsum)
}


#' @title Tabulate bootstrap results.
#'
#' @description
#' Produces brief summarry of bootstrap results
#' 
#' @param est output from \code{\link{est.hmltm}}.
#' @param cv \code{$cv} element of output from \code{\link{bs.hmltm}}. 
#' 
#' @return Returns a data frame with columns as follows: \cr
#' Stratum: stratum number \cr
#' n: original number of detections \cr
#' n.L: original encounter rate \cr
#' CV.n.L: percentage coefficient of variation of encounter rate\cr
#' N.grp original estimate of group abundance\cr
#' CV.N.grp: percentage coefficient of variation of group abundance\cr
#' N.grp.lo: lower bound of group abundance confidence interval\cr
#' N.grp.hi: upper bound of group abundance confidence interval\cr
#' E.s: original mean group size estimate\cr
#' CV.E.s: percentage coefficient of variation of mean group size\cr
#' N: original individual abundance estimate\cr
#' CV.N: percentage coefficient of individual abundance\cr
#' N.lo: lower bound of individual abundance confidence interval\cr
#' N.hi: upper bound of individual abundance confidence interval\cr
strat.estable=function(est,cv){
  min=est[,"n"]
  N=est[,"N"]
  Ngrp=est[,"Ngroups"]
  Es=est[,"mean.size"]
  cv.N=cv[,"N"]
  cv.Ngp=cv[,"Ngroups"]
  cv.Es=cv[,"mean.size"]
  N.ci=lnci.nmin(min,N,cv.N)
  Ngrp.ci=lnci.nmin(min,Ngrp,cv.N)
  Es.ci=lnci.nmin(0,Es,cv.Es)  
  outp=data.frame(Stratum=est$stratum,
                  n=est$n,
                  n.L=signif(est$n/est$L,3),
                  CV.n.L=round(cv[,3]*100),
                  N.grp=round(est$Ngroups),
                  CV.N.grp=round(cv$Ngroups*100),
                  N.grp.lo=round(Ngrp.ci$lower),
                  N.grp.hi=round(Ngrp.ci$upper),
                  Es=round(est$mean.size,1),
                  CV.Es=round(cv.Es*100),
                  Es.lo=round(Es.ci$lower),
                  Es.hi=round(Es.ci$upper),
                  N=round(est$N),
                  CV.N=round(cv$N*100),
                  N.lo=round(N.ci$lower),
                  N.hi=round(N.ci$upper)
  )
  return(outp)
}

#' @title Calculate CV of mean time available.
#'
#' @description
#' Calculates CV of mean time available and unavailable from an hmm.pars object. 
#'
#' @param hmm.pars object of class \code{hmm.pars}. 
#' @param B number of bootstrap replicates to do.
#' 
cv.avail=function(hhm.pars,B=1000){
  n=dim(hmm.pars$Et)[2]
  b.Et=array(rep(NA,B*2*n),dim=c(B,2,n),dimnames=list(1:B,State=c("Unavailable","Available"),Animal=1:n))
  cv.a=cv.u=rep(NA,n)
  for(i in 1:n){
    lN=Npars.from.lNpars(hmm.pars$Et[,i],hmm.pars$Sigma.Et[,,i]) # get normal parameters corresponding to lognormal(mu,Sigma)
    b.Et[,,i]=exp(mvrnorm(B,lN$mu,lN$Sigma)) # resample availability parameters on lognormal scale
    a=b.Et[,1,i]/(b.Et[,1,i]+b.Et[,2,i])
    u=b.Et[,2,i]/(b.Et[,1,i]+b.Et[,2,i])
    cv.a[i]=cv(a)
    cv.u[i]=cv(u)
  }
  return(list(cv.a=cv.a,cv.u=cv.u))
}

#' @title Bootstrap for hmltm model.
#'
#' @description
#' Stratified nonparameteric bootstrap of hidden Markov line transect model (hmltm) with transects 
#' as the sampling units. 
#'
#' @param hmltm.est output from \code{\link{est.hmltm}}. 
#' @param B number of bootstrap replicates to do.
#' @param hmm.pars.bs output \code{\link{hmmpars.boot}}, containing sets of refitted HMM parameter values.
#' @param bs.trace amount of reporting that \code{\link{optim}} should do while fitting models to
#' bootstrapped datasets. (\code{bs.trace}=0 is no reporting, which is fasters. See component 
#' \code{trace} of parameter \code{control.opt} of function \code{\link{optim}} for more details.)
#' @param report.by frequency with which to report count of number of bootstraps completed.
#' @param fixed.avail whether to treat the hidden Markov availability model as fixed 
#' (\code{fixed.avail}=TRUE) or random (\code{fixed.avail}=FALSE). If \code{fixed.avail}=FALSE, 
#' the availability model is also bootstrapped (see Details below).
#' 
#' @seealso \code{\link{bootsum}} summarises output from this function.
#' 
#' @details
#' If \code{fixed.avail}=FALSE, then: (1) IF \code{hmltm.est$hmltm.fit$fitpars$hmm.pars$Et} is not 
#' \code{NULL}, the availability process is bootstrapped by drawing pairs of mean times available 
#' (Ea) and unavailable (Eu) from a bivariate lognormal distribution with mean 
#' \code{hmltm.est$hmltm.fit$fitpars$hmm.pars$Et} and standard deviation 
#' \code{hmltm.est$hmltm.fit$fitpars$hmm.pars$Sigma.Et} and converting this to Markov Model transition
#' probability matrix parameters using \code{\link{makePi}}, ELSE IF \code{hmm.pars.bs} is not 
#' \code{NULL}, random samples with replacement, of HMM parameter sets are taken from 
#' \code{hmm.pars.bs}.
#' 
#' @seealso \code{\link{hmmpars.boot}}, which is is used to generate \code{hmm.pars.bs}.
#' 
bs.hmltm=function(hmltm.est,B,hmm.pars.bs=NULL,bs.trace=0,report.by=10,fixed.avail=FALSE){
  #------------------------------------------------------------------------------------------
  # Produces 3-dim array containing B sets of density and abundance estimates from
  #------------------------------------------------------------------------------------------
  # data extraction:
  dat1=hmltm.est$dat # extract original data frame from fitted object
  W.est=hmltm.est$W.est # extract estimation perp dist truncation
  survey.pars=hmltm.est$hmltm.fit$fitpars$survey.pars
  hmm.pars=hmltm.est$hmltm.fit$fitpars$hmm.pars
  control.fit=hmltm.est$hmltm.fit$fitpars$control.fit
  control.opt=hmltm.est$hmltm.fit$fitpars$control.opt
  control.opt$trace=bs.trace
  h.fun=hmltm.est$hmltm.fit$h.fun
  models=hmltm.est$hmltm.fit$models
  pars=hmltm.est$hmltm.fit$fit$par
  
  if(!fixed.avail) {
    # bootstrap availability parameters
    if(!is.null(hmm.pars$Et)){
      cat("Bootstrap with parametric resampling of mean times available and unavailable.\n")
      flush.console()
      n=dim(hmm.pars$Et)[2]
      ns=bsample(1:n,n,replace=TRUE) # sample animals to use
      b.Et=array(rep(NA,B*2*n),dim=c(B,2,n),dimnames=list(1:B,State=c("Unavailable","Available"),Animal=1:n))
      for(i in ns){
        lN=Npars.from.lNpars(hmm.pars$Et[,i],hmm.pars$Sigma.Et[,,i]) # get normal parameters corresponding to lognormal(mu,Sigma)
        b.Et[,,i]=exp(mvrnorm(B,lN$mu,lN$Sigma)) # resample availability parameters on lognormal scale
      }
    } else if(!is.null(hmm.pars.bs)) { # resample availability parameters
      cat("Bootstrap with parametric resampling of HMM parameters.\n")
      flush.console()
      nhmm=dim(hmm.pars.bs)[1]
      if(nhmm==1) stop("Only one set of hmm pars. need multiple sets of pars if not fixed avail. (i.e. if hmm.pars.bs!=NULL)")
      reps=bsample(1:nhmm,B,replace=TRUE)
    } else {
      warning("No availability parameters to resample from. Treating availability parameters as constant.")
      flush.console()
      fixed.avail=TRUE
    }
  } else {
    cat("Bootstrap treating HMM parameters as constant.\n")
    flush.console()
  }
  
  # get stuff needed to set up bootstrap datasets
  strat=unique(dat1$stratum) # unique stratum names
  nrows=table(dat1$stratum,dat1$transect) # number of rows required for each transect
  get.ntrans=function(trans) sum(trans>0)
  ntrans=apply(nrows,1,get.ntrans) # number of transects in each stratum
  nstrat=length(ntrans) # number of strata
  transects=as.numeric(colnames(nrows)) # vector of transect numbers
  
  # create bootstrap data, estimate, and store estimates
  estdim=dim(hmltm.est$point$ests) #### bitchange
  bestdim=c(estdim,B) #### bitchange
  bestdimnames=list(as.character(hmltm.est$point$ests$stratum),colnames(hmltm.est$point$ests),1:B)
  best=array(rep(NA,B*prod(estdim)),dim=bestdim,dimnames=bestdimnames)
  bn=rep(0,B) #### bitchange
  for(b in 1:B) { # do bootstraps
    # get hmm pars
    if(!fixed.avail) {
      if(!is.null(hmm.pars$Et)) {
        Pi=makePi(b.Et[b,1,],b.Et[b,2,]) # make Pi from resampled times
        delta=matrix(rep(NA,2*n),nrow=2)
        for(i in 1:n) delta[,i]=compdelta(Pi[,,i])
        b.hmm.pars=hmm.pars
        b.hmm.pars$pm=hmm.pars$pm
        b.hmm.pars$Pi=Pi
        b.hmm.pars$delta=delta
        b.hmm.pars$Et=b.Et[b,,]
      } else {
        b.hmm.pars=unvectorize.hmmpars(hmm.pars.bs[reps[b],]) # resampled HMM pars
      }
    } else {
      b.hmm.pars=hmm.pars # fixed availability pars
    }
    if(is.element("pm",names(b.hmm.pars))) names(b.hmm.pars)[which(names(b.hmm.pars)=="pm")]="pm" # to fix naming cock-up when creating hmmpars.bs
    # resample transects with replacement
    newtransind=matrix(rep(NA,nstrat*max(ntrans)),nrow=nstrat) # matrix of indices for resampled transects
    for(st in 1:nstrat) newtransind[st,1:ntrans[st]]=bsample(which(nrows[st,]>0),ntrans[st],replace=TRUE) # indices of matrix nrows for resampled transects
    newnrows=0;for(st in 1:nstrat) newnrows=newnrows+sum(na.omit(nrows[st,newtransind[st,]])) # calc number rows for bootstrap data frame
    bdat1=data.frame(matrix(rep(NA,dim(dat1)[2]*newnrows),nrow=newnrows)) # set up empty bootstrap data frame
    names(bdat1)=names(dat1)
    start=1;tno=1 # fill in new data frame from old:
    for(st in 1:nstrat) for(tr in 1:ntrans[st]) {
      addtrans=transects[newtransind[st,tr]]
      nadd=nrows[st,newtransind[st,tr]]
      newi=start:(start+nadd-1)
      oldi=which(dat1$stratum==strat[st] & dat1$transect==addtrans)
      bdat1[newi,]=dat1[oldi,]
      bdat1$transect[newi]=tno
      start=start+nadd
      tno=tno+1
    }
    bn[b]=length(bdat1$object[!is.na(bdat1$object)])
    if(b==1) cat("Sample sizes: ")
    if(b%%report.by==0) {
      cat(paste(bn[b],";",sep=""),"Iterations done:",b,"\n")
      if(b!=B) cat("Sample sizes: ")
    }
    else cat(paste(bn[b],";",sep=""))
    bhat=est.hmltm(bdat1,pars=pars,FUN=h.fun,models=models,survey.pars,b.hmm.pars,control.fit,control.opt,notrunc=TRUE,W.est=W.est) 
    for(st in 1:(nstrat+1)) best[st,,b]=as.numeric(bhat$point$ests[st,])
  }
  class(best)=c(class(best),"hmltm.bs")
  return(best)
}

#' @title Normal from logNormal parameters.
#'
#' @description
#' Returns mean and covariance matrix of multivariate normal random variables 
#' which, when logged generate lognormal random variables with mean mu and 
#' covariance matrix Sigma.
#' 
#' @param mu multivariate logNormal distribution mean.
#' @param Sigma multivariate logNormal distribution variace-covarice matrix.
#' 
Npars.from.lNpars=function(mu,Sigma){
  cv2=diag(Sigma)/mu^2 # Sigma is variance matrix
  logvar=log(cv2+1)
  logmu=log(mu)-logvar/2
  logSigma=Sigma*0
  np=length(mu)
  for(i in 1:(np-1)) {
    logSigma[i,i]=logvar[i] # variance of normal
    for(j in (i+1):np) {
      logSigma[i,j]=log(Sigma[i,j]/(exp(logmu[i]+logmu[j]+(logvar[i]+logvar[j])/2))+1) # covariance of normal
      logSigma[j,i]=logSigma[i,j]
    } 
  }
  logSigma[np,np]=logvar[np] # variance of normal
  return(list(mu=logmu,Sigma=logSigma))
}


#' @title Resample availability HMM.
#'
#' @description
#' Resamples a hidden Markov model (HMM) data from multiple observed availability time series 
#' (one per animal).
#' 
#' @param availhmm list with availability HMM paramters, as for the \code{hmm.pars} parameter of 
#' \code{\link{est.hmltm}}.
#' @param adat list of availability data time series. The ith element of the list must be named 
#' \code{$ai} and must be a vector of 0s and 1s, comprising the ith animal's time series, with 0 
#' corresponding to being unavailable and 1 to being available.
#' @param animals a vector of integers indicating which of the elements of \code{adat} are to be
#' resampled (e.g. animals=c(1,3,4) means the time series of animals 1, 3 and 4 only will be used).
#' @param seed random number seed (integer).
#' @param nperow number of iterations printed on one row if printprog=TRUE.
#' @param printprog if TRUE, prints progress through animals as it resamples.
#'
#' @details
#' Simulates a new series of availability observations (0s and 1s) using functions \code{dthmm} 
#' and \code{simulate} from library \code{\link{HiddenMarkov}}, then fits a HMM to these data
#' using function \code{BaumWelch} from library \code{\link{HiddenMarkov}}. Constructs a new 
#' hmm.pars object from the fitted HMM parameters.
#' 
#' @return A list with the same format as the input \code{availhmm}
#' 
#' @examples
#' data(bowhead.hmm.pars)
#' data(bowhead.adat)
#' animals=1:8
#' resamp=resample.hmmpars(bowhead.hmm.pars,bowhead.adat,animals)
#' 
resample.hmmpars=function(availhmm,adat,animals,seed=NULL,nperow=10,printprog=TRUE){
  b.availhmm=availhmm # initialise bootstrap HMM parameters
  # filter out animals not chosen
  b.availhmm$Pi=b.availhmm$Pi[,,animals]
  b.availhmm$pm=b.availhmm$pm[,animals]
  b.availhmm$delta=b.availhmm$delta[,animals]
  # simulate and fit new HMMs
  newanimal=1
  for(animal in animals) {
    hmmobj=dthmm(adat[[animal]],availhmm$Pi[,,animal],availhmm$delta[,animal],"binom",pm=list(prob=availhmm$pm[,animal]),pn=list(size=rep(1,length(adat[[animal]]))),discrete = TRUE)
    sdat=simulate(hmmobj, nsim=length(hmmobj$x),seed=seed)
    fitanimal=BaumWelch(sdat,control=list(maxiter=500,tol=1e-05,prt=FALSE,posdiff=TRUE,converge = expression(diff < tol)))
    b.availhmm$pm[,newanimal]=fitanimal$pm$prob
    b.availhmm$Pi[,,newanimal]=fitanimal$Pi
    b.availhmm$delta[,newanimal]=compdelta(fitanimal$Pi)
    if(printprog) {
      if(animal==1) cat("Individuals done:")
      if((animal%/%nperow)*nperow==animal) { # got an exact multiple of 10
        if(animal>1) cat(": Total=",animal,"\nIndividuals done:")
      }
    }
    cat(" *")
    newanimal=newanimal+1
  }
  cat("\n Total=",animal," individuals done\n")
  return(b.availhmm)  
}


#' @title Reformat hmm.pars object as a vector.
#'
#' @description
#' Reformats hmm.pars object as a vector - so that it can easily be written to file.
#' 
#' @param hmmpars list with availability HMM paramters, as for the \code{hmm.pars} parameter of 
#' \code{\link{est.hmltm}}.
#' 
#' @return A numeric vector.
vectorize.hmmpars=function(hmmpars) {
  Pi=hmmpars$Pi
  pm=hmmpars$pm
  delta=hmmpars$delta
  dims=dim(Pi)
  if(length(dims)!=3) stop("Pi must be a 3-D array")
  if(dims[1]!=dims[2]) stop("First two dimensions of Pi unequal")
  if(dim(pm)[1]!=dims[1] | dim(pm)[2]!=dims[3]) stop("dim(pm) inconsistent with dim(Pi)")
  if(dim(delta)[1]!=dims[1] | dim(delta)[2]!=dims[3]) stop("dim(delta) inconsistent with dim(Pi)")
  return(c(as.vector(dim(Pi)),as.vector(Pi),as.vector(pm),as.vector(delta)))
}

#' @title Reformat vector as hmm.pars object.
#'
#' @description
#' Reformats vector that was vectorised using \code{\link{vectorize.hmmpars}}, as a hmm.pars object.
#' 
#' @param hv vector that was vectorised using \code{\link{vectorize.hmmpars}}.
#' 
#' @return A hmm.pars object (a list).
unvectorize.hmmpars=function(hv) {
  m3d=hv[1:3]
  m2d=c(m3d[1],m3d[3])
  m3size=prod(m3d)
  m2size=prod(m2d)
  Pi=array(hv[(3+1):(3+m3size)],dim=m3d)
  pm=array(hv[(3+m3size+1):(3+m3size+m2size)],dim=m2d)
  delta=array(hv[(3+m3size+m2size+1):(3+m3size+m2size+m2size)],dim=m2d)
  return(list(pm=pm,Pi=Pi,delta=delta))
}

#' @title Bootstrap availability HMM.
#'
#' @description
#' Bootstraps hidden Markov model (HMM) data from multiple observed availability time series 
#' (one per animal).
#' 
#' @param availhmm list with availability HMM paramters, as for the \code{hmm.pars} parameter of 
#' \code{\link{est.hmltm}}.
#' @param adat list of availability data time series. The ith element of the list must be named 
#' \code{$ai} and must be a vector of 0s and 1s, with 0 corresponding to being unavailable and 1 to 
#' being available.
#' @param animals a vector of integers indicating which of the elements of \code{adat} are to be
#' resampled.
#' @param B number of bootstraps to perform (integer).
#' @param seed random number seed (integer).
#' @param printprog if TRUE prints progress through animals as it resamples.
#'
#' @details
#' Simulates a new series of availability observations (0s and 1s) using functions \code{dthmm} 
#' and \code{simulate} from library \code{\link{HiddenMarkov}}, then fits a HMM to these data
#' using function \code{BaumWelch} from library \code{\link{HiddenMarkov}}. Constructs a new 
#' hmm.pars object from the fitted HMM parameters. 
#' 
#' Does the above \code{B} times, each time reformtting the hmm.pars object as a vector using 
#' \code{\link{vectorize.hmmpars}} and then entering this as the next row in a matrix of dimension 
#' \code{B}xT, where T=\code{3+length(animals)*(nstate^2+nstate*2)} and \code{nstate} is the 
#' number of states in the HMM.
#' 
#' @return A matrix of dimension \code{B}xT, where T is as described above.
#' 
#' @examples
#' data(bowhead.hmm.pars)
#' data(bowhead.adat)
#' animals=1:8
#' bs=hmmpars.boot(bowhead.hmm.pars,bowhead.adat,animals,B=3)
#' 
hmmpars.boot=function(availhmm,adat,animals,seed=NULL,B,printprog=TRUE){
  # initialise matrix of correct dimensions:
  na=length(adat)
#  ncol=1
#  for(i in 1:na) if(length(adat[[i]])>ncol) ncol=length(adat[[i]])
  onelength=length(vectorize.hmmpars(list(Pi=availhmm$Pi[,,1,drop=FALSE],
              pm=availhmm$pm[,1,drop=FALSE],delta=availhmm$delta[,1,drop=FALSE])))
  ncol=3+(onelength-3)*length(animals) # 3 numbers specify vectorized length parameters
  hmmat=matrix(rep(NA,ncol*B),ncol=ncol)
  # do the bootstrapping
  for(b in 1:B) {
    hmp=resample.hmmpars(availhmm,adat,animals,seed=NULL)
    hmmvec=vectorize.hmmpars(hmp)
#    nobs=length(hmmvec)
#    hmmat[b,nobs]=hmmvec
    hmmat[b,]=hmmvec
    if(printprog) cat("Resamples done: ",b,"\n")
  }
  cat("\n Total=",b," Resamples done\n")
  return(hmmat)
}


#' @title Detection probability bootstrap with availability process times.
#'
#' @description
#' Nonparametric bootstrap of detection data with estimation of detection probabilities. If 
#' \code{fixed.avail}=FALSE, does parametric resampling of mean times available and unavailable for 
#' every resample of detection data, else treats these mean times as fixed.
#' 
#' @param dat detection data frame constructed by removing all rows with no detections from a 
#' data frame of the sort passed to \code{\link{est.hmltm}}.
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars availability hmm parameters, as for \code{\link{est.hmltm}}. Must have elements
#' \code{$Et} and \code{$Sigma.Et}
#' @param control.fit list controlling fit, as for \code{\link{est.hmltm}}.
#' @param control.opt list controlling function \code{\link{optim}}, as for \code{\link{est.hmltm}}.
#' @param fixed.avail if TRUE, hmm.pars is treated as fixed, else element \code{$Et} is parametrically 
#' resampled.
#' @param B number of bootstrap replicates.
#' 
#' @details
#' The rows of data frame \code{dat} are resampled with replacement to create new data frames with as 
#' many detections as were in \code{dat}. If \code{fixed.avail}=TRUE, then a pair of new mean times 
#' available and unavailable (\code{$Et}s) are generated for each resampled data frame, by resampling 
#' parametrically from a logNormal distribution with mean \code{hmm.pars$Et} and variance-covariance
#' matrix \code{hmm.pars$Sigma.Et}.
#' 
#' Function \code{\link{fit.hmltm}} is called to estimate detection probabilities and related things
#' for every bootstrap resample.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item{callist:}{ input reflection: everything passed to the function, bundled into a list}
#' \item{bs:}{ a list containing (a) a Bxn matrix \code{$phats} in which each row is the estimated 
#' detection probabilities for each of the n bootstrapped detections, (b) a Bxn matrix \code{$pars} 
#' in which each row is the estimated detection hazard parameters, (c) the following vectors 
#' of length B with estimates from each bootstrap: \code{$p0} (mean estimated p(0) over all 
#' detections), \code{$phat} (mean estimated detection probability over all detections), and (d) 
#' a Bx2 matrix \code{$b.Et} in which each row is the mean times unavailable and available.
#' }
#' }
#' 
bootstrap.p.with.Et=function(dat,pars,hfun,models,survey.pars,hmm.pars,
                           control.fit,control.opt,fixed.avail=FALSE,B=999){
  n=length(dat$x)
  npar=length(pars)
  b.p0=b.phat=rep(NA,B)
  b.pars=matrix(rep(NA,B*npar),ncol=npar)
  b.phats=matrix(rep(NA,B*n),ncol=n)# matrix for detection probs of each individual.
  # bootstrap availability parameters
  if(!fixed.avail) {
#    b.Et=mvrnorm(B,hmm.pars$Et,hmm.pars$Sigma.Et) # resample availability parameters
    lN=Npars.from.lNpars(hmm.pars$Et,hmm.pars$Sigma.Et) # get normal parameters corresponding to lognormal(mu,Sigma)
    b.Et=exp(mvrnorm(B,lN$mu,lN$Sigma)) # resample availability parameters on lognormal scale
  }
  # resample sightings data and re-estimate, using a resample of availability paramters:
  for(nb in 1:B){
    # resample detection locations
    samp.ind=bsample(1:n,size=n,replace=TRUE) # resample sightings data indices with replacement
    b.dat=dat[samp.ind,]# get resampled data
    # create new hmm.pars object with resampled availability parameters
    if(!fixed.avail) {
      pi21=1/b.Et[nb,2]
      pi12=1/b.Et[nb,1]
      Pi=matrix(c((1-pi12),pi12,pi21,(1-pi21)),nrow=2,byrow=TRUE)
      delta=compdelta(Pi)
      pm=c(0.0,1.0)
      b.hmm.pars=list(pm=pm,Pi=Pi,delta=delta,b.Et[nb],hmm.pars$Sigma.Et)
    }else {
      b.hmm.pars=hmm.pars
    }
    # refit model
    b.fit=fit.hmltm(b.dat,pars=pars,FUN=hfun,models=models,survey.pars=survey.pars,hmm.pars=b.hmm.pars,control.fit=control.fit,control.optim=control.opt)
    #    b.p0[nb]=b.fit$p[1]
    b.p0[nb]=b.fit$pzero
    b.phat[nb]=b.fit$phat
    b.phats[nb,]=b.fit$phats
    b.pars[nb,]=b.fit$fit$par
    cat(paste("done",nb,"\n"))
    flush.console()
  }
  # package results and return
  callist=list(dat=dat,pars=pars,hmm.pars=hmm.pars,hfun=hfun,models=models,
               survey.pars=survey.pars,control.fit=control.fit,control.optim=control.opt)
  bs=list(phats=b.phats,pars=b.pars,p0=b.p0,phat=b.phat,b.Et=b.Et)
  return(list(callist=callist,bs=bs))
}



#' @title Detection probability bootstrap with availability HMM.
#'
#' @description
#' Nonparametric bootstrap of detection data with re-estimation of detection probabilities. If 
#' \code{fixed.avail}=FALSE, does nonparametric resampling of availability HMM parameters contained
#' in hmm.pars.bs for every resample of detection data.
#' 
#' @param dat detection data frame constructed by removing all rows with no detections from a 
#' data frame of the sort passed to \code{\link{est.hmltm}}.
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars.bs multiple sets of availability hmm parameters, as output by 
#' \code{\link{hmmpars.boot}}. 
#' @param control.fit list controlling fit, as for \code{\link{est.hmltm}}.
#' @param control.opt list controlling function \code{\link{optim}}, as for \code{\link{est.hmltm}}.
#' @param fixed.avail if TRUE, hmm.pars is treated as fixed, else element \code{$Et} is parametrically 
#' resampled.
#' @param B number of bootstrap replicates.
#' @param silent argument of function \code{\link{try}}, controlling error message reporting.
#' 
#' @details
#' The rows of data frame \code{dat} are resampled with replacement to create new data frames with as 
#' many detections as were in \code{dat}. If \code{fixed.avail}=TRUE, then a new set of availability 
#' HMM parameters is obtainded by sampling iwth replacement from \code{hmm.pars.bs}.
#' 
#' Function \code{\link{est.hmltm}} is called to estimate detection probabilities and related things
#' for every bootstrap resample.
#' 
#' @return 
#' A list with the following elements:
#' \itemize{
#' \item{callist:}{ input reflection: everything passed to the function, bundled into a list}
#' \item{bs:}{ a list containing (a) a Bxn matrix \code{$phats} in which each row is the estimated 
#' detection probabilities for each of the n bootstrapped detections, (b) a Bxn matrix \code{$pars} 
#' in which each row is the estimated detection hazard parameters, and (c) the following vectors 
#' of length B with estimates from each bootstrap: \code{$Et} (mean times available and unavailable), 
#' \code{$p0} (mean estimated p(0) over all detections), \code{$phat} (mean estimated detection
#' probability over all detections), \code{$convergence} convergence diagnostic from \code{optim}.}
#' }
#' 
bootstrap.p.with.hmm=function(dat,pars,hfun,models,survey.pars,hmm.pars.bs,
                            control.fit,control.opt,fixed.avail=FALSE,B=999,silent=FALSE){
  n=length(dat$x)
  npar=length(pars)
  conv=b.p0=b.phat=rep(NA,B)
  b.pars=matrix(rep(NA,B*npar),ncol=npar)
  b.phats=matrix(rep(NA,B*n),ncol=n) # matrix for detection probs of each individual.
  # bootstrap availability parameters
  if(!fixed.avail) { # resample availability parameters
    nhmm=dim(hmm.pars.bs)[1]
    if(nhmm==1) stop("Only one set of hmm pars. need multiple sets of pars if not fixed avail.")
    reps=bsample(1:nhmm,B,replace=TRUE)
  }
  for(nb in 1:B) {
    if(!fixed.avail)  {
      b.hmm.pars=unvectorize.hmmpars(hmm.pars.bs[reps[nb],])
    }else {
      b.hmm.pars=hmm.pars.bs
    }
    if(is.element("pm",names(b.hmm.pars))) names(b.hmm.pars)[which(names(b.hmm.pars)=="pm")]="pm" # to fix naming cock-up when creating hmmpars.bs
    # resample detection locations
    samp.ind=bsample(1:n,size=n,replace=TRUE) # resample sightings data indices with replacement
    b.dat=dat[samp.ind,,drop=FALSE]# get resampled data
    names(b.dat)=names(dat)
    # refit model
    b.fit=try(fit.hmltm(b.dat,pars=pars,FUN=hfun,models=models,survey.pars=survey.pars,
                        hmm.pars=b.hmm.pars,control.fit=control.fit,control.optim=control.opt),
              silent=silent)
    if((class(b.fit)=="try-error")) {
      conv[nb]=-999
      b.p0[nb]=-999
      b.phat[nb]=-999
      b.phats[nb,]=rep(-999,n)
      b.pars[nb,]=rep(-999,length(pars))
    }else {
      conv[nb]=b.fit$fit$convergence
      b.p0[nb]=b.fit$p[1]
      b.phat[nb]=b.fit$phat
      b.phats[nb,]=b.fit$phats
      b.pars[nb,]=b.fit$fit$par
      cat(paste("done",nb,"\n"))
    }
    flush.console()
  }
  # package results and return
  callist=list(dat=dat,pars=pars,hmm.pars=hmm.pars.bs,hfun=hfun,survey.pars=survey.pars,
               control.fit=control.fit,control.optim=control.opt)
#  bs=list(hmm.pars=hmm.pars,p0=b.p0,phat=b.phat,phats=b.phats,pars=b.pars,convergence=conv)
  bs=list(phats=b.phats,pars=b.pars,p0=b.p0,phat=b.phat,convergence=conv)
  return(list(callist=callist,bs=bs))
}



#' @title Summarise detection probability bootstrap results.
#'
#' @description
#' Uses bootstrap results from \code{\link{bootstrap.p.with.Et}} or 
#' \code{\link{bootstrap.p.with.hmm}} to work out bootstrap means, variance estimates, CVs and
#' confidence intervals.
#' 
#' @param bs output from \code{\link{bootstrap.p.with.Et}} or \code{\link{bootstrap.p.with.hmm}}.
#' @param probs lower and upper percentile points for confidence interval reporting.
#' @param pcut minimum estimated detection probability to use. This is a quick and dirty method to 
#' robustify against small detection probability estimates skewing the distribution of \code{1/phat} 
#' badly for small samples. It is ad-hoc. If you use it, do histogram of $bs$phat to see if there is 
#' a reasonable cutpoint.
#' 
#' @return 
#' Returns a list with elements
#' \itemize{
#' \item{nboot:}{ number of bootstrap estimates used in constructing bootstrap statistics.}
#' \item{nbad:}{ number of bad estimates excluded from results.}
#' \item{parcov:}{ parameter estimate variance-covariance matrix.}
#' \item{parcorr:}{ parameter estimate correlation matrix.}
#' \item{bests:}{ bootstrap estimate statistics, comprising meand, standard error, percentage CV and 
#' confidence interval limits for: (a) estimated mean detection probability, (a) 1/(estimated mean 
#' detection probability), (c) estimated p(0), ad (d) detection hazard function parameters.}
#' }
#' 
bootsum.p=function(bs,probs=c(0.025,0.975),pcut=0){
  #--------------------------------------------------------------------------------
  # bs is output from bootstrap.p.with.Et() or bootstrap.with.hmm()
  # pcut is a quick and dirty min phat to allow - robustifies 1/phat for small 
  # samples, although it is ad-hoc. Do hist of $bs$phat to see if there is a 
  # reasonable cutpoint.
  #--------------------------------------------------------------------------------
  nboot=length(bs$bs$p0)
  if(is.null(bs$bs$convergence)) keep=which(bs$bs$p0>=0 & bs$bs$phat>pcut)
  else keep=which(bs$bs$p0>=0 & bs$bs$convergence==0 & bs$bs$phat>pcut)
  nbad=nboot-length(keep)
  npar=dim(bs$bs$par)[2]
  cinames=paste(as.character(probs*100),"%",sep="")  
  bests=matrix(rep(NA,(3+npar)*5),ncol=5)
  colnames(bests)=c("mean","std.err.","%CV",cinames)
  parnames=paste("par",as.character(1:npar),sep="")
  rownames(bests)=c("1/phat","phat","p(0)",parnames)
  # relative density
  bests[1,1]=mean(1/bs$bs$phat[keep])
  bests[1,2]=sd(1/bs$bs$phat[keep])
  bests[1,3]=sd(1/bs$bs$phat[keep])/mean(1/bs$bs$phat[keep])*100
  bests[1,4:5]=quantile(1/bs$bs$phat[keep],probs=probs)
  # detection probability
  bests[2,1]=mean(bs$bs$phat[keep])
  bests[2,2]=sd(bs$bs$phat[keep])
  bests[2,3]=sd(bs$bs$phat[keep])/mean(bs$bs$phat[keep])*100
  bests[2,4:5]=quantile(bs$bs$phat[keep],probs=probs)
  # p(0)
  bests[3,1]=mean(bs$bs$p0[keep])
  bests[3,2]=sd(bs$bs$p0[keep])
  bests[3,3]=sd(bs$bs$p0[keep])/mean(bs$bs$p0[keep])*100
  bests[3,4:5]=quantile(bs$bs$p0[keep],probs=probs)
  # parameters
  for(i in 1:npar) {
    bests[3+i,1]=sd(bs$bs$par[keep,i])
    bests[3+i,2]=mean(bs$bs$par[keep,i])
    bests[3+i,3]=sd(bs$bs$par[keep,i])/mean(bs$bs$par[keep,i])*100
    bests[3+i,4:5]=quantile(bs$bs$par[keep,i],probs=probs)
  }
  parcov=cov(bs$bs$par)
  parcorr=cov2cor(parcov)
  rownames(parcov)=rownames(parcorr)=colnames(parcov)=colnames(parcorr)=parnames
  return(list(nboot=nboot,nbad=nbad,bests=bests,parcov=parcov,parcorr=parcorr))
}

#-------------------------- End Bootstrap functions ------------------------
#                          ------------------------------


#                       --------------------------
#---------------------- Start Simulation functions ----------------------

#' @title Calculate CDF of forward distance.
#'
#' @description
#' Calculates CDF(y, given x)-u for single forward distance y, perpendicular distance x, and 
#' \code{u~runif(1,0,1)}. Currently only for models with no covariates.
#' 
#' @param x perpendicular distance.
#' @param y forward distance.
#' @param hfun detection hazard function name.
#' @param b likelihood function parameters.
#' @param pm state-dependent Bernoulli parameters.
#' @param Pi Markov model transition probability matrix.
#' @param delta Markov model stationary distribution.
#' @param ymax forward distance by which detection hazard has fallen to zero.
#' @param dy Markov model distance step size.
#' @param theta.f REDUNDANT must = 0.
#' @param theta.b REDUNDANT must = 90.
#' @param u uniform random variable (scalar).
#' 
cdfy.u=function(y,x,hfun,b,pm,Pi,delta,ymax,dy,theta.f,theta.b,u)
  #----------------------------------------------------------------------------
# For use with function uniroot.
# Calculates CDF(y)-u for single x,y and u~runif(1,0,1)
# Calls p.xy and calculates difference from u.
#----------------------------------------------------------------------------
{
  cdf=p.xy(x,y,hfun,b,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=TRUE)
  return(cdf-u)
}


#' @title Tabulate point and interval estimates.
#'
#' @description
#' Combines point estimates from \code{est.hmltm} and bootstrap estimates from \code{bs.hmltm} 
#' to produce a table ready for insertion into a paper. Does not include any strata with no
#' detections.
#' 
#' @param est output from \code{est.hmltm}.
#' @param bs output  from \code{bootsum} using the same model that created \code{est}.
#' 
estable=function(est,bs){
  est=est$point$ests
  nonzeros=which(est$n>0)
  if(!keepzeros) est=est[nonzeros,]
  # replace stratum.Area with p
  trow=length(est[,1])
  Nc=est$Ngroups[trow]*sum(est$covered.area[-trow])/sum(est$stratum.Area[-trow])
  p=est$n/(est$covered.area*est$Dgroups)
  p[trow]=est$n[trow]/Nc
  nL=signif(est$n/est$L,3)
  cv.nL=round(100*bs$cv[,"n/L"],1)
  cv.p=round(100*bs$cv[,"p"],1)
  cv.Ngrp=round(100*bs$cv[,"Ngroups"],1)
  cv.E.s=round(100*bs$cv[,"mean.size"],1)
  cv.N=round(100*bs$cv[,"N"],1)
  lcl.nL=bs$lower[,"n/L"]
  lcl.p=bs$lower[,"p"]
  lcl.Ngrp=bs$lower[,"Ngroups"]
  lcl.E.s=bs$lower[,"mean.size"]
  lcl.N=bs$lower[,"N"]
  ucl.nL=bs$upper[,"n/L"]
  ucl.p=bs$upper[,"p"]
  ucl.Ngrp=bs$upper[,"Ngroups"]
  ucl.E.s=bs$upper[,"mean.size"]
  ucl.N=bs$upper[,"N"]
  ci.N=paste("(",round(lcl.N),"; ",round(ucl.N),")",sep="")
  ci.E.s=paste("(",signif(lcl.E.s,3),"; ",signif(ucl.E.s,3),")",sep="")
  ci.Ngrp=paste("(",round(lcl.Ngrp),"; ",round(ucl.Ngrp),")",sep="")
  out=data.frame(Strat=est$stratum,A=round(est$stratum.Area),n=est$n,L=round(est$L,1),
                 a=round(est$covered.area,1),nL=signif(nL,3),cv.nl=cv.nL,p=signif(p,3),cv.p=cv.p,
                 Ngrp=round(est$Ngroups),cv.Ngrp=cv.Ngrp,ci.Ngrp=ci.Ngrp,
                 E.s=signif(est$mean.size,3),cv.E.s=cv.E.s,ci.E.s=ci.E.s,
                 N=round(est$N),cv.N=cv.N,ci.N=ci.N)
  return(out)
}




#' @title Simulate line transect survey data using availability HMM.
#'
#' @description
#' Simulates line transect survey data in which animal availability is simulated by generating
#' time series of availability data from the hidden Markov model(s) with parameters specified in
#' argument \code{hmm.pars}. Currently only for models with no covariates.
#' 
#' @param nw number of animals to simulate.
#' @param hfun detection hazard function name.
#' @param pars detection hazard function parameters.
#' @param hmm.pars availability hidden Markov model parameters (as per \code{\link{est.hmltm}}). 
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param print.progress if TRUE prints progress through simulations.
#'
#' @details
#' Simulates (x,y) observations from a population of nw animals using detection
#' hazard specified by hfun & pars, and HMM availability process specified by hmm.pars.
#' Survey parameters are in survey.pars.
#'
#' @return
#' Data frame with two elements: \code{$x} and \code{$y}, being perpendicular and forward distances
#' of detected animals, respectively.
#' 
simhmltm=function(nw,hfun,pars,hmm.pars,survey.pars,print.progress=TRUE){
  pm=hmm.pars$pm
  Pi=hmm.pars$Pi
  delta=hmm.pars$delta
  theta.f=survey.pars$theta.f
  theta.b=survey.pars$theta.b
  W=survey.pars$W
  ymax=survey.pars$ymax
  dy=survey.pars$dy
  b=tfm(pars,hfun)
  
  tmax=ceiling(ymax/dy) # max time units
  ymax=tmax*dy # ymax adjusted to be divisible by dy
  starty=rep(ymax,nw)-runif(nw,0,dy) # randomise start location within inteval dy of ymax
  xx=runif(nw,0,W) # generate perp dists of animals
  px=p.xy(xx,NULL,hfun,rep(b,nw),pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=TRUE) # calculate p(see|x)
  u=runif(nw,0,1)# random cdf(y) values for nw animals
  seen=which(px>=u) # identify those detected by the time they leave visible area
  nseen=length(seen)
  if(print.progress) cat(" Seen ",nseen,"\n","---------\n")
  if(nseen>0){ 
    pseen=px[seen] # record p(see) for detected animals
    x=xx[seen] # x-values of detected animals
    #    mincdf=p.xy(x,rep(ymax,nseen),hfun,b,pm,Pi,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=TRUE) # cdf value at ymax for all x values
    y=rep(NA,nseen) # initialise y
    ui=u[seen] # cdf(y) values for detected animals
    rootol=max(pseen)/1000 # tolerance for uniroot()
    for(i in 1:nseen){
      yi=starty[i]-(0:(tmax+1))*dy # set up disctete y-locations for time animal in view
      yi=yi[yi>=0] # exclude those that are behind abeam
      if(yi[length(yi)]>0) yi=c(yi,0) # need 0 in for consistency with pseen above (which calculates all the way to y=0)
      #      p=p.xy(rep(xi[i],length(yi)),yi,hfun,b=tfm(pars,hfun),pm,Pi,ymax,dy,theta.f,theta.b,cdf=TRUE) # calc prob(seen by y) for all y>y[i]
      ##      if(max(p)<pseen[i]) p=c(pseen[i],p) # add p(see by y=0) to p if it is not there (because of discretizatoin of yi)
      #      seenaty=max(which(p<=ui[i])) # index of closest y at which p(y)<=cdf(y) (vector y starts at ymax and gets smaller)
      mincdf=p.xy(x=x[i],y=range(yi)[2],hfun,b,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=TRUE)
      if(ui[i]<mincdf) {
        y[i]=ymax
        warning(paste("Observation y[",i,"]= ",y[i],">ymax put at ymax= ",ymax,".\n",sep=""))
      } else {
        cdf1=p.xy(x=x[i],y=range(yi)[1],hfun,b,pm,Pi,delta,ymax,dy,theta.f,theta.b,ally=FALSE,cdf=TRUE)
        uu=ui[i]
        if((cdf1<uu & mincdf<uu) | (cdf1>uu & mincdf>uu)) {
          cat(" (x,y) = (",x[i],",",range(yi),"); \n minCDF=",mincdf,"\n u=",uu,"\n maxCDF=",cdf1,"\n")
          #          cat("i, mincdf[i]=",i,", ",mincdf[i],"\n")
          #          cat("mincdf[i]=",mincdf,"\n")
        }
        cdf=uniroot(cdfy.u,interval=range(yi),tol=rootol,x=x[i],hfun=hfun,b=b,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,theta.f=theta.f,theta.b=theta.b,u=ui[i])
        y[i]=cdf$root
        if(is.null(y[i])) {
          cat("cdf(y)= ",ui[i],"\n")
          cat("y= ",y[i],"\n")
          cat("p= ",pseen[i],"\n")
          stop("Error in generating y.")
        }
      }
      if(print.progress) cat("Done ",i,"\n")
    }
  }
  return(data.frame(x=x,y=y))
}




#' @title Simulate line transect survey data using availability time series.
#'
#' @description
#' Simulates line transect survey data in which animal availability is simulated by sampling from
#' time series of availability data (from tagged animals, for example). Currently only for models 
#' with no covariates.
#' 
#' @param adat list of m>=1 vectors containing availability binary time series or depths time series.
#' @param xmax maximum perpendicular distance to simulate.
#' @param ymax maximum forward distance to simulate (must be at or beyond point that detection hazard
#' function is effectively zero).
#' @param spd speed observer is moving.
#' @param animals a vector of up to m integers specifying which members of the list \code{adat} are to 
#' be used to simulate availability. (e.g. \code{animals=c(1,3,5)} says the first, third and fifth time
#' series in \code{adat} are to be used).
#' @param hfun detection hazard function name. 
#' @param pars detection hazard function parameters. 
#' @param N number of animals within distance \code{xmax} of the transect line and hence subject to 
#' being detected (or possibly not).
#' @param dmax if not NULL, then it is the depth above which animals are considered to be available, and
#' \code{adat} is taken to be time series of depths, i.e. all \code{adat[[i]]<=dmax} are taken as times
#' animals are available. If \code{dmax} is NULL, then adat must be binary and all 1's are taken as times
#' animals are available.
#' @param seed random number seed.
#' @param poiss if TRUE, the availability time series are randomly permuted, hence generating the 
#' equivalent of a Poisson availability process.
#'
#' @details
#' Simulates (x,y) observations from a population of N animals using detection
#' hazard specified by hfun & pars, and resampling from the availability time series in adat for
#' animal availability.
#' Survey parameters are in survey.pars.
#'
#' @return
#' Data frame with two elements: \code{$x} and \code{$y}, being perpendicular and forward distances
#' of detected animals, respectively.
#' 
simhmltm.w=function(adat,xmax,ymax,spd,animals,hfun,pars,N,dmax=NULL,seed=NULL,poiss=FALSE)
{
  h=match.fun(hfun)
  b=tfm(pars,hfun) # this is inefficient ('cause just back-transform in h()) but makes for uniform handling outside this function
  if(!is.null(seed)) set.seed(seed)
  ytmax=ceiling(ymax/spd) # needs to be integer cause is number of time units
  wst=wa=adata=adat[animals] # select which time series to use
  if(is.null(dmax)) {# in this case adat is binary availability data
    for(i in 1:length(animals)) {
      # Determine the times when the animals are available and store in a list
      wa[[i]]=adata[[i]]>=0.5 # TRUE for 1s, FALSE for 0s
      # Determine the last possible starting times of sequences of length (2*ymax+1) for animal 
      wst[i]=length(adata[[i]])-ytmax
    }
  } else {
    for(i in 1:length(animals)) {
      # Determine the times when the animals are available and store in a list
      wa[[i]]=adata[[i]]<=dmax # TRUE for shallower than dmax, FALSE for deeper
      # Determine the last possible starting times of sequences of length (2*ymax+1) for animal 
      wst[i]=length(adata[[i]])-ytmax
    }
  }
#  # Determine the times when the animals are available and store in a list
#  wa<-list(adat$w1<=dmax, adat$w2<=dmax, adat$w3<=dmax, adat$w4<=dmax, adat$w5<=dmax, adat$w6<=dmax, adat$w7<=dmax,adat$w8<=dmax)    
#  # Determine the last possible starting times of sequences of length (2*ymax+1) for animal 
#  wst<-c(length(adat$w1)-ytmax, length(adat$w2)-ytmax, length(adat$w3)-ytmax, length(adat$w4)-ytmax, length(adat$w5)-ytmax,length(adat$w6)-ytmax,length(adat$w7)-ytmax,length(adat$w8)-ytmax)
  # Determine the number of times available and the number of unavailable
  if(length(animals)>1) iw=bsample(animals,N,replace=TRUE) else iw=rep(animals,N)
  xobs<-yobs<-rep(NULL,N) # initialise variables holding locations of detections
  y0=(0:ytmax)*spd # unadjusted y locations
  y1=runif(N,0,spd) # randomise closest location between 0 and 1/spd for all animals
  nseen=0
  for(i in 1:N) { # run through the animals
    y=y0+y1[i] # y locations for this animal (with randomised closest y)
    y=y[y<=ymax]
    surf<-wa[[iw[i]]]
    # If poisson==TRUE shuffle the surfacings (for comparison)
    if(poiss){surf<-bsample(surf)}
    # Select the starting indices for the run
    runstart <-bsample(1:wst[iw[i]],1)
    # Generate the x-value
    x<-runif(1,0,xmax)
    # Fly over the runs record positions of the detections.
    up=which(surf[runstart:(runstart+ytmax)])
    nup=length(up)
    if(nup>0) {
      p=h(rep(x,nup),y[up],b)
      u=runif(nup)
      saw=which(p>u)
      if(length(saw)>0) {
        nseen=nseen+1
        xobs[nseen]=x
        yobs[nseen]=max(y[saw])
      }
    }
  }
  if(nseen>0) dat=data.frame(x=xobs[1:nseen],y=yobs[1:nseen])
  else dat=NULL
  return(dat)
}



#' @title Simulation test hmltm performance
#'
#' @description
#' Simulates line transect survey with stochastic animal availability, estimates detection probability
#' and related parameters, and reports bias, variance of estimators. Currently only for models 
#' with no covariates.
#' 
#' @param simethod either "hmm" (for availability simulation with HMM) or "animals" (for availability 
#' simulation by resampling availability or depth dime series).
#' @param shfun detection hazard function name for simulation.
#' @param spars detection hazard function parameters for simulation.
#' @param ehfun detection hazard function name for estimation.
#' @param parstart detection hazard function parameter start vallues for estimation.
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param shmm.pars hmm.pars list (as for \code{\link{est.hmltm}}) to be used in simulating.
#' @param ehmm.pars hmm.pars list (as for \code{\link{est.hmltm}}) to be used in estimating.
#' @param animal which animals' availability data to use (only need if simethod="animals").
#' @param adat availability data (only need if simethod="animals"): list of m>=1 vectors 
#' containing availability binary time series or depths time series.
#' @param control.opt optimization parameters passed to \code{\link{optim}}, as for \code{\link{est.hmltm}}.
#' @param control.fit fit parameters, as for \code{\link{est.hmltm}}.
#' @param nsim number of simulations to do.
#' @param report.progress if TRUE, reports progress as does simulations.
#' @param En Mean sample size for simulations.
#' @param doplots if TRUE does plots of simulaton results.
#' @param varest if TRUE, estimates variance (using Hessian matrix).
#' @param hmmpars.bs output from \code{\link{hmmpars.boot}}; required if simethod="hmm" and the 
#' availability HMM parameters are not to be treated as fixed, in which case the
#' availability HMM parameters for each simulation are obtained by sampling with replacement from 
#' \code{hmmpars.bs}. If hmmpars.bs is NULL, the availability HMM parameters are treated as fixed
#' and equal to ehmm.pars.
#' @param print.n if TRUE, prints sample size for each simulation.
#' @param silent parameter of \code{\link{try}}, controlling error reporting.
#' @param nx number of perpendicular distances intervals to use in evaluating detection 
#' probability, p(x).
#'
#' @return
#' list with these elements:
#' \itemize{
#' \item{N:}{ population size.}
#' \item{esw:}{ true effective strip half-width.}
#' \item{p0:}{ true p(0).}
#' \item{p.:}{ true mean detection probability.}
#' \item{meann:}{ mean sample size across simulations.}
#' \item{biasNhat:}{ bias of estimated N.}
#' \item{biaseswhat:}{ bias of estiamted esw.}
#' \item{biasp0hat:}{ bias of estimated p(0).}
#' \item{biasphat:}{ bias of estimated mean detection probability.}
#' \item{biasinvphat:}{ bias of estimated inverse mean detection probability.}
#' \item{n:}{ sample sizes.}
#' \item{parest:}{ parameter estimates from each simulation.}
#' \item{p0hat:}{ p(0) estimates from each simulation.}
#' \item{phat:}{ mean detection probability estimates from each simulation.}
#' \item{invphat:}{ inverse mean detection probability estimates from each simulation.}
#' \item{Nhat:}{ abunance estimates from each simulation.}
#' \item{invpse:}{ estimated standard error of inverse mean detection probability estimates from each 
#' simulation (if varest=TRUE).}
#' \item{p0se:}{ estimated standard error of p(0) estimates from each simulation (if varest=TRUE).}
#' }
#' 
simest=function(simethod="hmm",shfun,spars,ehfun=shfun,parstart=spars,survey.pars,
                shmm.pars,ehmm.pars=shmm.pars,animal=NULL,adat=NULL,
                control.opt,control.fit,nsim=50,report.progress=TRUE,En=100,
                doplots=FALSE,varest=FALSE,hmmpars.bs=NULL,print.n=FALSE,
                silent=FALSE,nx=100)
{
  if(simethod!="animals" & simethod!="hmm" & simethod!="Et") stop("simethod must be 'animals' or 'hmm' or 'Et'")
  if(simethod=="animals" & (is.null(animal) | is.null(adat))) stop("Need to pass objects animal & adat with availability data") 
  if(varest) control.fit$hessian=TRUE
  
  Nhat=n=phat=p0hat=eswhat=invpse=p0se=rep(NA,nsim)
  Nhat.x=phat.x=p0hat.x=eswhat.x=invpse.x=p0se.x=rep(NA,nsim)
  parest=parest.x=matrix(rep(NA,nsim*length(parstart)),nrow=nsim)
  W=survey.pars$W
  #    pp=plot.pxfy0(shfun,spars,survey.pars,shmm.pars,doplots=FALSE)
  #    p.=estp(pars=spars,hfun=shfun,survey.pars=survey.pars,hmm.pars=shmm.pars)
  #    esw=p.*W
  xs=seq(0,W,length=nx)
  # non-Chris change:
  p=hmltm.px(x=xs,pars=spars,hfun=shfun,models=NULL,cov=NULL,survey.pars=survey.pars,
             hmm.pars=shmm.pars)
#  p=hmltm.px(x=xs,pars=spars,hfun=shfun,survey.pars=survey.pars,hmm.pars=shmm.pars)
  p0=p[1]
  esw=sintegral(p,xs)
  p.=esw/W
  N=round(En/p.)
  #    p0=pp$p0
  simestw=0
  #  simestw=simestw.x=0
  if(!is.null(hmmpars.bs) & simethod=="hmm") nhmm=dim(hmmpars.bs)[1]    
  if(!is.null(hmmpars.bs) & simethod=="Et") b.Et=mvrnorm(nsim,hmmpars.bs$Et,hmmpars.bs$Sigma.Et) # resample availability parameters
  
  for(i in 1:nsim) {
    if(!is.null(hmmpars.bs) & simethod=="hmm") {
      rep=bsample(1:nhmm,1)
      ehmm.pars=unvectorize.hmmpars(hmmpars.bs[1,])
      if(is.element("pm",names(ehmm.pars))) names(ehmm.pars)[which(names(ehmm.pars)=="pm")]="pm" # to fix naming cock-up when creating hmmpars.bs
    }      
    if(!is.null(hmmpars.bs) & simethod=="Et") {
      pi21=1/b.Et[i,2]
      pi12=1/b.Et[i,1]
      Pi=matrix(c((1-pi12),pi12,pi21,(1-pi21)),nrow=2,byrow=TRUE)
      delta=compdelta(Pi)
      pm=c(0.0,1.0)
      ehmm.pars=list(pm=pm,Pi=Pi,delta=delta,hmmpars.bs$Et,hmmpars.bs$Sigma.Et)
    }else {
      ehmm.pars=hmmpars.bs
    }
    
    if(simethod=="animals") { # simulate by resampling observed availability pattern(s) in adat
      dat<-simhmltm.w(adat,xmax=W,ymax=survey.pars$ymax,spd=survey.pars$spd,animals=animal,hfun=shfun,spars,N,dmax=2,seed=NULL,poiss=FALSE)
      n[i]=length(dat$x); if(print.n) cat("n=",n[i],"\n")
      simestw=try(fit.hmltm(xy=dat,pars=parstart,FUN=shfun,models=NULL,survey.pars=survey.pars,ehmm.pars,control.fit=control.fit,control.optim=control.opt),silent=silent)
    } else if(simethod=="hmm" | simethod=="Et") { # simulate using HMM specified in hmm.pars
      dat<-simhmltm(N,shfun,spars,shmm.pars,survey.pars,print.progress=FALSE)
      n[i]=length(dat$x); if(print.n) cat("n=",n[i],"\n")
      simestw=try(fit.hmltm(xy=dat,pars=parstart,FUN=shfun,models=NULL,survey.pars=survey.pars,ehmm.pars,control.fit=control.fit,control.optim=control.opt),silent=silent)
    } 
    if((class(simestw)=="try-error")) {
      parest[i,]=-999
      p0hat[i]=-999
      phat[i]=-999
      eswhat[i]=-999
      Nhat[i]=-999
      if(varest) {
        invpse[i]=-999
        p0se[i]=-999
      }
    }else {
      parest[i,]=simestw$fit$par
      p0hat[i]=simestw$p[1]
      phat[i]=simestw$phat
      eswhat[i]=phat[i]*W
      Nhat[i]=n[i]/phat[i]
      if(varest) {
        if(is.null(simestw$fit$hessian)) {
          if(i==1) warning("No Hessian so can't get analytic variance estimate.")
          invpse[i]=NA
          p0se[i]=NA
        } else {
          invpse[i]=invHessvar("invp",simestw)$se
          p0se[i]=invHessvar("p0",simestw)$se
        }
      }
    }
    
    if(report.progress) print(paste("done",i))
  }
  
  mnN=mean(Nhat)
  biasNhat=(mnN-N)/N
  #    mnN.x=mean(Nhat.x)
  #    biasNhat.x=(mnN.x-N)/N
  
  eswhat=phat*W
  mnesw=mean(eswhat)
  biaseswhat=(mnesw-esw)/esw
  #    eswhat.x=phat.x*W
  #    mnesw.x=mean(eswhat.x)
  #    biaseswhat.x=(mnesw.x-esw)/esw
  
  mnphat=mean(phat)
  biasphat=(mnphat-p.)/p.
  #    mnphat.x=mean(phat.x)
  #    biasphat.x=(mnphat.x-p.)/p.
  
  invphat=1/phat
  mninvphat=mean(invphat)
  biasinvphat=(mninvphat-1/p.)/(1/p.)
  #    invphat.x=1/phat.x
  #    mninvphat.x=mean(invphat.x)
  #    biasinvphat.x=(mninvphat.x-1/p.)/(1/p.)
  
  mnp0=mean(p0hat)
  biasp0hat=(mnp0-p0)/p0   
  #    mnp0.x=mean(p0hat.x)
  #    biasp0hat.x=(mnp0.x-p0)/p0
  
  meann=mean(n)
  
  simlist=list(N=N,esw=esw,p0=p0,p.=p.,
               meann=meann,
               biasNhat=biasNhat,biaseswhat=biaseswhat,biasp0hat=biasp0hat,
               biasphat=biasphat,biasinvphat=biasinvphat,
               n=n,
               parest=parest,p0hat=p0hat,phat=phat,invphat=invphat,Nhat=Nhat,invpse=invpse,p0se=p0se)
  
  if(doplots){
    par(mfrow=c(3,2))
    
    hist(Nhat,xlab="Estimated N",xlim=range(N,Nhat))
    lines(rep(N,2),c(0,nsim),col="red",lwd=3)
    lines(rep(mnN,2),c(0,nsim),col="blue",lwd=3,lty=2)
    #      hist(Nhat.x,xlab="Estimated N",xlim=range(N,Nhat.x))
    #      lines(rep(N,2),c(0,nsim),col="red",lwd=3)
    #      lines(rep(mnN.x,2),c(0,nsim),col="blue",lwd=3,lty=2)
    
    hist(eswhat,xlab="Estimated ESW",xlim=range(esw,eswhat))
    lines(rep(esw,2),c(0,nsim),col="red",lwd=3)
    lines(rep(mnesw,2),c(0,nsim),col="blue",lwd=3,lty=2)
    #      hist(eswhat.x,xlab="Estimated ESW",xlim=range(esw,eswhat.x))
    #      lines(rep(esw,2),c(0,nsim),col="red",lwd=3)
    #      lines(rep(mnesw.x,2),c(0,nsim),col="blue",lwd=3,lty=2)
    
    hist(p0hat,xlab="Estimated p(0)",xlim=range(p0,p0hat))
    lines(rep(p0,2),c(0,nsim),col="red",lwd=3)
    lines(rep(mnp0,2),c(0,nsim),col="blue",lwd=3,lty=2)     
    #      hist(p0hat.x,xlab="Estimated p(0)",xlim=range(p0,p0hat.x))
    #      lines(rep(p0,2),c(0,nsim),col="red",lwd=3)
    #      lines(rep(mnp0.x,2),c(0,nsim),col="blue",lwd=3,lty=2)
  }
  
  return(simlist)
}


#' @title Summarise simulation results.
#'
#' @description
#' Summarises (and prints) output from \code{\link{simest}}.
#'  
#' @param simout output from \code{\link{simest}}.
#' @param brief if TRUE, prints briefer summary.
#' 
#sumsim=function(simout,p=NULL,p0=NULL,brief=FALSE){
sumsim=function(simout,brief=FALSE){
  p=simout$p.
  p0=simout$p0
  esw=simout$esw
  N=simout$N
  
  nsimtot=length(simout$n)
  keep=(simout$phat>=0) # these had estimation problems
  nsimbad=nsimtot-sum(keep)
  nsim=nsimtot-nsimbad
  
  mean.n=mean(simout$n[keep])
  sd.n=sd(simout$n[keep])
  cv.n=sd.n/mean.n
  
  bias.p0=mean(simout$p0hat[keep])/p0-1
  sd.p0=sd(simout$p0hat[keep])
  cv.p0=sd.p0/mean(simout$p0hat[keep])
  rmse.p0=sqrt(mean((simout$p0hat-p0)^2))
  ci.p0.bias=(mean(simout$p0hat[keep])+c(-1,1)*1.96*sd.p0/sqrt(nsim)-p0)/p0
  
  bias.p=mean(simout$phat[keep])/p-1
  sd.p=sd(simout$phat[keep])
  cv.p=sd.p/mean(simout$phat[keep])
  rmse.p=sqrt(mean((simout$phat-p)^2))
  ci.p.bias=(mean(simout$phat[keep])+c(-1,1)*1.96*sd.p/sqrt(nsim)-p)/p
  
  invp=1/p
  bias.invp=mean(1/simout$phat[keep])/invp-1
  sd.invp=sd(1/simout$phat[keep])
  cv.invp=sd.invp/mean(1/simout$phat[keep])
  rmse.invp=sqrt(mean((1/simout$phat-1/p)^2))
  ci.invp.bias=(mean(1/simout$phat[keep])+c(-1,1)*1.96*sd.invp/sqrt(nsim)-invp)/invp
  
  lower.invp=1/simout$phat[keep]-1.96*simout$invpse[keep]
  upper.invp=1/simout$phat[keep]+1.96*simout$invpse[keep]
  nci=length(lower.invp)
  cover.invp=sum(lower.invp<=invp & invp<=upper.invp)/nci*100
  
  lower.p0=simout$p0hat[keep]-1.96*simout$p0se[keep]
  upper.p0=simout$p0hat[keep]+1.96*simout$p0se[keep]
  nci=length(lower.p0)
  cover.p0=sum(lower.p0<=p0 & p0<=upper.p0)/nci*100
  
  if(brief) {
    cat("Simulation Summary\n")
    cat("--------------------------------------------------\n")
    cat("Number converged simulations:",nsim,"\n")
    cat("% bad simulations      :",100*nsimbad/nsimtot,"\n")
    cat("1/p Summary:\n")
    cat("-------------\n")
    cat("%Bias                           :",100*bias.invp,"\n")
    cat("Std Err 1/p                     :",sd.invp,"\n")
    cat("RMSE 1/p                        :",rmse.invp,"\n")
    cat("--------------------------------------------------\n")
  }else {
    cat("Simulation Summary\n")
    cat("--------------------------------------------------\n")
    cat("Number converged simulations:",nsim,"\n")
    cat("Number bad simulations      :",nsimbad,"\n")
    cat("Number simulations in total :",nsimtot,"\n")
    cat("(Results below exclude bad simulations)\n\n")
    cat("Mean sample size            :",mean.n,"\n")
    cat("Std dev of sample size      :",sd.n,"\n")
    cat("%CV of sample size          :",cv.n,"\n\n")
    cat("p(0) Summary:\n")
    cat("-------------\n")
    cat("Truth                           :",p0,"\n")
    cat("Simulation mean                 :",mean(simout$p0hat),"\n")
    cat("%Bias                           :",100*bias.p0,"\n")
    cat("95% Bias CI                     :",100*ci.p0.bias,"\n")
    cat("Std Err p0                      :",sd.p0,"\n")
    cat("% CV p0                         :",cv.p0*100,"\n")
    cat("RMSE p0                         :",rmse.p0,"\n")
    cat("Mean Inf matrix Std Err estimate:",mean(simout$p0se),"\n")
    cat("Inf matrix coverage prob        :",cover.p0,"\n\n")
    cat("1/p Summary:\n")
    cat("-------------\n")
    cat("Truth                           :",invp,"\n")
    cat("Simulation mean                 :",mean(1/simout$phat),"\n")
    cat("%Bias                           :",100*bias.invp,"\n")
    cat("95% Bias CI                     :",100*ci.invp.bias,"\n")
    cat("Std Err 1/p                     :",sd.invp,"\n")
    cat("%CV 1/p                         :",cv.invp*100,"\n")
    cat("RMSE 1/p                        :",rmse.invp,"\n")
    cat("Mean Inf matrix Std Err estimate:",mean(simout$invpse),"\n")
    cat("Inf matrix coverage prob        :",cover.invp,"\n")
    cat("--------------------------------------------------\n")
  }
}

#---------------------- End Simulation functions ----------------------
#                       ------------------------

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



#natparvar=function(hmmlt){
#  #-------------------------------------------------------------------------------------
#  # NOTE: Does not work with covariates - 'cause with the way covariates are implemened,
#  #      have more natural pars than fitting pars.
#  #-------------------------------------------------------------------------------------
#  if(!is.nullmodel(hmmlt$models)) stop("This function only implemented for models without covariates.")
#  # extract things we need:
#  hfun=hmmlt$h.fun
#  par=hmmlt$fit$par
#  b=tfm(hmmlt$fit$par,hfun) 
#  d.dpar=dinvt.dpar(b,hfun) # derivative of natural with respect to fitting pars
#  # use inverse Hessian to calculate variance
#  invHess=solve(hmmlt$fit$hessian)
#  VAR=t(diag(d.dpar))%*%invHess%*%diag(d.dpar)
#  SE=sqrt(diag(VAR))
#  CV=SE/par
#  
#  return(list(invHess=invHess,d.dpar=d.dpar,par=par,se=SE,cv=CV))
#}


#------------------------- End Analytic Variance Estimation functions -----------------------
#                          ------------------------------------------


########################################### End of Functions ################################################

########################################### Start of Data ################################################

#' @name aerial.survey
#' @title Simulated bowhead whale aerial survey dataset.
#' @docType data
#' @description Simulated data representing an aerial survey of bowhead whales from an aircraft 
#' flying at 46.3 m/s. These data can used together with the dataset \code{\link{bowhead.hmm.pars}} 
#' to estimate bowhead whale abundance.
#' @usage aerial.survey
#' @format An mrds data frame (with 2 rows per detection) with 86 observations on the following 11 
#' variables.
#'  \describe{
#'    \item{\code{stratum}:}{ stratum number (a numeric vector).}
#'    \item{\code{area}:}{ stratum surface area (a numeric vector).}
#'    \item{\code{transect}:}{ transect number (a numeric vector).}
#'    \item{\code{L}:}{ transect length (a numeric vector).}
#'    \item{\code{size}:}{ size of each detected group (NA if no detection).}
#'    \item{\code{object}:}{ unique identifier for each detection (a numeric vector; NA if no detection).}
#'    \item{\code{side}:}{ the side of the plane from which each detection was made (NA if no detection).}
#'    \item{\code{obs}:}{ the observer who made each detection (NA if no detection).}
#'    \item{\code{bf}:}{ Beaufort sea state at the time of each detection (NA if no detection).}
#'    \item{\code{x}:}{ perpendicular distances to detections (NA if no detection).}
#'    \item{\code{y}:}{ perpendicular distances to detections (NA if no detection).}
#'  }
#' @source Simulated
#' @examples
#'  data(aerial.survey)
NULL


#' @name beaked.ship
#' @title Beaked whale shipboard survey dataset from Alboran sea.
#' @docType data
#' @description Data from the 2008 & 2009 shipboard survey of bowhead whales in the Alboran sea. Sightings
#' are real, strata and transects are made up for illustration.
#' @usage beaked.ship
#' @format A data frame with 81 observations on the following 12 variables.
#'  \describe{
#'    \item{\code{stratum}:}{ stratum number (a numeric vector). In this dataset this is a dummy 
#'    value - just there to have something in this compulsory variable.}
#'    \item{\code{area}:}{ stratum surface area (a numeric vector). In this dataset this is a dummy 
#'    value - just there to have something in this compulsory variable.}
#'    \item{\code{transect}:}{ transect number (a numeric vector). In this dataset this is a dummy 
#'    value - just there to have something in this compulsory variable.}
#'    \item{\code{L}:}{ transect length (a numeric vector). In this dataset this is a dummy 
#'    value - just there to have something in this compulsory variable.}
#'    \item{\code{x}:}{ perpendicular distances to detections (NA if no detection).}
#'    \item{\code{y}:}{ perpendicular distances to detections (NA if no detection).}
#'    \item{\code{size}:}{ size of each detected group (NA if no detection).}
#'    \item{\code{bf}:}{ Beaufort sea state at the time of each detection (NA if no detection).}
#'    \item{\code{ht}:}{ Observation platform height in m (NA if no detection).}
#'    \item{\code{object}:}{ unique identifier for each detection (a numeric vector; NA if no detection).}
#'  }
#' @details Test dataset that contains only detections (transects without detections have been omitted). 
#' It is one of the datasets analysed in Borchers et al. (2013).
#' @source Ana Canadas.
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' @examples
#'  data(beaked.ship)
NULL


#' @name beaked.hmm.pars
#' @title Alboran sea beaked whale availability HMM parameters.
#' @docType data
#' @description Hidden Markov model (HMM) parameter estimates for beaked whale availability obtained 
#' from mean times available and unavailable observed by Ana Canadas (pers commn.) in thge Alboran sea. 
#' This dataset was used in the analyses of Borchers et al. (2013). 
#' @usage beaked.hmm.pars
#' @format A list with the following five elements.
#'  \describe{
#'    \item{\code{pm}:}{ a 2x1 matrix containing the vector of state-dependent Bernoulli availability 
#'    parameters. \code{pm[1,i]} is the probability of whale i being available given state 1 (the less 
#'    available behavoural state), and \code{pm[1,i]} is the probability of whale i being available given 
#'    state 1 (the more available behavoural state).}
#'    \item{\code{Pi}:}{ a 2x2x1 array containg the transition probability matrix. States can be 
#'    interpreted as behavioural states, one of which being a state in which the animal is more 
#'    likely to be available than when in the other state.}
#'    \item{\code{delta}:}{ a 2x1 matrix containing the stationary distribution of \code{Pi} for each 
#'    whale. So \code{delta[1,i]} is the probability that whale i is in behavioural state 1 when 
#'    observation starts, and \code{delta[2,i]} is the probability that it is in behavioural state 2 when 
#'    observation starts.}
#'    \item{\code{Et}:}{ a 2x1 matrix containing the expected times animals are available and unavailable (in seconds).}
#'    \item{\code{Sigma.Et}:}{ a 2x2 matrix containing variance-covariance matrix of the expected times animals are available and unavailable.}
#'  }
#'  
#' @source Canadas (pers commn.).
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' @examples
#'  data(beaked.hmm.pars)
NULL


#' @name bowhead.hmm.pars
#' @title Greenland bowhead whale availability HMM parameters.
#' @docType data
#' @description Hidden Markov model (HMM) parameter estimates for bowhead whale availability obtained 
#' from fitting HMMs (using library \code{\link{HiddenMarkov}}) to data from electronic depth-recording 
#' tags that were attached to eight bowhead whales from the West Greenland population. The tags 
#' generated 8 time series of durations between 2.6 and 53 hours, with depths recorded every second 
#' (see Laidre et al., 2007). Following previous practice (Heide-Jorgensen et al., 2007), animals were 
#' considered to be available for detection only when within 2 m of the surface. The time series were 
#' accordingly converted into binary availability time series and HMMs were fitted to these. This 
#' dataset was used in the analyses of Borchers et al. (2013) - albeit with survey data that had no
#' forward distances. 
#' @usage bowhead.hmm.pars
#' @format A list with the following three elements.
#'  \describe{
#'    \item{\code{Pi}:}{ a 2x2x8 array containg the 8 HMM transition probability matrices (one for each 
#'    tagged whale). States can be interpreted as behavioural states, one of which being a state in 
#'    which the animal is more likely to be available than when in the other state.}
#'    \item{\code{pm}:}{ a 2x8 matrix containing the 8 vectors of state-dependent Bernoulli availability 
#'    parameters. \code{pm[1,i]} is the probability of whale i being available given state 1 (the less 
#'    available behavoural state), and \code{pm[1,i]} is the probability of whale i being available given 
#'    state 1 (the more available behavoural state).}
#'    \item{\code{delta}:}{ a 2x8 matrix containing the stationary distribution of \code{Pi} for each 
#'    whale. So \code{delta[1,i]} is the probability that whale i is in behavioural state 1 when 
#'    observation starts, and \code{delta[2,i]} is the probability that it is in behavioural state 2 when 
#'    observation starts.}
#'  }
#'  
#'  @seealso The depth time series data are in object \code{\link{bowhead.depths}} and the binary 
#'  presence/absence data obtained from these are in object \code{\link{bowhead.adat}}.
#'  
#' @source Laidre et al. (2007).
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' Heide-Jørgensen, M. P., Laidre, K., Borchers, D. L., Samarrra,F., and Stern, H. 2007. Increasing 
#' abundance of bowhead whales in west greenland. Biology Letters 3, 577–580.
#' 
#' Laidre, K., Heide-Jørgensen, M. P., and Nielsen, T. 2007. Role of bowhead whale as a predator in 
#' West Rreenland. Marine Ecology Progress Series 346, 285–297.
#' 
#' @examples
#'  data(bowhead.hmm.pars)
NULL


#' @name bowhead.hmm.pars.bs
#' @title Bootstrapped Greenland bowhead whale availability HMM parameters.
#' @docType data
#' @description Bootstrapped Hidden Markov model (HMM) parameter estimates for bowhead whale 
#' availability, obtained using by passing \code{\link{bowhead.hmm.pars}} and 
#' \code{\link{bowhead.adat}} to \code{\link{hmmpars.boot}}.
#' @format A matrix in which each row is a set of HMM parameters converted to vector format by 
#' \code{\link{vectorize.hmmpars}}. Each row can be converted back to the format required by
#' \code{\link{est.hmltm}} using \code{\link{unvectorize.hmmpars}}
#' @usage bowhead.hmm.pars.bs
NULL


#' @name bowhead.adat
#' @title Greenland bowhead whale availability data time series.
#' @docType data
#' @description A list with elements \code{$a1} to \code{$a8}. Element \code{$ai} is a vector for whale
#' i containing 0s and 1s, with 0s indicating that the whale was deeper than 2m and 1s indicating it was not. These 
#' observations are 1 second apart.
#' @format A list of 8 numeric vectors. Their lengths are 191,048, 73,871, 24,673, 9,385, 37,689, 
#' 47,222, 15,490 and 39,989, respectively.
#' @usage bowhead.adat
NULL


#' @name bowhead.depths
#' @title Greenland bowhead whale depth data time series.
#' @docType data
#' @description A list with elements \code{$a1} to \code{$a8}. Element \code{$ai} has the depths of
#' tagged bowhead whale i each second.
#' @format A list of 8 numeric vectors. Their lengths are 191,048, 73,871, 24,673, 9,385, 37,689, 
#' 47,222, 15,490 and 39,989, respectively.
#' @usage bowhead.depths
NULL


#' @name porpoise.hmm.pars
#' @title Harbour porpoise availability HMM parameters.
#' @docType data
#' @description Hidden Markov model (HMM) parameter estimates for harbour porpoise availability obtained 
#' the dive durations and proportons of time available of 7 tagged harbour porpoise, as reported by
#' Westgate et al. (1995). HMM parameters were obtained using funtion \code{\link{make.hmm.pars.from.Et}}.
#' @usage porpoise.hmm.pars
#' @format A list with the following three elements.
#'  \describe{
#'    \item{\code{Pi}:}{ a 2x2x7 array containg the 7 HMM transition probability matrices (one for each 
#'    tagged porpoise). States can be interpreted as behavioural states, one of which being a state in 
#'    which the animal is more likely to be available than when in the other state.}
#'    \item{\code{pm}:}{ a 2x7 matrix containing the 8 vectors of state-dependent Bernoulli availability 
#'    parameters. pm[1,i] is the probability of porpoise i being available given state 1 (the less 
#'    available behavoural state), and pm[1,i] is the probability of porpoise i being available given 
#'    state 1 (the more available behavoural state).}
#'    \item{\code{delta}:}{ the stationary distribution of Pi for each porpoise. So delta[1,i] is the 
#'    probability that porpoise i is in behavioural state 1 when observation starts, and delta[2,i] is 
#'    the probability that it is in behavioural state 2 when observation starts.}
#'  }
#' @source Westgate et al. (1995).
#' @references 
#' Westgate, A. J., Read, A. J., Berggren, P., Koopman, H. N., and Gaskin, D. E. 1995. Diving behaviour 
#' of harbour porpoises, Phocoena phocoena. Canadian Journal of Fisheries and Aquatic Sciences 52, 
#' 1064-1073.
#' @examples
#'  data(porpoise.hmm.pars)
NULL
