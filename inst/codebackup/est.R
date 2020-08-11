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