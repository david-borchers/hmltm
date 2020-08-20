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


##' @title Tabulate point and interval estimates.
##'
##' @description
##' Combines point estimates from \code{est.hmltm} and bootstrap estimates from \code{bs.hmltm} 
##' to produce a table ready for insertion into a paper. Does not include any strata with no
##' detections.
##' 
##' @param est output from \code{est.hmltm}.
##' @param bs output  from \code{bootsum} using the same model that created \code{est}.
##' 
# estable=function(est,bs,keepzeros=TRUE){
#   est=est$point$ests
#   nonzeros=which(est$n>0)
#   if(!keepzeros) est=est[nonzeros,]
#   # replace stratum.Area with p
#   trow=length(est[,1])
#   Nc=est$Ngroups[trow]*sum(est$covered.area[-trow])/sum(est$stratum.Area[-trow])
#   p=est$n/(est$covered.area*est$Dgroups)
#   p[trow]=est$n[trow]/Nc
#   nL=signif(est$n/est$L,3)
#   cv.nL=round(100*bs$cv[,"n/L"],1)
#   cv.p=round(100*bs$cv[,"p"],1)
#   cv.Ngrp=round(100*bs$cv[,"Ngroups"],1)
#   cv.E.s=round(100*bs$cv[,"mean.size"],1)
#   cv.N=round(100*bs$cv[,"N"],1)
#   lcl.nL=bs$lower[,"n/L"]
#   lcl.p=bs$lower[,"p"]
#   lcl.Ngrp=bs$lower[,"Ngroups"]
#   lcl.E.s=bs$lower[,"mean.size"]
#   lcl.N=bs$lower[,"N"]
#   ucl.nL=bs$upper[,"n/L"]
#   ucl.p=bs$upper[,"p"]
#   ucl.Ngrp=bs$upper[,"Ngroups"]
#   ucl.E.s=bs$upper[,"mean.size"]
#   ucl.N=bs$upper[,"N"]
#   ci.N=paste("(",round(lcl.N),"; ",round(ucl.N),")",sep="")
#   ci.E.s=paste("(",signif(lcl.E.s,3),"; ",signif(ucl.E.s,3),")",sep="")
#   ci.Ngrp=paste("(",round(lcl.Ngrp),"; ",round(ucl.Ngrp),")",sep="")
#   out=data.frame(Strat=est$stratum,A=round(est$stratum.Area),n=est$n,L=round(est$L,1),
#                  a=round(est$covered.area,1),nL=signif(nL,3),cv.nl=cv.nL,p=signif(p,3),cv.p=cv.p,
#                  Ngrp=round(est$Ngroups),cv.Ngrp=cv.Ngrp,ci.Ngrp=ci.Ngrp,
#                  E.s=signif(est$mean.size,3),cv.E.s=cv.E.s,ci.E.s=ci.E.s,
#                  N=round(est$N),cv.N=cv.N,ci.N=ci.N)
#   return(out)
# }




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
