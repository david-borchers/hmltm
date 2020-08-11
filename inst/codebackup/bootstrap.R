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
