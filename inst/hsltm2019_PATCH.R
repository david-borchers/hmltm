# ===============================================================================
# This is a massive patch for package hsltm, to deal with density estimation only
# in some strata. Need it for bowhead analysis December 2018, because can't get 
# get hsltm to compile.
#
# The patch contains more functions than are needed - I just included all 
# functions in any file for which a patch to at least one function was needed, 
# so it contains replicates of hsltm functions that did not need to be patched.
# ===============================================================================

# ====================================
# Patch from file density_abundance.R
# ====================================
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
#' @param Dstrat vector containing the numbers of strata in which density is to be estimated 
#' (must be a subset of \code{unique(dat$stratum)}.)
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
#' 
#' @export
#' @useDynLib hsltm
est.hmltm <- function(dat,pars,FUN,models=list(y=NULL,x=NULL),survey.pars,hmm.pars,
                      control.fit,control.opt,twosit=FALSE,notrunc=FALSE,W.est=NULL,
                      groupfromy=NULL,Dstrat=unique(dat$stratum))
{
  ## # convert mrds data to cds data format (combining detections)
  ## if(twosit) 
  ##   dat1 <- make.onesit(dat) 
  ## else 
  ##   dat1 <- dat
  
  if(twosit) {
    # add the observer's id to the x-model, and to the y-model if necessary
    if(!is.null(models$x))
      models$x <- update(models$x, ~ . + id)
    else
      models$x <- ~id
    
    if(!is.null(models$y))
      models$y <- update(models$y, ~ . + id)
    else if(!is.null(dat$y))
      models$y <- ~id
  }
  
  # data truncation:
  if(!notrunc) {
    if(min(dat$x,na.rm=TRUE) < survey.pars$Wl | survey.pars$W < max(dat$x,na.rm=TRUE))
      dat <- truncdat(dat,minval=survey.pars$Wl,maxval=survey.pars$W,twosit=FALSE)
  }
  
  srows <- !is.na(dat$object)
  sdat <- dat[srows,]
  
  # fit the model:
  hmltm.fit <- fit.hmltm(sdat,pars,FUN,models,survey.pars,hmm.pars,control.fit,control.opt,
                         groupfromy=groupfromy)
  
  # estimate density, etc:
  if(is.null(W.est)) 
    W.est <- survey.pars$W
  
  if(!is.null(survey.pars$Wl)) 
    W.est <- survey.pars$W-survey.pars$Wl # since in this case data all shifted left by $Wl
  
  point <- NDest(dat,hmltm.fit,W.est,Dstrat,survey.pars)
  hmltm.obj <- list(hmltm.fit=hmltm.fit,point=point,dat=dat,W.est=W.est)
  class(hmltm.obj) <- c(class(hmltm.obj),"hmltm")
  
  return(hmltm.obj)
}

#' @title Estimates density and abundance.
#'
#' @description
#' Horvitz-Thompson like estimation of density and abundance of groups and of individuals, as well as
#' of group size (estimated as the ratio of individual density and group density estimates). Produces 
#' estimates by stratum and over all strata.
#' This version allows specification of strata to be used for density and abundance estimation
#' 
#' @param dat MRDS data frame. Must have cols "stratum","area","transect","L","size","object","x","y" 
#' (and possibly others).
#' @param hmltm.fit output from \code{\link{fit.hmltm}}.
#' @param W perpendicular truncation distance for estimation.
#' @param survey.pars hmltm survey.pars object (a list)
#' 
#' @export
NDest <- function(dat,hmltm.fit,W,Dstrat=Dstrat,survey.pars){
  
  if(inherits(hmltm.fit,"hmmlt")) fit = hmltm.fit
  else if(inherits(hmltm.fit,"hmltm")) fit = hmltm.fit$hmltm.fit
  else stop("Invalid hmltm.fit object, must be class `hmmlt` or `hmltm`.")
  
  dat <- truncdat(dat,minval=survey.pars$Wl,maxval=survey.pars$W,twosit=FALSE)
  maxx <- max(na.omit(dat$x))
  
  if(maxx>W) {
    cat("Maximum perp. dist=",maxx,"is greater than W=",W,"\n")
    stop("You need a bigger W.")
  }
  
  # Add 1/p column
  dat$invp <- rep(NA,dim(dat)[1])
  invp <- fitted.invp1(fit,W=W)
  
  for(i in 1:length(invp$object)) {
    row <- which(dat$stratum==invp$stratum[i] & dat$transect==invp$transect[i] & dat$object==invp$object[i])
    
    if(length(row)>1) {
      cat("Target stratum:",invp$stratum[i],"\n")
      cat("Target transect:",invp$transect[i],"\n")
      cat("Target sighting:",invp$object[i],"\n")
      cat("Found >1: at rows",row,"\n")
      stop("")
    }
    
    dat$invp[row] <- invp$invp[i]
  }
  
  # Reduce data to data for selected strata:
  dat = dat[is.element(dat$stratum,Dstrat),]
  
  # Calculate density and abundance by stratum
  m2km <- 1/1000
  strat <- unique(dat$stratum)
  nstrat <- length(strat)
  n <- L <- a <- A <- Dg <- D <- Ng <- N <- sbar <- rep(0,nstrat+1)
  stratname <- rep("",nstrat+1)
  
  for(i in 1:nstrat){
    stratname[i] <- as.character(strat[i])
    vdat <- dat[dat$stratum==strat[i],]
    trans <- unique(vdat$transect)
    L.tr <- 0
    
    for(tr in 1:length(trans))
      L.tr <- L.tr+vdat$L[min(which(vdat$transect==trans[tr]))]
    
    L[i] <- L.tr
    a[i] <- L[i]*2*W*m2km
    A[i] <- vdat$area[1]
    svdat <- vdat[!is.na(vdat$object),]
    n[i] <- length(svdat$invp)
    Dg[i] <- sum(svdat$invp)/a[i]
    D[i] <- sum(svdat$size*svdat$invp)/a[i]
    sbar[i] <- D[i]/Dg[i]
    Ng[i] <- Dg[i]*A[i]
    N[i] <- D[i]*A[i]
  }
  
  stratname[nstrat+1] <- "Total"
  Ng[nstrat+1] <- sum(Ng[1:nstrat])
  N[nstrat+1] <- sum(N[1:nstrat])
  A[nstrat+1] <- sum(A[1:nstrat])
  Dg[nstrat+1] <- Ng[nstrat+1]/sum(A[1:nstrat])
  D[nstrat+1] <- N[nstrat+1]/sum(A[1:nstrat])
  n[nstrat+1] <- sum(n[1:nstrat])
  L[nstrat+1] <- sum(L[1:nstrat])
  a[nstrat+1] <- sum(a[1:nstrat])
  sbar[nstrat+1] <- D[nstrat+1]/Dg[nstrat+1]
  
  # add transect frequency:
  tfreq <- apply(table(dat$stratum,dat$transect)>0,1,sum)
  k <- c(tfreq,sum(tfreq))
  
  return(list(invp=invp,
              ests=data.frame(stratum=stratname,n=n,k=k,L=L,covered.area=a,stratum.Area=A,
                              Dgroups=signif(Dg,3),Ngroups=signif(Ng,3),mean.size=round(sbar,1),
                              D=signif(D,5),N=round(N,1))
  )
  )
}


# 11 Aug 2020 Patched to here

# ====================================
# Patch from file bootstrap.R
# ====================================
#                          ---------------------
#-------------------------- Bootstrap functions ------------------------


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
#' @return If \code{ests} is NULL, returns a list with elements as follows for each statistic 
#' in \code{bests}:
#' \itemize{
#'  \item{nbad} {number of bad estimates of the statistic in question that were excluded from the 
#'  summary. (Bad estimates occur, for example, when a bootstrap sample involves no detections and
#'  the estimation function tries to calculate mean group size.)}
#'  \item{means} {Bootstrap means of the statistic in question.}
#'  \item{cv} {Percentage CVs of the statistic in question.}
#'  \item{se} {SEs of the statistic in question.}
#'  \item{lower} {Lower \code{cilevel} percentiles of the statistic in question.}
#'  \item{upper} {Upper \code{cilevel} percentiles of the statistic in question.}
#' }
#' If \code{ests} is an object of class 'hmltm', returns a data frame with rows only for every 
#' stratum with detections (and the total) and columns as follows:
#' \itemize{
#' \item{Stratum} stratum number
#' \item{n} original number of detections
#' \item{n.L} original encounter rate
#' \item{CV.n.L} percentage coefficient of variation of encounter rate
#' \item{N.grp} original estimate of group abundance
#' \item{CV.N.grp} percentage coefficient of variation of group abundance
#' \item{N.grp.lo} lower bound of group abundance confidence interval
#' \item{N.grp.hi} upper bound of group abundance confidence interval
#' \item{Es} original mean group size estimate
#' \item{CV.Es} percentage coefficient of variation of mean group size
#' \item{Es.lo} lower bound of mean group size confidence interval
#' \item{Es.hi} upper bound of mean group size confidence interval
#' \item{N} original individual abundance estimate
#' \item{CV.N} percentage coefficient of individual abundance
#' \item{N.lo} lower bound of individual abundance confidence interval
#' \item{N.hi} upper bound of individual abundance confidence interval
#'}
#'
#' @export
bootsum <- function(bests,ests=NULL,cilevel=0.95,write.csvs=FALSE,dir=getwd()){
  if(cilevel<=0 | cilevel>=1) 
    stop("cilevel must be greater than 0 and less than 1.")
  if(!is.null(ests) & !inherits(ests,"hmltm")) 
    stop("ests must be of class 'hmltm'.")
  
  ns <- apply(bests[,2,],1,sum) # sum sample sizes by stratum
  keepstrat <- (ns>0) # only consider strata with some detections
  
  if(!is.null(ests)) 
    keepcols <- c("stratum","n","L","Ngroups","mean.size","N")
  else 
    keepcols <- colnames(bests)
  
  if(!is.null(ests)) 
    bsests <- bests[keepstrat,keepcols,]
  else 
    bsests <- bests[keepstrat,,]
  
  # replace L with encounter rate
  Lcol <- which(colnames(bests[,,1])=="L")
  bsests[,Lcol,] <- bsests[,"n",]/bsests[,"L",] 
  colnames(bsests)[Lcol] <- "n/L"
  
  # replace stratum.Area with p
  Acol <- which(colnames(bests[,,1])=="stratum.Area")
  bsests[,Acol,] <- bsests[,"n",]/(bsests[,"covered.area",]*bsests[,"Dgroups",]) 
  colnames(bsests)[Acol] <- "p"
  bdim <- dim(bsests)
  nstrat <- bdim[1]
  nests <- bdim[2]
  cv <- matrix(rep(NA,nstrat*nests),ncol=bdim[2])
  rownames <- dimnames(bsests)[[1]]
  rownames[length(rownames)] <- "Total"
  colnames <- dimnames(bsests)[[2]]
  dimnames(cv) <- list(rep("",length(rownames)),colnames)
  nbad <- se <- means <- lower <- upper <- cv
  B <- bdim[3]
  
  cat("Results from ",B," bootstrap replicates:\n",sep="")
  cat("----------------------------------------\n")
  
  for(i in 1:nstrat) {
    for(j in 1:nests){
      goodests <- na.omit(bsests[i,j,])
      nbad[i,j] <- B-length(goodests)
      means[i,j] <- mean(goodests)
      se[i,j] <- sd(goodests)
      cv[i,j] <- se[i,j]/means[i,j]
      perc <- quantile(goodests,probs = c((1-cilevel)/2,(1-(1-cilevel)/2)))
      lower[i,j] <- perc[1]
      upper[i,j] <- perc[2]
    }
  }
  
  # remove stats of stratum name, and make stratum a character
  nbad <- as.data.frame(nbad,row.names=1:dim(nbad)[1])
  means <- as.data.frame(means,row.names=1:dim(nbad)[1])
  se <- as.data.frame(se,row.names=1:dim(nbad)[1])
  cv <- as.data.frame(cv,row.names=1:dim(nbad)[1])
  lower <- as.data.frame(lower,row.names=1:dim(nbad)[1])
  upper <- as.data.frame(upper,row.names=1:dim(nbad)[1])
  nbad[,1] <- rownames
  means[,1] <- rownames
  se[,1] <- rownames
  cv[,1] <- rownames
  lower[,1] <- rownames
  upper[,1] <- rownames
  
  outsum <- list(nbad=nbad,mean=means,se=se,cv=cv,lower=lower,upper=upper)
  
  if(!is.null(ests)) 
    outsum <- strat.estable(ests$point$ests[keepstrat,],outsum$cv)
  
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
#' @return Returns a data frame with columns as follows:
#' \item{Stratum}{stratum number}
#' \item{n}{original number of detections}
#' \item{n.L}{original encounter rate}
#' \item{CV.n.L}{percentage coefficient of variation of encounter rate}
#' \item{N.grp}{original estimate of group abundance}
#' \item{CV.N.grp}{percentage coefficient of variation of group abundance}
#' \item{N.grp.lo}{lower bound of group abundance confidence interval}
#' \item{N.grp.hi}{upper bound of group abundance confidence interval}
#' \item{E.s}{original mean group size estimate}
#' \item{CV.E.s}{percentage coefficient of variation of mean group size}
#' \item{N}{original individual abundance estimate}
#' \item{CV.N}{percentage coefficient of individual abundance}
#' \item{N.lo}{lower bound of individual abundance confidence interval}
#' \item{N.hi}{upper bound of individual abundance confidence interval}
#' 
#' @export
strat.estable <- function(est,cv){
  min <- est[,"n"]
  N <- est[,"N"]
  Ngrp <- est[,"Ngroups"]
  Es <- est[,"mean.size"]
  cv.N <- cv[,"N"]
  cv.Ngp <- cv[,"Ngroups"]
  cv.Es <- cv[,"mean.size"]
  N.ci <- lnci.nmin(min,N,cv.N)
  Ngrp.ci <- lnci.nmin(min,Ngrp,cv.N)
  Es.ci <- lnci.nmin(0,Es,cv.Es)  
  
  outp <- data.frame(Stratum=est$stratum,
                     n=est$n,
                     n.L=signif(est$n/est$L,3),
                     CV.n.L=round(cv[,"n/L"]*100),
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
                     N.hi=round(N.ci$upper))
  
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
#' @export
cv.avail <- function(hmm.pars,B=1000){
  n <- dim(hmm.pars$Et)[2]
  b.Et <- array(rep(NA,B*2*n),dim=c(B,2,n),dimnames=list(1:B,State=c("Unavailable","Available"),Animal=1:n))
  cv.a <- cv.u <- rep(NA,n)
  
  for(i in 1:n){
    # get normal parameters corresponding to lognormal(mu,Sigma)
    lN <- Npars.from.lNpars(hmm.pars$Et[,i],hmm.pars$Sigma.Et[,,i])
    
    # resample availability parameters on lognormal scale
    b.Et[,,i] <- exp(mvrnorm(B,lN$mu,lN$Sigma)) 
    
    a <- b.Et[,1,i]/(b.Et[,1,i]+b.Et[,2,i])
    u <- b.Et[,2,i]/(b.Et[,1,i]+b.Et[,2,i])
    cv.a[i] <- cv(a)
    cv.u[i] <- cv(u)
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
#' @export
bs.hmltm=function(hmltm.est,B,hmm.pars.bs=NULL,bs.trace=0,report.by=10,fixed.avail=FALSE){
  #######################################################
  ## Produces 3-dim array containing B sets of density ##
  ## and abundance estimates from data extraction      ##
  #######################################################
  
  dat1 <- hmltm.est$dat # extract original data frame from fitted object
  W.est <- hmltm.est$W.est # extract estimation perp dist truncation
  survey.pars <- hmltm.est$hmltm.fit$fitpars$survey.pars
  hmm.pars <- hmltm.est$hmltm.fit$fitpars$hmm.pars
  control.fit <- hmltm.est$hmltm.fit$fitpars$control.fit
  control.opt <- hmltm.est$hmltm.fit$fitpars$control.opt
  control.opt$trace <- bs.trace
  h.fun <- hmltm.est$hmltm.fit$h.fun
  models <- hmltm.est$hmltm.fit$models
  pars <- hmltm.est$hmltm.fit$fit$par
  
  if(!fixed.avail) {
    # bootstrap availability parameters
    if(!is.null(hmm.pars$Et)){
      cat("Bootstrap with parametric resampling of mean times available and unavailable.\n")
      flush.console()
      
      n <- dim(hmm.pars$Et)[2]
      ns <- bsample(1:n,n,replace=TRUE) # sample animals to use
      b.Et <- array(rep(NA,B*2*n),dim=c(B,2,n),
                    dimnames=list(1:B,State=c("Unavailable","Available"),Animal=1:n))
      
      for(i in ns){
        # get normal parameters corresponding to lognormal(mu,Sigma)
        lN <- Npars.from.lNpars(hmm.pars$Et[,i],hmm.pars$Sigma.Et[,,i]) 
        
        # resample availability parameters on lognormal scale
        b.Et[,,i] <- exp(mvrnorm(B,lN$mu,lN$Sigma)) 
      }
    } else if(!is.null(hmm.pars.bs)) { # resample availability parameters
      cat("Bootstrap with parametric resampling of HMM parameters.\n")
      flush.console()
      
      nhmm <- dim(hmm.pars.bs)[1]
      
      if(nhmm==1) 
        stop(paste("Only one set of hmm pars. need multiple sets of pars if not fixed avail ",
                   "(i.e. if hmm.pars.bs!=NULL)"))
      
      reps <- bsample(1:nhmm,B,replace=TRUE)
    } else {
      warning("No availability parameters to resample from. Treating availability parameters as constant.")
      flush.console()
      fixed.avail <- TRUE
    }
  } else {
    cat("Bootstrap treating HMM parameters as constant.\n")
    flush.console()
  }
  
  # get stuff needed to set up bootstrap datasets
  strat <- unique(dat1$stratum) # unique stratum names
  nrows <- table(dat1$stratum,dat1$transect) # number of rows required for each transect
  get.ntrans <- function(trans) sum(trans>0)
  ntrans <- apply(nrows,1,get.ntrans) # number of transects in each stratum
  nstrat <- length(ntrans) # number of strata
  transects <- as.numeric(colnames(nrows)) # vector of transect numbers
  
  # create bootstrap data, estimate, and store estimates
  nDstrat = length(hmltm.est$point$ests$stratum)-1 # Number of strata for which we want estimates
  Dstrat = hmltm.est$point$ests$stratum[1:nDstrat] # Strata for which we want estimates
  estdim <- dim(hmltm.est$point$ests) #### bitchange
  bestdim <- c(estdim,B) #### bitchange
  bestdimnames <- list(as.character(hmltm.est$point$ests$stratum),colnames(hmltm.est$point$ests),1:B)
  best <- array(rep(NA,B*prod(estdim)),dim=bestdim,dimnames=bestdimnames)
  bn <- rep(0,B) #### bitchange
  
  for(b in 1:B) { # do bootstraps
    # get hmm pars
    if(!fixed.avail) {
      if(!is.null(hmm.pars$Et)) {
        Pi <- makePi(b.Et[b,1,],b.Et[b,2,]) # make Pi from resampled times
        delta <- matrix(rep(NA,2*n),nrow=2)
        
        for(i in 1:n) 
          delta[,i] <- compdelta(Pi[,,i])
        
        b.hmm.pars <- hmm.pars
        b.hmm.pars$pm <- hmm.pars$pm
        b.hmm.pars$Pi <- Pi
        b.hmm.pars$delta <- delta
        b.hmm.pars$Et <- b.Et[b,,]
      } else {
        b.hmm.pars <- unvectorize.hmmpars(hmm.pars.bs[reps[b],]) # resampled HMM pars
      }
    } else {
      b.hmm.pars <- hmm.pars # fixed availability pars
    }
    
    # to fix naming cock-up when creating hmmpars.bs
    if(is.element("pm",names(b.hmm.pars))) 
      names(b.hmm.pars)[which(names(b.hmm.pars)=="pm")] <- "pm"
    
    # resample transects with replacement
    newtransind <- matrix(rep(NA,nstrat*max(ntrans)),nrow=nstrat) # matrix of indices for resampled transects
    for(st in 1:nstrat) {
      # indices of matrix nrows for resampled transects
      newtransind[st,1:ntrans[st]] <- bsample(which(nrows[st,]>0),ntrans[st],replace=TRUE)
    }
    
    # calc number rows for bootstrap data frame
    newnrows <- 0
    for(st in 1:nstrat) 
      newnrows <- newnrows+sum(na.omit(nrows[st,newtransind[st,]]))
    
    # set up empty bootstrap data frame
    bdat1 <- data.frame(matrix(rep(NA,dim(dat1)[2]*newnrows),nrow=newnrows))
    names(bdat1) <- names(dat1)
    start <- 1
    tno <- 1 # fill in new data frame from old:
    
    for(st in 1:nstrat) 
      for(tr in 1:ntrans[st]) {
        addtrans <- transects[newtransind[st,tr]]
        nadd <- nrows[st,newtransind[st,tr]]
        newi <- start:(start+nadd-1)
        oldi <- which(dat1$stratum==strat[st] & dat1$transect==addtrans)
        bdat1[newi,] <- dat1[oldi,]
        bdat1$transect[newi] <- tno
        start <- start+nadd
        tno <- tno+1
      }
    
    bn[b] <- length(bdat1$object[!is.na(bdat1$object)])
    
    if(b==1) 
      cat("Sample sizes: ")
    
    if(b%%report.by==0) {
      cat(paste(bn[b],";",sep=""),"Iterations done:",b,"\n")
      if(b!=B) 
        cat("Sample sizes: ")
    }
    else 
      cat(paste(bn[b],";",sep=""))
    
    bhat <- est.hmltm(bdat1,pars=pars,FUN=h.fun,models=models,survey.pars,b.hmm.pars,control.fit,
                      control.opt,notrunc=TRUE,W.est=W.est,Dstrat=Dstrat) 
    
    for(st in 1:(nDstrat+1)) 
      best[st,,b] <- as.numeric(bhat$point$ests[st,])
  }
  
  class(best) <- c(class(best),"hmltm.bs")
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
#' @export
Npars.from.lNpars <- function(mu,Sigma){
  cv2 <- diag(Sigma)/mu^2 # Sigma is variance matrix
  logvar <- log(cv2+1)
  logmu <- log(mu)-logvar/2
  logSigma <- Sigma*0
  np <- length(mu)
  
  for(i in 1:(np-1)) {
    logSigma[i,i] <- logvar[i] # variance of normal
    for(j in (i+1):np) {
      # covariance of normal
      logSigma[i,j] <- log(Sigma[i,j]/(exp(logmu[i]+logmu[j]+(logvar[i]+logvar[j])/2))+1) 
      logSigma[j,i] <- logSigma[i,j]
    } 
  }
  
  logSigma[np,np] <- logvar[np] # variance of normal
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
#' @importFrom HiddenMarkov dthmm
#' @importFrom HiddenMarkov BaumWelch
#' 
#' @export
resample.hmmpars <- function(availhmm,adat,animals,seed=NULL,nperow=10,printprog=TRUE){
  b.availhmm <- availhmm # initialise bootstrap HMM parameters
  
  # filter out animals not chosen
  b.availhmm$Pi <- b.availhmm$Pi[,,animals]
  b.availhmm$pm <- b.availhmm$pm[,animals]
  b.availhmm$delta <- b.availhmm$delta[,animals]
  
  # simulate and fit new HMMs
  newanimal <- 1
  for(animal in animals) {
    hmmobj <- dthmm(adat[[animal]],availhmm$Pi[,,animal],availhmm$delta[,animal],"binom",
                    pm=list(prob=availhmm$pm[,animal]),pn=list(size=rep(1,length(adat[[animal]]))),
                    discrete = TRUE)
    
    sdat <- simulate(hmmobj, nsim=length(hmmobj$x),seed=seed)
    fitanimal <- BaumWelch(sdat,control=list(maxiter=500,tol=1e-05,prt=FALSE,posdiff=TRUE,
                                             converge = expression(diff < tol)))
    
    b.availhmm$pm[,newanimal] <- fitanimal$pm$prob
    b.availhmm$Pi[,,newanimal] <- fitanimal$Pi
    b.availhmm$delta[,newanimal] <- compdelta(fitanimal$Pi)
    
    if(printprog) {
      if(animal==1) 
        cat("Individuals done:")
      if((animal%/%nperow)*nperow==animal) # got an exact multiple of 10
        if(animal>1) 
          cat(": Total=",animal,"\nIndividuals done:")
    }
    
    cat(" *")
    newanimal <- newanimal+1
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
#' 
#' @export
vectorize.hmmpars <- function(hmmpars) {
  Pi <- hmmpars$Pi
  pm <- hmmpars$pm
  delta <- hmmpars$delta
  dims <- dim(Pi)
  
  if(length(dims)!=3) 
    stop("Pi must be a 3-D array")
  if(dims[1]!=dims[2]) 
    stop("First two dimensions of Pi unequal")
  if(dim(pm)[1]!=dims[1] | dim(pm)[2]!=dims[3]) 
    stop("dim(pm) inconsistent with dim(Pi)")
  if(dim(delta)[1]!=dims[1] | dim(delta)[2]!=dims[3]) 
    stop("dim(delta) inconsistent with dim(Pi)")
  
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
#' 
#' @export
unvectorize.hmmpars <- function(hv) {
  m3d <- hv[1:3]
  m2d <- c(m3d[1],m3d[3])
  m3size <- prod(m3d)
  m2size <- prod(m2d)
  Pi <- array(hv[(3+1):(3+m3size)],dim=m3d)
  pm <- array(hv[(3+m3size+1):(3+m3size+m2size)],dim=m2d)
  delta <- array(hv[(3+m3size+m2size+1):(3+m3size+m2size+m2size)],dim=m2d)
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
#' @export
hmmpars.boot <- function(availhmm,adat,animals,seed=NULL,B,printprog=TRUE){
  # initialise matrix of correct dimensions:
  na <- length(adat)
  #  ncol=1
  #  for(i in 1:na) if(length(adat[[i]])>ncol) ncol=length(adat[[i]])
  
  onelength <- length(vectorize.hmmpars(list(Pi=availhmm$Pi[,,1,drop=FALSE],
                                             pm=availhmm$pm[,1,drop=FALSE],
                                             delta=availhmm$delta[,1,drop=FALSE])))
  
  ncol <- 3+(onelength-3)*length(animals) # 3 numbers specify vectorized length parameters
  hmmat <- matrix(rep(NA,ncol*B),ncol=ncol)
  
  # do the bootstrapping
  for(b in 1:B) {
    hmp <- resample.hmmpars(availhmm,adat,animals,seed=NULL)
    hmmvec <- vectorize.hmmpars(hmp)
    #    nobs=length(hmmvec)
    #    hmmat[b,nobs]=hmmvec
    hmmat[b,] <- hmmvec
    
    if(printprog) 
      cat("Resamples done: ",b,"\n")
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
#' @importFrom MASS mvrnorm
#' @importFrom HiddenMarkov compdelta
#' 
#' @export
bootstrap.p.with.Et <- function(dat,pars,hfun,models,survey.pars,hmm.pars,
                                control.fit,control.opt,fixed.avail=FALSE,B=999){
  n <- length(dat$x)
  npar <- length(pars)
  b.p0 <- b.phat <- rep(NA,B)
  b.pars <- matrix(rep(NA,B*npar),ncol=npar)
  b.phats <- matrix(rep(NA,B*n),ncol=n)# matrix for detection probs of each individual.
  
  # bootstrap availability parameters
  if(!fixed.avail) {
    # get normal parameters corresponding to lognormal(mu,Sigma)
    lN <- Npars.from.lNpars(hmm.pars$Et,hmm.pars$Sigma.Et) 
    
    # resample availability parameters on lognormal scale
    b.Et <- exp(mvrnorm(B,lN$mu,lN$Sigma)) 
  }
  
  # resample sightings data and re-estimate, using a resample of availability paramters:
  for(nb in 1:B){
    # resample detection locations
    samp.ind <- bsample(1:n,size=n,replace=TRUE) # resample sightings data indices with replacement
    b.dat <- dat[samp.ind,]# get resampled data
    
    # create new hmm.pars object with resampled availability parameters
    if(!fixed.avail) {
      pi21 <- 1/b.Et[nb,2]
      pi12 <- 1/b.Et[nb,1]
      Pi <- matrix(c((1-pi12),pi12,pi21,(1-pi21)),nrow=2,byrow=TRUE)
      delta <- compdelta(Pi)
      pm <- c(0.0,1.0)
      b.hmm.pars <- list(pm=pm,Pi=Pi,delta=delta,b.Et[nb],hmm.pars$Sigma.Et)
    }else {
      b.hmm.pars <- hmm.pars
    }
    
    # refit model
    b.fit <- fit.hmltm(b.dat,pars=pars,FUN=hfun,models=models,survey.pars=survey.pars,
                       hmm.pars=b.hmm.pars,control.fit=control.fit,control.optim=control.opt)
    
    b.p0[nb] <- b.fit$pzero
    b.phat[nb] <- b.fit$phat
    b.phats[nb,] <- b.fit$phats
    b.pars[nb,] <- b.fit$fit$par
    
    cat(paste("done",nb,"\n"))
    flush.console()
  }
  
  # package results and return
  callist <- list(dat=dat,pars=pars,hmm.pars=hmm.pars,hfun=hfun,models=models,
                  survey.pars=survey.pars,control.fit=control.fit,control.optim=control.opt)
  bs <- list(phats=b.phats,pars=b.pars,p0=b.p0,phat=b.phat,b.Et=b.Et)
  
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
#' @export
bootstrap.p.with.hmm <- function(dat,pars,hfun,models,survey.pars,hmm.pars.bs,
                                 control.fit,control.opt,fixed.avail=FALSE,B=999,silent=FALSE){
  n <- length(dat$x)
  npar <- length(pars)
  conv <- b.p0 <- b.phat <- rep(NA,B)
  b.pars <- matrix(rep(NA,B*npar),ncol=npar)
  b.phats <- matrix(rep(NA,B*n),ncol=n) # matrix for detection probs of each individual.
  
  # bootstrap availability parameters
  if(!fixed.avail) { # resample availability parameters
    nhmm <- dim(hmm.pars.bs)[1]
    if(nhmm==1) stop("Only one set of hmm pars. need multiple sets of pars if not fixed avail.")
    reps <- bsample(1:nhmm,B,replace=TRUE)
  }
  for(nb in 1:B) {
    if(!fixed.avail)
      b.hmm.pars <- unvectorize.hmmpars(hmm.pars.bs[reps[nb],])
    else
      b.hmm.pars <- hmm.pars.bs
    
    # to fix naming cock-up when creating hmmpars.bs
    if(is.element("pm",names(b.hmm.pars))) 
      names(b.hmm.pars)[which(names(b.hmm.pars)=="pm")] <- "pm" 
    
    # resample detection locations
    samp.ind <- bsample(1:n,size=n,replace=TRUE) # resample sightings data indices with replacement
    b.dat <- dat[samp.ind,,drop=FALSE]# get resampled data
    names(b.dat) <- names(dat)
    
    # refit model
    b.fit <- try(fit.hmltm(b.dat,pars=pars,FUN=hfun,models=models,survey.pars=survey.pars,
                           hmm.pars=b.hmm.pars,control.fit=control.fit,control.optim=control.opt),
                 silent=silent)
    
    if((class(b.fit)=="try-error")) {
      conv[nb] <- -999
      b.p0[nb] <- -999
      b.phat[nb] <- -999
      b.phats[nb,] <- rep(-999,n)
      b.pars[nb,] <- rep(-999,length(pars))
    } else {
      conv[nb] <- b.fit$fit$convergence
      b.p0[nb] <- b.fit$p[1]
      b.phat[nb] <- b.fit$phat
      b.phats[nb,] <- b.fit$phats
      b.pars[nb,] <- b.fit$fit$par
      cat(paste("done",nb,"\n"))
    }
    flush.console()
  }
  
  # package results and return
  callist <- list(dat=dat,pars=pars,hmm.pars=hmm.pars.bs,hfun=hfun,survey.pars=survey.pars,
                  control.fit=control.fit,control.optim=control.opt)
  
  bs <- list(phats=b.phats,pars=b.pars,p0=b.p0,phat=b.phat,convergence=conv)
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
#' @export
bootsum.p <- function(bs,probs=c(0.025,0.975),pcut=0){
  ################################################################################
  ## bs is output from bootstrap.p.with.Et() or bootstrap.with.hmm()            ##
  ## pcut is a quick and dirty min phat to allow - robustifies 1/phat for small ##
  ## samples, although it is ad-hoc. Do hist of $bs$phat to see if there is a   ##
  ## reasonable cutpoint.                                                       ##
  ################################################################################
  
  nboot <- length(bs$bs$p0)
  
  if(is.null(bs$bs$convergence)) 
    keep <- which(bs$bs$p0>=0 & bs$bs$phat>pcut)
  else 
    keep <- which(bs$bs$p0>=0 & bs$bs$convergence==0 & bs$bs$phat>pcut)
  
  nbad <- nboot-length(keep)
  npar <- dim(bs$bs$par)[2]
  cinames <- paste(as.character(probs*100),"%",sep="")  
  bests <- matrix(rep(NA,(3+npar)*5),ncol=5)
  colnames(bests) <- c("mean","std.err.","%CV",cinames)
  parnames <- paste("par",as.character(1:npar),sep="")
  rownames(bests) <- c("1/phat","phat","p(0)",parnames)
  
  # relative density
  bests[1,1] <- mean(1/bs$bs$phat[keep])
  bests[1,2] <- sd(1/bs$bs$phat[keep])
  bests[1,3] <- sd(1/bs$bs$phat[keep])/mean(1/bs$bs$phat[keep])*100
  bests[1,4:5] <- quantile(1/bs$bs$phat[keep],probs=probs)
  
  # detection probability
  bests[2,1] <- mean(bs$bs$phat[keep])
  bests[2,2] <- sd(bs$bs$phat[keep])
  bests[2,3] <- sd(bs$bs$phat[keep])/mean(bs$bs$phat[keep])*100
  bests[2,4:5] <- quantile(bs$bs$phat[keep],probs=probs)
  
  # p(0)
  bests[3,1] <- mean(bs$bs$p0[keep])
  bests[3,2] <- sd(bs$bs$p0[keep])
  bests[3,3] <- sd(bs$bs$p0[keep])/mean(bs$bs$p0[keep])*100
  bests[3,4:5] <- quantile(bs$bs$p0[keep],probs=probs)
  
  # parameters
  for(i in 1:npar) {
    bests[3+i,1] <- sd(bs$bs$par[keep,i])
    bests[3+i,2] <- mean(bs$bs$par[keep,i])
    bests[3+i,3] <- sd(bs$bs$par[keep,i])/mean(bs$bs$par[keep,i])*100
    bests[3+i,4:5] <- quantile(bs$bs$par[keep,i],probs=probs)
  }
  
  parcov <- cov(bs$bs$par)
  parcorr <- cov2cor(parcov)
  rownames(parcov) <- rownames(parcorr) <- colnames(parcov) <- colnames(parcorr) <- parnames
  
  return(list(nboot=nboot,nbad=nbad,bests=bests,parcov=parcov,parcorr=parcorr))
}





# ====================================
# Patch from file DLB_utility.R
# ====================================
#                      ------------------------
#--------------------- DLB's utility functions --------------------

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
bsample <- function(x,size,replace=FALSE,prob=NULL) {
  if(length(x)==1) 
    return(x)
  else 
    return(sample(x,size,replace,prob))
}


#' @title Caclulates coefficient of variation.
#'
#' @description
#'  Utility function
#'  
#' @param x Random variable.
cv <- function(x) sd(x)/mean(x) # calculates coefficient of variation

#' @title Converts nm to metres.
#'
#' @description
#'  Utility function
#'  
#' @param x Distance in nautical miles.
nm2m <- function(x) return(x*1852) # converts nautical miles to metres


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
#' 
#' @export
histline <- function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                     xlab="x",ylab="y",...)
{
  n <- length(height)
  if(length(breaks)!=(n+1)) 
    stop("breaks must be 1 longer than height")
  
  if(outline) {
    y <- c(0,rep(height,times=rep(2,n)),0)
    x <- rep(breaks,times=rep(2,(n+1)))
  } else {
    y <- rep(0,4*n)
    x <- rep(0,4*n+2)
    
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)] <- c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)] <- c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x <- x[1:(4*n)]
  }
  
  if(lineonly) {
    if(!fill) 
      lines(x,y,...)
    else 
      polygon(x,y,...)
  } else {
    if(!fill) 
      plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
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
sintegral <- function (fx, x, n.pts=16, type="int") 
{
  #  if (class(fx) == "function") 
  #    fx = fx(x)
  #  n.x = length(x)
  #  if (n.x != length(fx)) 
  #    stop("Unequal input vector lengths")
  #  if (n.pts < 64) 
  #    n.pts = 64
  ap <- approx(x, fx, n = 2 * n.pts + 1)
  h <- diff(ap$x)[1]
  integral <- h*(ap$y[2*(1:n.pts)-1]+4*ap$y[2*(1:n.pts)]+ap$y[2*(1:n.pts)+1])/3
  
  if(type!="cdf") 
    return(sum(integral)) 
  else 
    return(list(x=ap$x[2*(1:n.pts)],y=cumsum(integral)))
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
#' 
#' # because hmm.pars is in units of time, not distance, you need to specify spd.
#' laake.a(hmm.pars,ymax=200,spd=4) 
#' 
#' @export

laake.a <- function(hmm.pars,ymax,spd=NULL){
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a <- rep(NA,nav)
  
  for(i in 1:nav)
    a[i] <- jeffa(hmm.pars$Pi[,,i],ymax,spd)
  
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
#' 
#' @export
jeffa <- function(Pi,w,spd=NULL,E=NULL)
{
  if(!is.null(spd)) 
    w <- w/spd # convert distance to time if spd given
  
  if(!is.null(E))
    warning("Used E, not Pi for calculations")
  else
    E <- makeE(Pi) # expected time up, expected time down
  
  a <- (E[2] + E[1]*(1-exp(-w/E[1])))/sum(E)    
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
#' 
#' @export
instant.a <- function(hmm.pars,Et=NULL){
  if(!is.null(Et)) {
    if(!is.null(hmm.pars)) 
      warning("Used Et, not Pi for calculations")
    
    if(is.vector(Et)) 
      Et <- matrix(Et,ncol=1)
    
    if(dim(Et)[1]!=2) 
      stop("1st dimension of Et must be 2.")
    
    nav=dim(Et)[2]
  } else {
    if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) 
      stop("1st two dimensions of hmm.pars$Pi must be 2.")
    
    if(length(dim(hmm.pars$Pi))==2) 
      hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
    
    nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  }
  
  a <- rep(NA,nav)
  
  for(i in 1:nav) {
    if(is.null(Et)) 
      a[i] <- simplea(hmm.pars$Pi[,,i],Et)
    else 
      a[i] <- simplea(NULL,Et[,i])
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
#' 
#' @export
simplea <- function(Pi,E=NULL)
{
  if(is.null(E)) 
    E <- makeE(Pi) # expected time up, expected time down
  a <- E[2]/sum(E)    
  
  return(a)
}

#' @title McLaren's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates McLaren's availability correction factor, for each of a set of m 2-state Markov model 
#' availability parameters passed to it and returns these and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element \code{\$Pi} (where m is number of availability parameter sets).
#' @param w is max forward distance things can be seen at (or max forward time). Must be scalar.
#' @param spd is observer speed; omit if w is max forward TIME.
#'
#' @references
#' McLaren, I.A. 1961. Methods of determining the numbers and availability of ringed seals in the 
#' eastern Canadian Arctic. Arctic 14:162--175.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' mclaren.a(hmm.pars,w=10,spd=4)
#' mclaren.a(hmm.pars,w=100,spd=4) # can be greater than 1 (!)
#' 
#' @export
mclaren.a <- function(hmm.pars,w,spd=1){
  if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) 
    stop("1st two dimensions of hmm.pars$Pi must be 2.")
  
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a <- rep(NA,nav)
  
  for(i in 1:nav){
    Eti <- makeE(hmm.pars$Pi[,,i])*spd # Pi is TIME; if w is DISTANCE need to multiply time by speed
    a[i] <- (Eti[2]+w)/sum(Eti)
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
#' 
#' @export
richard.a <- function(hmm.pars,w,spd=1){
  y <- na.omit(w)
  n <- length(y)
  
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  ina <- matrix(rep(instant.a(hmm.pars)$a,n),nrow=n,byrow=TRUE)
  mca <- matrix(rep(NA,n*nav),nrow=n)
  
  for(i in 1:n) 
    mca[i,] <- mclaren.a(hmm.pars,y[i],spd)$a
  
  sumfb <- apply(ina/mca,2,sum)
  a <- 1/((1/ina[1,])*(sumfb/n))
  
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
#' 
#' @export
makePi <- function(Eu,Ea)
{
  nav <- length(Eu)
  if(length(Ea)!=nav) 
    stop("Lengths of Eu and Ea must be the same")
  #  if(nav==1) {
  #    Pi=matrix(rep(0,4),nrow=2,
  #              dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available")))
  #    Pi[1,2]=1/Eu
  #    Pi[2,1]=1/Ea
  #    Pi[1,1]=1-Pi[1,2]
  #    Pi[2,2]=1-Pi[2,1]
  #  } else {
  Pi <- array(rep(0,2*2*nav),dim=c(2,2,nav),
              dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                            Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Pi[1,2,i] <- 1/Eu[i]
    Pi[2,1,i] <- 1/Ea[i]
    Pi[1,1,i] <- 1-Pi[1,2,i]
    Pi[2,2,i] <- 1-Pi[2,1,i]    
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
#' 
#' @export
makeE <- function(Pi){
  #----------------------------------------------------------
  # Returns expected time in states 1 and 2 for the 2x2 
  # probability transition matrix Pi for a 2-state
  # Markov process.
  #----------------------------------------------------------
  if(length(dim(Pi))==2){
    E <- c(1/Pi[1,2],1/Pi[2,1])
    names(E) <- c("Unavailable","Available")
  } else {
    nav <- dim(Pi)[3]
    E <- matrix(rep(NA,2*nav),nrow=2,
                dimnames=list(State=c("Unavailable","Available"),Animal=as.character(1:nav)))
    
    for(i in 1:nav)
      E[,i] <- c(1/Pi[1,2,i],1/Pi[2,1,i])
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
#' @export
make.hmm.pars.from.Et <- function(Ea,Eu,seEa,seEu,covEt=0,pm=NULL) {
  nav <- length(Ea)
  
  if(length(Eu)!=nav |length(seEa)!=nav |length(seEu)!=nav |length(covEt)!=nav) 
    stop("Lengths of Ea, Eu, seEa, seEu, covEt must all be the same.")
  
  if(is.null(pm)) 
    pm <- matrix(c(rep(0,nav),rep(1,nav)),nrow=2,byrow=TRUE)
  
  if(is.vector(pm)) {
    if(length(pm)!=2) 
      stop("pm must either be a vector of length 2 or a matrix of dimension length(Ea)x2.")
    
    pm <- matrix(c(pm[1],pm[2]),ncol=2)
  }
  
  if(dim(pm)[2]!=nav) 
    stop("Inconsistent dimensions of Ea and pm.")
  
  Pi <- Sigma.Et <- array(rep(NA,2*2*nav),dim=c(2,2,nav),
                          dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                                        Animal=as.character(1:nav)))
  
  Et <- delta <- newpm <- matrix(rep(NA,2*nav),ncol=nav,dimnames=list(State=c("Unavailable","Available"),
                                                                      Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Et[,i] <- c(Eu[i],Ea[i])
    Sigma.Et[,,i] <- diag(c(seEu[i],seEa[i])^2)
    #    cvEt=c(seEu[i]/Et[1,i],seEa[i]/Et[2,i])
    #    Sigma.Et[,,i]=diag((cvEt*Et)^2)
    Sigma.Et[1,2,i] <- Sigma.Et[2,1,i] <- covEt[i]
    Pi[,,i] <- makePi(Et[1,i],Et[2,i])
    delta[,i] <- compdelta(Pi[,,i])
    newpm[,i] <- pm[,i]
  }
  
  hmm.pars <- list(pm=newpm,Pi=Pi,delta=delta,Et=Et,Sigma.Et=Sigma.Et)
  return(hmm.pars)  
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
#'
#' @details Packs the above in a list suitable for passing as \code{survey.pars} to 
#' \code{\link{est.hmltm}}.
#' 
#' @export
make.survey.pars <- function(spd,W,ymax,Wl=0,dT=1){
  return(list(spd=spd,W=W,ymax=ymax,Wl=Wl,dT=dT,dy=spd*dT))
}

#' @title Decides if model is a null model.
#'
#' @description
#'  Logical function: true if \code{models} includes no covariates
#'  
#' @param models list of characters with elements \code{$y} and \code{$x} specifying y- and 
#' x-covariate models. Either \code{NULL} or regression model format (without response on left).
#' 
#' @export
is.nullmodel <- function(models){
  null <- TRUE
  
  for(i in 1:length(models)) 
    null <- null & is.null(models[[i]])
  
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
#' 
#' @export
make.covb <- function(b,FUN,models,dat)
{
  nfixed <- switch(FUN,
                   h.IP.0 = 2,
                   h.EP1.0 = 2,
                   h.EP2.0 = 3,
                   h.EP1x.0 = 3,
                   h.EP2x.0 = 4,
                   0)
  
  if(nfixed==0) 
    stop("Hazard model ",FUN," is not programmed (yet).")
  
  if(is.nullmodel(models)) {
    if(length(b) != nfixed) 
      stop("length of b inconsistent with model")
    
    n <- dim(dat)[1]
    #    covb=matrix(c(rep(b,rep(n,nfixed))),ncol=nfixed)
    covb <- rep(b,n)
  } else if(FUN=="h.EP2.0" | FUN=="h.EP1.0" | FUN=="h.IP.0") {
    X <- model.matrix(as.formula(models$y),data=dat)
    n <- dim(X)[1]
    nb <- dim(X)[2]
    
    nfixed <- nfixed-1
    if(length(b) != (nfixed+nb)) 
      stop("length of b inconsistent with model")
    
    covb <- matrix(c(rep(b[1:nfixed],rep(n,nfixed)),X%*%b[nfixed+(1:nb)]),ncol=(nfixed+1))
    covb <- as.vector(t(covb))
  } else if(FUN=="h.EP2x.0" | FUN=="h.EP1x.0") {
    if(!is.null(models$y)){
      X <- model.matrix(as.formula(models$y),data=dat)
      n <- dim(X)[1]
      nb <- dim(X)[2]   
      
      if(!is.null(models$x)) { # here if have covariates for x and y
        X.x <- model.matrix(as.formula(models$x),data=dat)
        nb.x <- dim(X.x)[2]
        nfixed <- nfixed-2
        
        if(length(b) != (nfixed+nb+nb.x)) 
          stop("length of b inconsistent with model and model.x")
        
        covb <- matrix(c(rep(b[1:nfixed],rep(n,nfixed)),X%*%b[nfixed+(1:nb)],X.x%*%b[nfixed+nb+(1:nb.x)]),
                       ncol=(nfixed+2))  
        covb <- as.vector(t(covb))
      } else { # here if have covariates for y only
        nfixed <- nfixed-1
        if(length(b) != (nfixed+nb)) 
          stop("length of b inconsistent with model and model.x")
        
        covb <- matrix(c(rep(b[1:(nfixed-1)],rep(n,(nfixed-1))),X%*%b[(nfixed-1)+(1:nb)],
                         rep(b[length(b)],n)),ncol=(nfixed+1))  #### bitchange 
        covb <- as.vector(t(covb))
      }
    } else { # here if have covariates for x only
      X.x <- model.matrix(as.formula(models$x),data=dat)
      n.x <- dim(X.x)[1]
      nb.x <- dim(X.x)[2]
      nfixed <- nfixed-1
      
      if(length(b) != (nfixed+nb.x)) 
        stop("length of b inconsistent with model and model.x")
      
      covb <- matrix(c(rep(b[1:nfixed],rep(n.x,nfixed)),X.x%*%b[nfixed+(1:nb.x)]),ncol=(nfixed+1))  
      covb <- as.vector(t(covb))        
    }
  } else {
    stop("Invalid FUN.")
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
#' 
#' @export
poiss.equiv <- function(availhmm,zero=0){
  Pois.availhmm <- availhmm
  Pi <- availhmm$Pi
  pcu <- FALSE
  
  if(is.element("pcu",names(availhmm))) {
    names(availhmm)[which(names(availhmm)=="pcu")] <- "pm"
    pcu <- TRUE
  }
  
  pm <- availhmm$pm
  delta <- availhmm$delta
  if(is.vector(pm)&!is.matrix(Pi) | !is.vector(pm)&is.matrix(Pi)) 
    stop("Single animal: pcu/pm is not a vector or Pi is not a matrix")
  
  if(is.vector(pm)) { # convert to matrix and array so can use loop below
    Pi <- array(Pi,dim=c(2,2,1))
    pm <- matrix(pm,ncol=1)
    delta <- matrix(delta,ncol=1)
  }
  
  nw <- dim(pm)[2]
  PiPoiss <- matrix(c(zero,1-zero,zero,1-zero),byrow=TRUE,nrow=2)
  
  for(i in 1:nw) {
    Pi[,,i] <- PiPoiss
    delta[,i] <- compdelta(PiPoiss)
  }
  
  E <- Estate <- matrix(rep(0,2*nw),nrow=2)
  
  if(nw>1) {
    for(w in 1:nw) {
      Estate[,w] <- c(1/availhmm$Pi[1,2,w],1/availhmm$Pi[2,1,w])
      events <- apply(Estate*pm,2,sum)
      duration <- apply(Estate,2,sum)
      eventrate <- rep(0,length(duration))
      eventrate[duration>0] <- events[duration>0]/duration[duration>0]
      Pois.availhmm$Pi <- Pi
      
      if(pcu) 
        Pois.availhmm$pcu <- matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      else 
        Pois.availhmm$pm <- matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      
      Pois.availhmm$delta <- delta
    }
  } else {
    Estate <- c(1/availhmm$Pi[1,2],1/availhmm$Pi[2,1])
    events <- sum(Estate*pm)
    duration <- sum(Estate)
    eventrate <- events/duration
    Pois.availhmm$Pi <- Pi[,,1]
    
    if(pcu) 
      Pois.availhmm$pcu <- c(0,eventrate)
    else 
      Pois.availhmm$pm <- c(0,eventrate)
    
    Pois.availhmm$delta <- as.vector(delta)
  }
  
  return(Pois.availhmm)
}

#' @title Logit function.
#'
#' @description
#'  Logit function
#'  
#' @param p probability (scalar or vector).
logit <- function(p) return(log(p/(1-p)))  # returns logit of p

#' @title Inverse logit function.
#'
#' @description
#'  Inverse logit function
#'  
#' @param x scalar or .
inv.logit <- function(x) return(1/(1+exp(-x))) # returns p from x=(logit of p)


