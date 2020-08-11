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
