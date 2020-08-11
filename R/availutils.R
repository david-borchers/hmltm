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
