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
