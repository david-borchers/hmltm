\name{hsltm-package}
\alias{hsltm-package}
\alias{hsltm}
\docType{package}
\title{
Hidden State Line Transect Models
}
\description{
Estimates density and abundance from line transect survey data using forward distances to detections and estimates of availability process parameters.
}
\details{
\tabular{ll}{
Package: \tab hsltm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-05\cr
License: \tab GNU General Public License Version 2 or later\cr
}
Package \code{hsltm} (for "Hidden State Line Transect Model") estimates density and abundance from single-platform line transect surveys on which animals are not continuously available for detection. Unlike conventional line transect methods, hsltm methods do not assume that all animals, or even all available animals, on the transect line are detected. They assume only that all available animals at radial distance zero are detected.\cr

In addition to the perpendicular distances to all detections required by conventional line transect methods, hsltm methods require (1) forward distances to detections and (2) estimates of the parameters of Markov models (MMs), hidden Markov models (HMMs) or Markov modulated Poisson process models (MMPPs - not yet implemented) describing animal availability patterns. More than one set of availability model paramters can be specified, in which case the estimator treats the availability model parameters as random effects.
}
\author{
David Borchers and Chris Jefferson

Maintainer: David Borchers <dlb@st-andrews.ac.uk>
}
\references{
Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics. DOI: 10.1111/biom.12049

Langrock, R., Borchers, D.L. and Skaug, H. Markov-modulated nonhomogeneous Poisson processes for unbiased estimation of marine mammal abundance. 2013. Journal of the American Statistical Association. DOI: 10.1080/01621459.2013.797356
}
\note{
  The datasets bowhead.adat and bowhead.depths are not to be used for research
   or publication without prior consent from Mads Peter Heide-Jorgensen, email:
  mhj<``at'' symbol>ghsdk.dk.
}
\keyword{Line transect, availability bias, Hidden markov model}
\seealso{
Some packages for fitting HMMs to availability data: \code{\link{HiddenMarkov}},  \code{msm}
}
\examples{
library(hsltm)

# get data and availability model parameters
data(aerial.survey) # get bowhead aerial survey line transect data
data(bowhead.hmm.pars) # get bowhead availability HMM parameters (8 sets of parameters)

# set survey and fitting parameters
survey.pars=make.survey.pars(spd=43.6,Wl=100,W=3000,ymax=2200) # specify survey parameters
control.fit=list(hessian=FALSE,nx=64) # fitting parameters
control.opt=list(trace=5,maxit=1000) # optimisation parameters

# specify model and starting parameter values and fit model
hfun="h.EP2x.0";models=list(y=NULL,x=NULL) # specify detection hazard model (no covariates)
pars=c(1.66, 0.63, 64, 877) # detection hazard parameter starting values
# Point estimation:
W.est=2500 # set perpendicular truncation distance for Horvitz-Thompson-like estimator
bwEP2x.null=est.hmltm(aerial.survey,pars,hfun,models,survey.pars,bowhead.hmm.pars,
                      control.fit,control.opt,W.est=W.est)

# display point estimates of density and abundance:
bwEP2x.null$point$ests
}