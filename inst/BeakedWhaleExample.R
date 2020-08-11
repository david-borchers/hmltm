library(hmltm)

# get data and availability model parameters
data(beaked.ship) # get beaked whale ship survey line transect data
data(beaked.hmm.pars) # get beaked availability HMM parameters (one set of parameters)

# set survey and fitting parameters
spd=1852*(6.577778)/(60*60)        # vessel speed in metres/second
#W.est=3500 # set perpendicular truncation distance for Horvitz-Thompson-like estimator
#survey.pars=make.survey.pars(spd=spd,Wl=100,W=W.est,ymax=5000) # specify survey parameters
W.est=2000 # set perpendicular truncation distance for Horvitz-Thompson-like estimator
survey.pars=make.survey.pars(spd=spd,Wl=0,W=W.est,ymax=4500) # specify survey parameters
control.fit=list(hessian=FALSE,nx=64) # fitting parameters
control.opt=list(trace=5,maxit=1000) # optimisation parameters

# specify model and starting parameter values and fit model
hfun="h.EP2x.0";models=list(y=NULL,x=NULL) # specify detection hazard model (no covariates)
pars=c(1.66, 0.63, 64, 877) # detection hazard parameter starting values
# Point estimation:
bkEP2x.null=est.hmltm(beaked.ship,pars,hfun,models,survey.pars,beaked.hmm.pars,
                      control.fit,control.opt,W.est=W.est)
# display point estimates of density and abundance:
bkEP2x.null$point$ests
# display probability of detection at zero distance
bkEP2x.null$hmltm.fit$pzero

# specify model and starting parameter values and fit model
hfun="h.IP.0";models=list(y=NULL,x=NULL) # specify detection hazard model (no covariates)
pars=c(2.482797, 31.280377); models=list(y=NULL,x=NULL)
bkIP.null=est.hmltm(beaked.ship,pars,hfun,models,survey.pars,beaked.hmm.pars,
                    control.fit,control.opt,W.est=W.est)
# display point estimates of density and abundance:
bkIP.null$point$ests

# display probability of detection at zero distance
bkIP.null$hmltm.fit$pzero


# Fit with platorm height
pars=c(3.124305, 75.493550, 1.092958); models=list(y="~ht",x=NULL)
bkIP.ht=est.hmltm(beaked.ship,pars,hfun,models,survey.pars,beaked.hmm.pars,
                    control.fit,control.opt,W.est=W.est)
# display point estimates of density and abundance:
bkIP.ht$point$ests

# display probability of detection at zero distance
bkIP.ht$hmltm.fit$pzero


# Fit with platorm height and group size
pars=c(3.124305, 75.493550, 1.092958, 1); models=list(y="~ht+size",x=NULL)
bkIP.ht.size=est.hmltm(beaked.ship,pars,hfun,models,survey.pars,beaked.hmm.pars,
                  control.fit,control.opt,W.est=W.est)
# display point estimates of density and abundance:
bkIP.ht.size$point$ests

# display probability of detection at zero distance
bkIP.ht.size$hmltm.fit$pzero

# Look at AICs of the above models
AIC.hmltm(bkIP.null,bkIP.ht,bkIP.ht.size,criterion="AIC")

# Do some plots of the best model, and get goodness-of-fit stats
ybreaks=c(0,175,500,1000,2000,3000,4200)
xbreaks=c(0,100,300,500,1000,1500,2000)

EP1bw.bfx=fxfit.plot(bkIP.ht.size,breaks=xbreaks,type="prob",allx=FALSE)
pl.EP1bw=fyfit.plot(bkIP.ht.size,breaks=ybreaks,allx=FALSE,nys=250)
EP1bw.gof.x=hmmlt.gof.x(bkIP.ht.size);EP1bw.gof.x$p.ks
EP1bw.gof.y=hmmlt.gof.y(bkIP.ht.size);EP1bw.gof.y$p.ks

