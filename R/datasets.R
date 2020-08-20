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