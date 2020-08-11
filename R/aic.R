
#' @title Summarise AIC and related for hmltm object(s).
#'
#' @description
#' Tabulates a summary of AIC and related statistcs for fitted models of class
#' \code{hmltm}.
#'
#' @param object An object of class "hmltm", as produced by \link{est.hmltm}.
#' @param ... Other object of class "hmltm".
#' @param sort Logical to sort output on \code{criterion} or not.
#' @param k Value to multiply sample size by for AIC calculation
#' @param dmax Models with \code{citerion} difference values greater than this
#' (scalar) will be given zero (AIC or AICc) weight
#' @param criterion The criterion to use when comparing models and reporting differences:
#' "AICc", or "AIC".
#'
#' @return Summary of models, number of parameters, log-liklihoods, AIC and AICc values,
#' difference in AIC or AICc values and associated model weights
#'
#' @details Hacked from Murray Efford's \code{AIC.secr} function
#'

#' @export
#'
AIC.hmltm = function(object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c("AICc", "AIC")) {
  allargs <- list(...)
  modelnames <- (c(as.character(match.call(expand.dots = FALSE)$object),
                   as.character(match.call(expand.dots = FALSE)$...)))
  allargs <- hmltmlist(object, allargs)
  # Patch to fix fact that I stupidly made constant model NULL rather than "~1""
  for(i in 1:length(allargs)) for(j in 1:2) {
    if(is.null(allargs[[i]]$hmltm.fit$models[[1]])) allargs[[i]]$hmltm.fit$models[[1]] = "~1"
    if(is.null(allargs[[i]]$hmltm.fit$models[[2]])) allargs[[i]]$hmltm.fit$models[[2]] = "~1"
  }
  names(allargs) <- modelnames
  AIC(allargs, sort = sort, k = k, dmax = dmax, criterion = criterion)
  #  AIC.hmltmlist(allargs, sort = sort, k = k, dmax = dmax, criterion = criterion)
}



# ==========================================================
# Undocumented and unexported functions used by AIC.hmltm:
# ==========================================================
hmltmlist = function (...)
{
  dots <- match.call(expand.dots = FALSE)$...
  allargs <- list(...)
  if (length(allargs) == 1 & inherits(allargs[[1]], "hmltmlist")) {
    return(allargs[[1]])
  }
  else {
    if (is.null(names(allargs))) {
      dots2 <- substitute(list(...))[-1]
      names(allargs) <- sapply(dots2, deparse)
    }
    allargs <- lapply(allargs, function(x) if (inherits(x,"hmltm"))
      list(x)
      else x)
    temp <- do.call(c, allargs)
    if (is.null(names(temp)))
      names(temp) <- paste("hmltm", 1:length(temp), sep = "")
    if (!all(sapply(temp, function(x) inherits(x, "hmltm"))))
      stop("objects must be of class 'hmltm' or 'hmltmlist'")
    class(temp) <- "hmltmlist"
    temp
  }
}



logLik.hmltm = function(object) object$hmltm.fit$fit$value



AIC.hmltmlist = function (object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c("AICc", "AIC"))
{
  if (k != 2)
    warning("k != 2 and AIC.hmltm output may be mis-labelled")
  if (length(list(...)) > 0)
    warning("... argument ignored in 'AIC.hmltmlist'")
  criterion <- match.arg(criterion)
  modelnames <- names(object)
  allargs <- object
  if (class(allargs) != "hmltmlist")
    stop("'object' must be 'hmltmlist'.")
  output <- data.frame(t(sapply(allargs, oneline.hmltm, k = k)),
                       stringsAsFactors = F)
  for (i in 3:6) output[, i] <- as.numeric(output[, i])
  output$delta <- output[, criterion] - min(output[, criterion])
  OK <- abs(output$delta) < abs(dmax)
  sumdelta <- sum(exp(-output$delta[OK]/2))
  output$wt <- ifelse(OK, round(exp(-output$delta/2)/sumdelta,
                                4), 0)
  row.names(output) <- modelnames
  if (sort)
    output <- output[order(output[, criterion]), ]
  names(output)[7] <- paste("d", criterion, sep = "")
  names(output)[8] <- paste(criterion, "wt", sep = "")
  if (nrow(output) == 1) {
    output[, 8] <- NULL
    output[, 7] <- NULL
  }
  output
}


oneline.hmltm = function (hmltm, k) {
  n = length(as.numeric(na.omit(hmltm$dat$x)))
  Npar <- length(hmltm$hmltm.fit$fit$par)
  AICval <- 2 * hmltm$hmltm.fit$fit$value + k * Npar
  AICcval <- ifelse((n - Npar - 1) > 0, 2 * (hmltm$hmltm.fit$fit$value +
                                               Npar) + 2 * Npar * (Npar + 1)/(n - Npar - 1), NA)
#  c(Dmodel = paste(as.character(hmltm$model[[1]])[c(2,1,3)],collapse=""),
#    pmodel = paste(as.character(hmltm$model[[2]])[c(2,1,3)],collapse=""),
  c(y.model = as.character(hmltm$hmltm.fit$models[[1]]),
    x.model = as.character(hmltm$hmltm.fit$models[[2]]),
    npar = Npar,
    logLik = signif(hmltm$hmltm.fit$fit$value,5), AIC = round(AICval, 3), AICc = round(AICcval,3))
}

