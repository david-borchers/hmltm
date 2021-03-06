% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sintegral}
\alias{sintegral}
\title{Integration using Simpson's method.}
\usage{
sintegral(fx, x, n.pts = 16, type = "int")
}
\arguments{
\item{fx}{Values at \code{x} (see below) of the function to be integrated.}

\item{x}{Values at which function fx has been evaluated.}

\item{n.pts}{The number of points to be used in integration. If \code{x} contains more than
\code{n.pts} then \code{n.pts} will be set to \code{length(x)}.}

\item{type}{if equal to `\code{cdf}' a list comprising x-values and cdf values at each x-value 
is returned, else the integral over the range of \code{x} is returned.}
}
\description{
Modified version of function sintegral from library Bolstad.
}
