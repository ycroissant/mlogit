% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpar.R
\name{rpar}
\alias{rpar}
\alias{print.rpar}
\alias{summary.rpar}
\title{random parameter objects}
\usage{
rpar(x, par = NULL, norm = NULL, ...)

\method{print}{rpar}(
  x,
  digits = max(3, getOption("digits") - 2),
  width = getOption("width"),
  ...
)

\method{summary}{rpar}(object, ...)
}
\arguments{
\item{x, object}{a \code{mlogit} or a \code{rpar} object,}

\item{par}{the name or the index of the parameters to be extracted
; if \code{NULL}, all the parameters are selected,}

\item{norm}{the coefficient used for normalization if any,}

\item{...}{further arguments.}

\item{digits}{the number of digits}

\item{width}{the width of the printed output}
}
\value{
a \code{rpar} object, which contains:
\itemize{
\item dist: the name of the distribution,
\item mean: the first parameter of the distribution,
\item sigma: the second parameter of the distribution,
\item name: the name of the parameter.
}
}
\description{
\code{rpar} objects contain the relevant information about estimated
random parameters. The homonymous function extract on \code{rpar} object
from a \code{mlogit} object.
}
\details{
\code{mlogit} objects contain an element called \code{rpar} which contain a
list of \code{rpar} objects, one for each estimated random
parameter. The \code{print} method prints the name of the distribution
and the parameter, the \code{summary} behave like the one for numeric
vectors.
}
\seealso{
\code{\link[=mlogit]{mlogit()}} for the estimation of a random parameters logit
model.
}
\author{
Yves Croissant
}
\keyword{regression}
