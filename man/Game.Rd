% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Game}
\alias{Game}
\alias{Game2}
\title{Ranked data for gaming platforms}
\format{
A dataframe containing :
\itemize{
\item ch.Platform: where \code{platform} is one of \code{Xbox},
\code{PlayStation}, \code{PSPortable}, \code{GameCube},
\code{GameBoy} and \code{PC}. This variables contain the ranking of
the platforms from 1 to 6,
\item own.Platform: these 6 variables are dummies which indicate
whether the given plaform is already owned by the respondent,
\item age: the age of the respondent,
\item hours: hours per week spent on gaming.,
}
}
\source{
\href{https://journaldata.zbw.eu/journals/jae}{Journal of Applied Econometrics data archive}.
}
\description{
A sample of 91 Dutch individuals
}
\details{
The data are also provided in long format (use in this case
\code{data(Game2)}. In this case, the alternative and the choice
situation are respectively indicated in the \code{platform} and
\code{chid} variables.
}
\references{
\insertRef{FOK:PAAP:VAND:12}{mlogit}
}
\keyword{datasets}
