% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rjags_fcts.R
\name{t.plot}
\alias{t.plot}
\title{Diagnostics plots}
\usage{
\method{t}{plot}(out)
}
\arguments{
\item{out}{Output from clustNSD}
}
\description{
This function uses the traceplot to display a plot of iterations vs. sampled values for each variable in the chain with a separate plot per variable.
}
\examples{
data(Christian_rjags)
t.plot(Christian_rjags)
}
\keyword{t.plot}
\keyword{traceplot}

