\name{get.test}
\alias{get.test}

\title{
Get the names of the test samples to be included in the output data.
}

\description{
The sample names will be required in the output data when carrying out downstream analyes. This function takes the samples file ('samples.txt') and extracts the required names.
}

\usage{
get.test(samples.fn = "samples.txt", test = TRUE)
}

\arguments{
  \item{samples.fn}{the name of the file containing a map of sample names between microarray and HTS data}
  \item{test}{set to TRUE to get names of test samples}
}

\value{
Returns a vector of test sample names.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{get.names}}, \code{\link{get.train}}
}

\examples{
## Not run:
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
sample.names <- get.test( samples.fn )
## End(Not run)
}

\keyword{names}
\keyword{test}
