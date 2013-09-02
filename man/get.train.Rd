\name{get.train}
\alias{get.train}

\title{
Get the names of the training samples to be included in the output data.
}

\description{
The sample names will be required in the output data when carrying out downstream analyes. This function takes the samples file ('samples.txt') and extracts the required names.
}

\usage{
get.train(samples.fn = "samples.txt", test = FALSE)
}

\arguments{
  \item{samples.fn}{the name of the file containing a map of sample names between microarray and HTS data}
  \item{test}{set to FALSE to get names of training samples}
}

\value{
Returns the names of training samples.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{get.names}}, \code{\link{get.test}}
}

\examples{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
sample.names <- get.train( samples.fn )
}

\keyword{names}
\keyword{training}
