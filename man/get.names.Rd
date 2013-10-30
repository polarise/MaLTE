\name{get.names}
\alias{get.names}

\title{
Get the names of the training/test samples to be included in the output data.
}

\description{
The sample names will be required in the output data when carrying out downstream analyes. This function takes the samples file ('samples.txt') and extracts the required names.
}

\usage{
get.names(samples.fn = "samples.txt", test = TRUE)
}

\arguments{
  \item{samples.fn}{the name of the file containing a map of sample names between microarray and HTS data}
  \item{test}{if test=TRUE then get names of test samples; otherwise, return the names of training samples}
}

\value{
Returns a vector of either training or test sample names.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{get.train}}, \code{\link{get.test}}
}

\examples{
\dontrun{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
sample.names <- get.names( samples.fn )
}
}

\keyword{names}
