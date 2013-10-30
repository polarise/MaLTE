\name{TT.Ready.Gene}
\alias{TT.Ready.Gene}

\title{
Constructor for objects of class TT.Ready.Gene.
}

\description{
Each instance of a TT.Ready.Gene object holds the training and test data for a single gene. Objects of this class are required used by the run method to predict gene expression.
}

\usage{
TT.Ready.Gene(m)
}

\arguments{
  \item{m}{
\code{m} is a list having the following slots: \code{gene.id}, \code{no.train}, \code{no.test}, \code{no.probes}, \code{hts.train}, \code{hts.test}, \code{probes.train}, \code{probes.test}
}
}

\details{
Here follows a detailed description on the slots in \code{m}:
  \enumerate{
  \item{\code{gene.id}}{a unique gene identifier}
  \item{\code{no.train}}{the number of training samples}
  \item{\code{no.test}}{the number of test samples}
  \item{\code{no.probes}}{the number of probes on the microarray that overlap this genes}
  \item{\code{hts.train}}{a vector of length no.train of HTS expression estimates}
  \item{\code{hts.test}}{a vector of length no.test of HTS expression estimates (optional)}
  \item{\code{probes.train}}{a vectorised matrix of length no.train*no.probes of microarray probe intensities; remind me to tell you whether each row is a probe or sample}
  \item{\code{probes.test}}{a vectorised matrix of length no.test*no.probes of microarray probe intensities}
 }
}

\value{
Returns an object of class TT.Ready.Gene having the slots described above.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Ready-class}}, \code{\link{TT.Ready.Tx}}
}

\examples{
## Not run:
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
raw_ma.fn <- paste( system.file( package="MaLTE" ), "data", "raw_ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )
## End(Not run)
}

\keyword{gene}
\keyword{ready}
\keyword{data}
