\name{read.data}
\alias{read.data}

\title{
Read and collate training and test data into lists of gene objects
}

\description{
Reads in the training and test data obtained from the raw dataa and converts them into appropriate format for subsequent analysis.
}

\usage{
read.data(train.fn, test.fn)
}

\arguments{
  \item{train.fn}{the name of a file containing training data having both HTS (RNA-Seq) and microarray data for supervised learning}
  \item{test.fn}{the name of the file containg the test data; it may contain HTS data for verification}
}

\value{
Returns a list of objects of class \code{TT.Ready.Gene}.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{read.txs.data}}
}

\examples{
## Not run:
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )
## End(Not run)
}

\keyword{gene}
\keyword{data}
