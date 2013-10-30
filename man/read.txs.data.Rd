\name{read.txs.data}
\alias{read.txs.data}

\title{
Read and collate training and test data into lists of transcript objects
}

\description{
Reads in the training and test data obtained from the raw dataa and converts them into appropriate format for subsequent analysis.
}

\usage{
read.txs.data(train.fn, test.fn)
}

\arguments{
  \item{train.fn}{the name of a file containing training data having both HTS (RNA-Seq) and microarray data for supervised learning}
  \item{test.fn}{the name of the file containg the test data; it may contain HTS data for verification}
}

\value{
Returns a list of objects of class \code{TT.Ready.Tx}.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{read.data}}
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
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

# locate of other input files
hts.txs.fn <- paste( system.file( package="MaLTE" ), "data", "hts_txs_data.txt.gz", sep="/" )
g2tx.fn <- paste( system.file( package="MaLTE" ), "data", "gene_transcripts.txt.gz", sep="/" )

# prepare training and test data for transcripts isoforms
prepare.txs.data( samples.fn=samples.fn, train.fn="train_data.txt.gz", test.fn="test_data.txt.gz", hts.txs.fn=hts.txs.fn, g2tx.fn=g2tx.fn )
tt.ready.txs <- read.txs.data( train.fn="train_txs_data.txt.gz", test.fn="test_txs_data.txt.gz" )
## End(Not run)
}

\keyword{transcript}
\keyword{data}
