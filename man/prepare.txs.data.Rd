\name{prepare.txs.data}
\alias{prepare.txs.data}

\title{
Prepare data required for MaLTE transcript isoform processing
}

\description{
Runs a Python script to combine the microarray and high-throughput-sequencing (HTS) data into training and test sets according to the sample identifiers provided. All analysis is performed in the current working directory.
}

\usage{
prepare.txs.data(samples.fn = "samples.txt", train.fn = "train_data.txt", test.fn = "test_data.txt", hts.txs.fn = "hts_txs_data.txt", g2tx.fn = "gene_transcripts.txt")
}

\arguments{
  \item{samples.fn}{the name of a file containing a map of sample names used in the headers of the microarray and HTS data files}
  \item{train.fn}{the name of the file containing gene training data; the output from running \code{prepare.data} then \code{read.data}}
  \item{test.fn}{the name of the file containing gene test data}
  \item{hts.txs.fn}{the name of the file containing HTS (RNA-Seq) transcript isoform expression estimates}
  \item{g2tx.fn}{a map of genes to transcript identifiers; may be downloaded from BioMart}
}

\value{
Void but creates two files: \code{train_txs_data.txt.gz} and \code{test_txs_data.txt.gz} containing formatted training and test data.
}

\author{
Paul K. Korir
}

\note{
\code{prepare.data} must be run before running \code{prepare.txs.data} because the output of the former is required by the latter. This script calls an underlying Python script. For this reason, MaLTE presently only runs in GNU/Linux and UNIX-like systems.
}

\seealso{
\code{\link{prepare.data}}
}

\examples{
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
}

\keyword{prepare}
\keyword{data}
\keyword{transcript}
\keyword{training}
\keyword{test}
