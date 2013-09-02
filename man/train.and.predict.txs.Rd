\name{train.and.predict.txs}
\alias{train.and.predict.txs}

\title{
Function that performs training and prediction for each transcript isoform.
}

\description{
Training is currently implemented using the conditional random forest algorithm but is open to improvements in machine learning algorithms. For all transcript isoforms within each gene, this function learns a model based on HTS expression estimates and microarray probe intensities then either predicts the HTS expression or an out-of-bag estimate for a completely new set of probes for the transcript isoforms of the same gene. The function \code{run} is a member function of TT.Ready.Tx class and calles \code{train.and.predict.txs}.
}

\usage{
train.and.predict.txs(object, params.object)
}

\arguments{
  \item{object}{an object of class \code{TT.Seq.Tx}}
  \item{params.object}{an object of class \code{TT.Params}}
}

\value{
Returns an object of class \code{TT.Seq.Tx}.
}

\author{
Paul K. Korir
}

\note{
This function is used by \code{array2seq} expression prediction function.
}

\seealso{
\code{\link{train.and.predict}}, \code{\link{array2seq}}
}

\examples{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

hts.txs.fn <- paste( system.file( package="MaLTE" ), "data", "hts_txs_data.txt.gz", sep="/" )
g2tx.fn <- paste( system.file( package="MaLTE" ), "data", "gene_transcripts.txt.gz", sep="/" )

# prepare training and test transcript data
prepare.txs.data( samples.fn=samples.fn, train.fn="train_data.txt.gz", test.fn="test_data.txt.gz", hts.txs.fn=hts.txs.fn, g2tx.fn=g2tx.fn )

# read in the training and test data
tt.ready.txs <- read.txs.data( train.fn="train_txs_data.txt.gz", test.fn="test_txs_data.txt.gz" )

# default parameters
params.object <- TT.Params()

# train and predict
a.txs <- train.and.predict.txs( tt.ready.txs[[1]], params.object )
a.txs
predictions( a.txs )
sample.names <- get.names( samples.fn )
}

\keyword{train}
\keyword{predict}
\keyword{transcript}
