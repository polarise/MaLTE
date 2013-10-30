\name{tx.id}
\alias{tx.id}

\title{
Get the transcript identifier(s) of objects of class \code{TT.Seq.Tx}.
}

\description{
Accessor function for transcript identifier for \code{TT.Seq.Tx} objects.
}

\usage{
tx.id(object)
}

\arguments{
  \item{object}{an object of class \code{TT.Seq.Tx}}
}

\value{
Returns a vector of transcript identifiers.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{gene.id}}
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

hts.txs.fn <- paste( system.file( package="MaLTE" ), "data", "hts_txs_data.txt.gz", sep="/" )
g2tx.fn <- paste( system.file( package="MaLTE" ), "data", "gene_transcripts.txt.gz", sep="/" )

# prepare training and test data for transcripts
prepare.txs.data( samples.fn=samples.fn, train.fn="train_data.txt.gz", test.fn="test_data.txt.gz", hts.txs.fn=hts.txs.fn, g2tx.fn=g2tx.fn )

# read in the training and test data
tt.ready.txs <- read.txs.data( train.fn="train_txs_data.txt.gz", test.fn="test_txs_data.txt.gz" )

# set params
tt.params <- TT.Params( )

# train and predict
tt.seq.txs <- array2seq( tt.ready.txs, tt.params )

id <- tx.id( tt.seq.txs[[1]] )
ids <- lapply( tt.seq.txs, tx.id )
## End(Not run)
}

\keyword{transcript}
\keyword{identifier}
\keyword{ID}
