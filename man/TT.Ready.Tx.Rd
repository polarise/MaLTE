\name{TT.Ready.Tx}
\alias{TT.Ready.Tx}

\title{
Constructor for objects of class TT.Ready.Tx
}

\description{
Each instance of a TT.Ready.Tx object holds the training and test data for the set of transcript isoforms of a single gene. Objects of this class are used by the run method to predict transcript isoform expression.
}

\usage{
TT.Ready.Tx(m)
}

\arguments{
  \item{m}{
\code{m} is a list with the following slots: \code{gene.id}, \code{tx.id}, \code{no.txs}, \code{no.train}, \code{no.test}, \code{no.probes}, \code{hts.train}, \code{hts.test}, \code{probes.train}, and  \code{probes.test}.
}
}

\details{
Here follows a detailed description of the slots:
  \enumerate{
    \item{gene.id}{a unique gene identifier}
    \item{tx.id}{a vector of transcript identifiers}
    \item{no.txs}{the number of transcript isoforms in this gene for which we wish to predict expression}
    \item{no.train}{the number of training samples}
    \item{no.test}{the number of test samples}
    \item{no.probes}{the number of probes on the microarray that overlap this genes}
    \item{hts.train}{a vectorised matrix of length no.train*no.txs of HTS expression estimates}
    \item{hts.test}{a vectorised matrix of length no.test*no.txs of HTS expression estimates (optional)}
    \item{probes.train}{a vectorised matrix of length no.train*no.probes of microarray probe intensities; remind me to tell you whether each row is a probe or sample}
    \item{probes.test}{a vectorised matrix of length no.test*no.probes of microarray probe intensities}
  }
}

\value{
Returns an object of class TT.Ready.Tx with the slots above.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Ready-class}}, \code{\link{TT.Ready.Gene}}
}

\examples{
\dontrun{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

hts.txs.fn <- paste( system.file( package="MaLTE" ), "data", "hts_txs_data.txt.gz", sep="/" )
g2tx.fn <- paste( system.file( package="MaLTE" ), "data", "gene_transcripts.txt.gz", sep="/" )

prepare.txs.data( samples.fn=samples.fn, train.fn="train_data.txt.gz", test.fn="test_data.txt.gz", hts.txs.fn=hts.txs.fn, g2tx.fn=g2tx.fn )
tt.ready.txs <- read.txs.data( train.fn="train_txs_data.txt.gz", test.fn="test_txs_data.txt.gz" )
}
}

\keyword{transcript}
\keyword{data}
