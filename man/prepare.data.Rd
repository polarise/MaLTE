\name{prepare.data}
\alias{prepare.data}

\title{
Prepare data required for MaLTE gene processing
}

\description{
Runs a Python script to combine the microarray and high-throughput-sequencing (HTS) data into training and test sets according to the sample identifiers provided. All analysis is performed in the current working directory.
}

\usage{
prepare.data(samples.fn = "samples.txt", ma.fn = "ma_data.txt", hts.fn = "hts_data.txt", g2p.fn = "gene_probesets.txt", raw = FALSE, PCs=FALSE)
}

\arguments{
  \item{samples.fn}{the name of a file containing a map of sample names used in the headers of the microarray and HTS data files}
  \item{ma.fn}{the name of the file containing microarray probe intensities}
  \item{hts.fn}{the name of the file containing HTS (RNA-Seq) gene expression estimates}
  \item{g2p.fn}{the map of genes to probe sets; a probe set typically consists of about four probes}
  \item{raw}{logical; if TRUE treat the microarray probe intensities as the raw output from Affymetrix Power Tools; this implies different column positions for the probe and probe set identifiers}
  \item{PCs}{logical; are principal components present in the \code{samples.txt} file?}
}

\value{
Void but creates two files: \code{train_data.txt.gz} and \code{test_data.txt.gz} containing formatted training and test data.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{prepare.txs.data}}
}

\note{
This script calls an underlying Python script. For this reason, MaLTE presently only runs in GNU/Linux and UNIX-like systems.
}

\examples{
\dontrun{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )
}
}

\keyword{prepare}
\keyword{data}
\keyword{gene}
\keyword{training}
\keyword{test}
