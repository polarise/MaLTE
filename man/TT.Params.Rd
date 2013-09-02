\name{TT.Params}
\alias{TT.Params}

\title{
Constructor for objects of class TT.Params
}

\description{
TT.Params objects store the parameters needed for training and prediction. These are used by the run method for objects of class TT.Ready.
}

\usage{
TT.Params(mtry = 2, ntree = 1000, feature.select = TRUE, min.probes = 15, cor.thresh = 0, OOB = FALSE)
}

\arguments{
  \item{mtry}{the number of predictors selected at random to try and effect a split on a node}
  \item{ntree}{the number of trees in the forest}
  \item{feature.select}{should we perform feature selection? default: TRUE}
  \item{min.probes}{the minimum number of probes used to effect feature selection}
  \item{cor.thresh}{feature selection is performed using the correlation of predictors with the response; the minimum correlation at which to accept predictors; default: zero}
  \item{OOB}{logical that determines whether test data will be used or use training data to perform out-of-bag predictions on the training data}
}

\value{Returns an object of class TT.Params}

\author{
  Paul K. Korir
}

\examples{
# create a TT.Params object using default parameters
tt.params <- TT.Params()
}

\keyword{ ~params }
\keyword{ ~parameters }
