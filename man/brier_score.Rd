\name{brier.score}
\alias{brier.score}
\title{Brier Score: assessment of prediction accuracy}
\description{The Brier Score is a probabilistic number of errors that takes
	     the predicted probabilities into consideration. A small Brier
	     Score indicates high prediction accuracy. This function assumes
             2-class data, and requires the true class labels to be known.}
\usage{brier.score (predictedArr, truthArr)}

\arguments{
\item{predictedArr}{a vector consisting of the predicted probabilities
                    that the test sample belongs to class 1.}
\item{truthArr}{a zero-one vector indicating the known class labels of
                the test samples. We assume this vector has the same
		length as \code{predictedArr}.}
}

\details{The Brier Score computes the sum of squares of the differences 
         between the true class and the predicted probability over all 
         test samples. If the predicted probabilities are constrained 
         to equal to 0 or 1, the Brier Score is equal to the total 
         number of classification errors.}

\value{A non-negative real number.}

\references{
Brier, G.W. (1950) 
Verification of forecasts expressed in terms of probability. Monthly Weather Review 78: 1-3.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{}

\seealso{\code{\link{bma.predict}},
         \code{\link{iterateBMAglm.train.predict}}
}

\examples{
library (Biobase)
library (BMA)
library (iterativeBMA)
data(trainData)
data(trainClass)
data (testData)

ret.vec <- iterateBMAglm.train.predict (train.expr.set=trainData, test.expr.set=testData, trainClass, p=100)

## compute the Brier Score
data (testClass)
brier.score (ret.vec, testClass)

}
\keyword{classif}



