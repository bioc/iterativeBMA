\name{iterativeBMA-package}
\alias{iterativeBMA-package}
\alias{iterativeBMA}
\docType{package}
\title{
The Iterative Bayesian Model Averaging (BMA) algorithm
}
\description{
The iterative Bayesian Model Averaging (BMA) algorithm is a variable 
selection and classification algorithm with an application of 
classifying 2-class microarray samples, as described in Yeung, 
Bumgarner and Raftery (Bioinformatics 2005, 21: 2394-2402).
}
\details{
\tabular{ll}{
Package: \tab iterativeBMA\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2005-12-30\cr
License: \tab GPL version 2 or higher\cr
}
The function \code{iterateBMAglm.train} selects relevant variables by
iteratively applying the \code{bic.glm} function from the \code{BMA} package.
The data is assumed to consist of two classes. 
The function \code{iterateBMAglm.train.predict} combines the training
and prediction phases, and returns the predicted posterior probabilities
that each test sample belongs to class 1.
The function \code{iterateBMAglm.train.predict.test} combines the training,
prediction and test phases, and returns a list consisting of the
numbers of selected genes and models using the training data, the number
of classification errors and the Brier Score on the test set.
}
\author{
Ka Yee Yeung, University of Washington, Seattle, WA, with
contributions from Adrian Raftery and Ian Painter

Maintainer: Ka Yee Yeung <kayee@u.washington.edu>
}
\references{
Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\keyword{multivariate}
\keyword{classif}

\seealso{\code{\link{iterateBMAglm.train.predict}},  
	 \code{\link{iterateBMAglm.train.predict.test}},
	 \code{\link{bma.predict}},
	 \code{\link{brier.score}}
}

\examples{
library (Biobase)
library (BMA)
library (iterativeBMA)
data(trainData)
data(trainClass)

## training phase: select relevant genes
ret.bic.glm <- iterateBMAglm.train (train.expr.set=trainData, trainClass, p=100)

## get the selected genes with probne0 > 0
ret.gene.names <- ret.bic.glm$namesx[ret.bic.glm$probne0 > 0]

data (testData)

## get the subset of test data with the genes from the last iteration of bic.glm
curr.test.dat <- t(exprs(testData)[ret.gene.names,])

## to compute the predicted probabilities for the test samples
y.pred.test <- apply (curr.test.dat, 1, bma.predict, postprobArr=ret.bic.glm$postprob, mleArr=ret.bic.glm$mle)

## compute the Brier Score if the class labels of the test samples are known
data (testClass)
brier.score (y.pred.test, testClass)
}
