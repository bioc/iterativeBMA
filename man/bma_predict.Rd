\name{bma.predict}
\alias{bma.predict}
\title{Predicted Probabilities from Bayesian Model Averaging}
\description{This function computes the predicted posterior probability 
             that each test sample belongs to class 1.  It assumes
             2-class data, and requires the true class labels to be known.}
\usage{bma.predict (newdataArr, postprobArr, mleArr)}

\arguments{
\item{newdataArr}{a vector consisting of the data from a test sample.}
\item{postprobArr}{a vector consisting of the posterior probability
                   of each BMA selected model.}
\item{mleArr}{matrix with one row per model and one column per variable giving 
              the maximum likelihood estimate of each coefficient for each 
              BMA selected model.}
}

\details{Let Y be the response variable (class labels for samples in our
         case).  In Bayesian Model Averaging (BMA), the posterior 
	 probability of Y=1 given the training set is the weighted 
	 average of the posterior probability of Y=1 given the training set 
	 and model M multiplied by the posterior probability of model M 
	 given the training set, summing over a set of models M.}

\value{A real number between zero and one, representing the predicted
       posterior probability.}

\references{
Raftery, A.E. (1995). 
Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{}

\seealso{\code{\link{brier.score}},
         \code{\link{iterateBMAglm.train}}
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
\keyword{classif}



