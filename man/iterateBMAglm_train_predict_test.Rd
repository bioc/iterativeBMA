\name{iterateBMAglm.train.predict.test}
\alias{iterateBMAglm.train.predict.test}
\title{Iterative Bayesian Model Averaging: training, prediction and testing}
\description{Classification and variable selection on microarray data.
	     This is a multivariate technique to select a small number
	     of relevant variables (typically genes) to classify
	     microarray samples.  This function performs the training,
	     prediction and testing steps.  The data is assumed to consist of
	     two classes, and the classes of the test data is assumed to be
	     known. Logistic regression is used for classification.}
\usage{iterateBMAglm.train.predict.test (train.expr.set, test.expr.set, train.class, test.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1)}

\arguments{
\item{train.expr.set}{an \code{ExpressionSet} object.
		 We assume the rows in the expression data represent variables (genes), 
		 while the columns  represent 
		 samples or experiments. This training data is used to
		 select relevant genes (variables) for classification.}
\item{test.expr.set}{an \code{ExpressionSet} object.
		 We assume the rows in the expression data represent variables (genes), 
		 while the columns  represent samples or experiments.
		  The variables selected using the
		training data is used to classify samples on this test data.}
\item{train.class}{class vector for the observations (samples or 
                   experiments) in the training data. 
		   Class numbers are assumed to start from 0,
		   and the length of this class vector should be equal
		   to the number of rows in train.dat.
		   Since we assume 2-class data, we expect the class vector
		   consists of zero's and one's.}
\item{test.class}{class vector for the observations (samples or 
                  experiments) in the test data. 
		  Class numbers are assumed to start from 0,
		  and the length of this class vector should be equal
		  to the number of rows in test.dat.
		  Since we assume 2-class data, we expect the class vector
		  consists of zero's and one's.}
\item{p}{a number indicating the maximum number of top univariate genes
	 used in the iterative BMA algorithm.  This number is assumed to be
	 less than the total number of genes in the training data.
	 A larger p usually requires longer computational time as more
	 iterations of the BMA algorithm are potentially applied.
	 The default is 100.}
\item{nbest}{a number specifying the number of models of each size 
             returned to \code{bic.glm} in the \code{BMA} package. 
	     The default is 10.}
\item{maxNvar}{a number indicating the maximum number of variables used in
	       each iteration of \code{bic.glm} from the \code{BMA} package.
	       The default is 30.}
\item{maxIter}{a number indicating the maximum of iterations of 
               \code{bic.glm}. The default is 20000.}
\item{thresProbne0}{a number specifying the threshold for the posterior
                    probability that each variable (gene) is non-zero (in
		    percent).  Variables (genes) with such posterior 
		    probability less than this threshold are dropped in
		    the iterative application of \code{bic.glm}.  The default
		    is 1 percent.}
}

\details{This function consists of the training phase, prediction phase,
	 and the testing phase.  The training phase consists of first
	 ordering all the variables (genes) by a univariate measure
	 called between-groups to within-groups sums-of-squares (BSS/WSS)
	 ratio, and then iteratively applying the \code{bic.glm} algorithm
	 from the \code{BMA} package.  The prediction phase uses the variables
	 (genes) selected in the training phase to classify the samples
	 in the test set.  The testing phase assumes that the class labels
	 of the samples in the test set are known, and computes the number of 
	 classification errors and the Brier Score.}

\value{A list consisting of 4 elements are returned:
\item{num.genes}{The number of relevant genes selected using the training
		 data.}
\item{num.model}{The number of models selected using the training data.}
\item{num.err}{The number of classification errors produced when the
	       the predicted class labels of the test samples are compared
	       to the known class labels.}
\item{brierScore}{The Brier Score computed using the predicted and known
		  class labels of the test samples.  The Brier Score
		  represents a probabilistic number of errors. A small
		  Brier Score implies high prediction accuracy.}
}

\references{
Raftery, A.E. (1995). 
Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{The \code{BMA} and \code{Biobase} packages are required.}

\seealso{\code{\link{iterateBMAglm.train}},  
	 \code{\link{iterateBMAglm.train.predict}}
}

\examples{
library (Biobase)
library (BMA)
library (iterativeBMA)
data(trainData)
data(trainClass)
data (testData)
data (testClass)

iterateBMAglm.train.predict.test (train.expr.set=trainData, test.expr.set=testData, trainClass, testClass, p=100)

}
\keyword{multivariate}
\keyword{classif}



