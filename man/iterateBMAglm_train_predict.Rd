\name{iterateBMAglm.train.predict}
\alias{iterateBMAglm.train.predict}
\title{Iterative Bayesian Model Averaging: training and prediction}
\description{Classification and variable selection on microarray data.
	     This is a multivariate technique to select a small number
	     of relevant variables (typically genes) to classify
	     microarray samples.  This function performs the training,
	     and prediction steps.  The data is assumed to consist of
	     two classes. Logistic regression is used for classification.}
\usage{iterateBMAglm.train.predict (train.expr.set, test.expr.set, train.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1)}

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

\details{This function consists of the training phase and the prediction 
         phase.  The training phase consists of first
	 ordering all the variables (genes) by a univariate measure
	 called between-groups to within-groups sums-of-squares (BSS/WSS)
	 ratio, and then iteratively applying the \code{bic.glm} algorithm
	 from the \code{BMA} package.  The prediction phase uses the variables
	 (genes) selected in the training phase to classify the samples
	 in the test set. }

\value{A vector consisting of the predicted probability that each test 
       sample belongs to class 1 is returned.}

\references{
Raftery, A.E. (1995). 
Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{The \code{BMA} and \code{Biobase} packages are required.}

\seealso{\code{\link{iterateBMAglm.train}},  
	 \code{\link{iterateBMAglm.train.predict.test}},
	 \code{\link{brier.score}}
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
\keyword{multivariate}
\keyword{classif}



