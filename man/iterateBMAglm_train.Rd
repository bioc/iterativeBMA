\name{iterateBMAglm.train}
\alias{iterateBMAglm.train}
\title{Iterative Bayesian Model Averaging: training step}
\description{Classification and variable selection on microarray data.
	     This is a multivariate technique to select a small number
	     of relevant variables (typically genes) to classify
	     microarray samples.  This function performs the training phase.
	     The data is assumed to consist of
	     two classes. Logistic regression is used for classification.}
\usage{iterateBMAglm.train (train.expr.set, train.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1)}

\arguments{
\item{train.expr.set}{an \code{ExpressionSet} object.
		 We assume the rows in the expression data represent variables (genes), 
		 while the columns  represent 
		 samples or experiments. This training data is used to
		 select relevant genes (variables) for classification.}
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

\details{The training phase consists of first
	 ordering all the variables (genes) by a univariate measure
	 called between-groups to within-groups sums-of-squares (BSS/WSS)
	 ratio, and then iteratively applying the \code{bic.glm} algorithm
	 from the \code{BMA} package.  In the first application of
	 the \code{bic.glm} algorithm, the top \code{maxNvar} univariate
	 ranked genes are used.  After each application of the \code{bic.glm}
	 algorithm, the genes with \code{probne0} < \code{thresProbne0}
	 are dropped, and the next univariate ordered genes are added
	 to the BMA window.}

\value{An object of class \code{bic.glm} returned by the last iteration
       of \code{bic.glm}.  The object is a list consisting of
       the following components:
\item{namesx}{the names of the variables in the last iteration of 
              \code{bic.glm}.}
\item{postprob}{the posterior probabilities of the models selected.}
\item{deviance}{the estimated model deviances.}
\item{label}{labels identifying the models selected.}
\item{bic}{values of BIC for the models.}
\item{size}{the number of independent variables in each of the models.}
\item{which}{a logical matrix with one row per model and one column per 
             variable indicating whether that variable is in the model.}
\item{probne0}{the posterior probability that each variable is non-zero 
               (in percent).}
\item{postmean}{the posterior mean of each coefficient (from model averaging).}
\item{postsd}{the posterior standard deviation of each coefficient 
              (from model averaging).}
\item{condpostmean}{the posterior mean of each coefficient conditional on 
                    the variable being included in the model.}
\item{condpostsd}{the posterior standard deviation of each coefficient 
                  conditional on the variable being included in the model.}
\item{mle}{matrix with one row per model and one column per variable giving 
           the maximum likelihood estimate of each coefficient for each model.}
\item{se}{matrix with one row per model and one column per variable giving 
          the standard error of each coefficient for each model.}
\item{reduced}{a logical indicating whether any variables were dropped 
               before model averaging.}
\item{dropped}{a vector containing the names of those variables dropped 
               before model averaging.}
\item{call}{the matched call that created the bma.lm object.}
}

\references{
Raftery, A.E. (1995). 
Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{The \code{BMA} and \code{Biobase} packages are required.}

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

## show the posterior probabilities of selected models
ret.bic.glm$postprob

data (testData)

## get the subset of test data with the genes from the last iteration of bic.glm
curr.test.dat <- t(exprs(testData)[ret.gene.names,])

## to compute the predicted probabilities for the test samples
y.pred.test <- apply (curr.test.dat, 1, bma.predict, postprobArr=ret.bic.glm$postprob, mleArr=ret.bic.glm$mle)

## compute the Brier Score if the class labels of the test samples are known
data (testClass)
brier.score (y.pred.test, testClass)
}
\keyword{multivariate}
\keyword{classif}



