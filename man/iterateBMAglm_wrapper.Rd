\name{iterateBMAglm.wrapper}
\alias{iterateBMAglm.wrapper}
\title{Iterative Bayesian Model Averaging}
\description{This function repeatedly calls \code{bic.glm} from the
	     \code{BMA} package until all variables are exhausted.
	     The data is assumed to consist of
	     two classes. Logistic regression is used for classification.}
\usage{iterateBMAglm.wrapper (sortedA, y, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1) }

\arguments{
\item{sortedA}{data matrix where columns are variables and rows are 
	       observations.  The variables (columns) are assumed to
	       be sorted using a univariate measure.
               In the case of gene expression data, the columns (variables) 
	       represent genes, while the rows (observations) represent 
	       samples or experiments.}
\item{y}{class vector for the observations (samples or 
         experiments) in the training data. 
	 Class numbers are assumed to start from 0,
	 and the length of this class vector should be equal
	 to the number of rows in sortedA.
	 Since we assume 2-class data, we expect the class vector
	 consists of zero's and one's.}
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

\details{In this function, the variables are assumed to be sorted, and
	 \code{bic.glm} is called repeatedly.  In the first application of
	 the \code{bic.glm} algorithm, the top \code{maxNvar} univariate
	 ranked genes are used.  After each application of the \code{bic.glm}
	 algorithm, the genes with \code{probne0} < \code{thresProbne0}
	 are dropped, and the next univariate ordered genes are added
	 to the BMA window.
	 The function \code{iterateBMAglm.train} calls \code{BssWssFast} before
	 calling this function.
	 Using this function, users can experiment with alternative 
	 univariate measures.}

\value{If all variables are exhausted, an object of class 
       \code{bic.glm} returned by the last iteration
       of \code{bic.glm}. Otherwise, -1 is returned. 
       The object of class \code{bic.glm} is a list consisting 
       of the following components:
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

\seealso{\code{\link{iterateBMAglm.train}},  
	 \code{\link{iterateBMAglm.train.predict}},
	 \code{\link{iterateBMAglm.train.predict.test}},
	 \code{\link{BssWssFast}}
}

\examples{
library (Biobase)
library (BMA)
library (iterativeBMA)
data(trainData)
data(trainClass)

## Use the BSS/WSS ratio to rank all genes in the training data
sorted.vec <- BssWssFast (t(exprs(trainData)), trainClass, numClass = 2)
## get the top ranked 50 genes
sorted.train.dat <- t(exprs(trainData[sorted.vec$ix[1:50], ]))
 
## run iterative bic.glm
ret.bic.glm <- iterateBMAglm.wrapper (sorted.train.dat, y=trainClass)

## The above commands are equivalent to the following 
ret.bic.glm <- iterateBMAglm.train (train.expr.set=trainData, trainClass, p=50)

}
\keyword{multivariate}
\keyword{classif}



