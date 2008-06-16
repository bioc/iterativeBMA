\name{BssWssFast}
\alias{BssWssFast}
\title{Between-groups sum-of-squares to within-groups sum-of-squares
       ratio}
\description{This is a univariate technique to select relevant
	     genes in classification of microarray data. In classifying
	     samples of microarray data, this ratio is computed for each
	     gene.  A large between-groups to within-groups sum-of-squares 
	     ratio indicates a potentially relevant gene.}
\usage{BssWssFast (X, givenClassArr, numClass = 2)}

\arguments{
\item{X}{data matrix where columns are variables and rows are observations.
         In the case of gene expression data, the columns (variables) represent
	 genes, while the rows (observations) represent samples or experiments.}
\item{givenClassArr}{class vector for the observations (samples or 
                     experiments). Class numbers are assumed to start from 0,
		     and the length of this class vector should be equal
		     to the number of rows in X.
		     In the case of 2-class data, we expect the class vector
		     consists of zero's and one's.}
\item{numClass}{number of classes. The default is 2.}
}

\details{This function is called by \code{iterateBMAglm.2class}.}

\value{A list of 2 elements are returned:
\item{x}{A vector containing the BSS/WSS ratios in descending order.}
\item{ix}{A vector containing the indices corresponding to the sorted ratios.}
}

\references{
Dudoit, S., Fridlyand, J. and Speed, T.P. (2002) Comparison of discrimination methods for the classification of tumors using gene expression data. Journal of the American Statistical Association 97: 77-87.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{}

\seealso{\code{\link{iterateBMAglm.train}},  \code{\link{trainData}},
	 \code{\link{trainClass}}}

\examples{
data(trainData)
data(trainClass)

ret.bsswss <- BssWssFast (X=t(exprs(trainData)), givenClassArr=trainClass, numClass = 2)

}
\keyword{univar}


