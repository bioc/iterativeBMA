\name{imageplot.iterate.bma}
\alias{imageplot.iterate.bma}
\title{An image plot visualization tool}
\description{Create a visualization of the models and variables selected by the
             iterative BMA algorithm.}
\usage{imageplot.iterate.bma (bicreg.out, color="default", ...)}

\arguments{
\item{bicreg.out}{An object of type 'bicreg', 'bic.glm' or 'bic.surv'}
\item{color}{The color of the plot. The value "default" uses the 
             current default R color scheme for image. The value 
             "blackandwhite" produces a black and white image.}
\item{...}{Other parameters to be passed to the image and axis functions.}
}

\details{This function is a modification of the \code{imageplot.bma}
         function from the \code{BMA} package. The difference is that
	 variables (genes) with \code{probne0} equal to 0 are removed
	 before plotting. The arguments of this function is identical
	 to those in \code{imageplot.bma}.}

\value{An heatmap-style image, with the BMA selected variables on the vertical
       axis, and the BMA selected models on the horizontal axis. The variables
       (genes) are sorted in descreasing order of the posterior probability 
       that the variable is not equal to 0 (\code{probne0}) from top to 
       bottom.  The models are sorted in descreasing order of the
       model posterior probability (\code{postprob}) from left to right.
       }

\references{
Clyde, M. (1999) 
Bayesian Model Averaging and Model Search Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.

Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005) 
Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data. 
Bioinformatics 21: 2394-2402.
}
\note{The \code{BMA} and \code{Biobase} packages are required.}

\seealso{\code{\link{iterateBMAglm.train}}
}

\examples{
library (Biobase)
library (BMA)
library (iterativeBMA)
data(trainData)
data(trainClass)

## training phase: select relevant genes
ret.bic.glm <- iterateBMAglm.train (train.expr.set=trainData, trainClass, p=100)

## produce an image plot to visualize the selected genes and models
imageplot.iterate.bma (ret.bic.glm)

}
\keyword{classif}



