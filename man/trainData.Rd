\name{trainData}
\docType{data}
\alias{trainData}
\title{Sample Training Data for the Iterative BMA Algorithm}
\description{
This is an adapted leukemia (ALL, AML) \code{ExpressionSet} from Golub et al. (1999).
This \code{ExpressionSet} consists of the expression levels from 38 ALL or AML samples
(rows), and 100  genes (columns). This dataset is used
as an example training data in our examples.
}
\usage{data(trainData)}
\format{The \code{ExpressionSet} is called \code{trainData}. Each entry
in the exprs matrix represents the expression level of one gene from an ALL
or AML sample.
}
\details{
For illustration purposes, a subset of  100 genes from the package \code{golubEsets}
 is included in this package.
} 
\source{
The \code{golubEsets} bioconductor data package, or
\url{http://www.broad.mit.edu/cgi-bin/cancer/datasets.cgi}.
}
\references{
Golub, T.R., Slonim, D.K., Tamayo, P., Huard, C., Gaasenbeek, M., Mesirov, J.P., Coller, H., Loh, M.L., Downing, J.R., Caligiuri, M.A., et al. (1999) 
Molecular classification of cancer: class discovery and class prediction by gene expression monitoring. Science 286: 531-7.
}
\keyword{datasets}

