\name{iterativeBMA-internal}
\alias{iterativeBMA-internal}
\alias{iterateBMAglm}
\alias{iterateBMAinit}
\alias{maxExpValueKY}
\alias{minExpValueKY}
\alias{zero.threshold}
\alias{imageplot.bma.mod}
\alias{convertSingleName}
\alias{convertModelName}
\alias{bma.punct.string}

\title{Internal functions for iterativeBMA}
\description{Internal functions for \code{iterativeBMA}, not meant to be called directly.}

\usage{
iterateBMAglm (x, y, curr.mat, stopVar=0, nextVar, nbest=10, thresProbne0 = 1, maxNvar = 30, maxIter=20000)
iterateBMAinit (x, maxNvar = 30)
imageplot.bma.mod (bicreg.out, color = "default", ...)
convertSingleName (curr.name, orig.expr.set)
convertModelName (curr.model, orig.expr.set)
maxExpValueKY
minExpValueKY
zero.threshold
bma.punct.string
}
\keyword{internal}

