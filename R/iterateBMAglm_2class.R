
## iterative bic.glm in R: 
zero.threshold <- 0.0001

iterateBMAglm <- function (x, y, curr.mat, stopVar=0, nextVar, nbest=10, thresProbne0 = 1, maxNvar = 30, maxIter=20000) {

   ## iterative bic.glm
   currIter <- 0
   while (stopVar == 0 && currIter < maxIter) {
   	## run bic.glm from the BMA package
   	ret.bic.glm <- bic.glm (x=curr.mat, y=y, maxCol=(maxNvar+1), glm.family="binomial", nbest=nbest)

   	## get a logical vector for which probne0 < thresProbne0
   	rmVector <- which (ret.bic.glm$probne0 < thresProbne0)

	## now, check if we could add more variables
	## We want to run bic.glm with as many variables as possible
	## so that we are considering more models simultaneously AND
	## for efficiency reasons
	if (length(ret.bic.glm$probne0) < maxNvar) {
	  ## we have extra space, so we want to add variables
	  add.slots <- maxNvar - length(ret.bic.glm$probne0)
	  ## stuff up curr.mat with junk, to be removed later
	  num.var.to.add <- ncol(x) - nextVar + 1
	  if (num.var.to.add > length(rmVector)) {
	     ## we only add if we have more variables left to add than
	     ## space from removing low threshold ones
	     num.spots.needed <- num.var.to.add - length(rmVector)
	     junk.ncol <- min (num.spots.needed, add.slots)
	     junk.mat <- matrix(0, nrow=nrow(curr.mat), ncol=junk.ncol)
	     curr.mat <- cbind (curr.mat, junk.mat)
	     rmVector.old <- rmVector
	     rmVector <- c(rmVector.old, (length(ret.bic.glm$probne0)+1):(length(ret.bic.glm$probne0)+junk.ncol))
	  }
	}

	## only increase thresProbne0 if there is no space to add in variables
	if (length(rmVector) == 0) {
	  ## we only end up here if length(ret.bic.glm$probne0) >= maxNvar
	  ## no gene to swap in!!, increase threshold
	  currMin <- min (ret.bic.glm$probne0)
	  print (paste("no gene to swap! Min probne0 = ", currMin, sep=""))
          newThresProbne0 <- currMin + 1
	  print (paste("new probne0 threshold = ", newThresProbne0, sep=""))
	  rmVector <- which (ret.bic.glm$probne0 < newThresProbne0)
	  ## print (rmVector)
	}

  	## now, guaranteed there is at least 1 variable to swap
	## or there is space to add in new variables

	if (nextVar <= ncol(x)) {
	   ## set up new X
	   ## print ("set up new X")
	   lastVar <- length(rmVector) + nextVar - 1
	   ## to make sure lastVar <= ncol(x)
	   if (lastVar > ncol(x)) {
	     rmVector <- rmVector [1: (ncol(x) -nextVar + 1)]
	     lastVar <- ncol(x)
           }
	   ## print (paste(nextVar, lastVar, sep=" "))
	   ##print (dimnames(x)[[2]][nextVar:lastVar])
	   
	   ## check for singularity
	   curr.mat[, rmVector] <- x[, nextVar:lastVar]
	   ## change the colume name as well!!
	   dimnames(curr.mat)[[2]][rmVector] <- dimnames(x)[[2]][nextVar:lastVar]
	   svd.out <- svd (curr.mat)

           if (any(svd.out$d < zero.threshold) == TRUE) {
		## singularity, try to add one variable at a time
		print ("inside svd")
		curr.mat <- curr.mat[, -(rmVector)]
		for (curr.var in nextVar:lastVar) {
		   orig.dimnames <- dimnames (curr.mat)[[2]]
		   curr.mat <- cbind (curr.mat, x[, curr.var])
		   dimnames(curr.mat)[[2]] <- c(orig.dimnames, dimnames(x)[[2]] [curr.var])
		   svd.out <- svd (curr.mat)
		   if (any(svd.out$d < zero.threshold) == TRUE) {
			## NOT okay to add this variable, rm the last column
			curr.mat <- curr.mat[, -(ncol(curr.mat))]
		   }
		}
      	   }

	   nextVar <- lastVar + 1

       	} else {
	   ## there is no variable to be removed OR exhausted all data
           stopVar <- 1
        }
	currIter <- currIter + 1
   }
   print (paste(currIter, ": explored up to variable ## ", nextVar-1, sep=""))
  
   ## print out selected genes if iterateBMA is over
   if (stopVar == 1) {
	selectedGenes <- which(ret.bic.glm$probne0 >= thresProbne0)
   }
   
   list(curr.mat=curr.mat, stopVar=stopVar, nextVar=nextVar)
}


iterateBMAinit <- function (x, maxNvar = 30) {

   maxNvar <- min (maxNvar, ncol(x))
   curr.mat <- x[, 1:maxNvar]
   stopVar <- 0
   nextVar <- maxNvar + 1

   ## make sure the curr.mat is not singular
   svd.out <- svd (curr.mat)

   ## successively add a column
   ## if a variable causes singularity --> a linear combination of others
   ##           --> just ignore this variable
   if (any(svd.out$d < zero.threshold) == TRUE) {
	## start by adding 1 variable at a time
	curr.mat <- as.matrix(x[, 1])
	dimnames (curr.mat)[[2]] <- as.vector(dimnames(x)[[2]][1], mode="list") 
	curr.var <- 2
	while ((ncol(curr.mat) < maxNvar) && (curr.var <= ncol(x))) {
	  orig.dimnames <- dimnames (curr.mat)[[2]]
	  curr.mat <- cbind (curr.mat, x[, curr.var])
	  dimnames(curr.mat)[[2]] <- c(orig.dimnames, dimnames(x)[[2]] [curr.var])
	  print (paste("iterateBMAinit: add 1 column --> total ", ncol(curr.mat), sep=""))
	  svd.out <- svd (curr.mat)
	  if (any(svd.out$d < zero.threshold) == TRUE) {
		## NOT okay to add this variable, rm the last column
	  	orig.dimnames.2 <- dimnames (curr.mat)[[2]]
		orig.dimnames.1 <- dimnames (curr.mat)[[1]]
		col.to.rm <- ncol(curr.mat)
		curr.mat <- as.matrix (curr.mat[, -(col.to.rm)])
		dimnames(curr.mat) <- list(orig.dimnames.1, orig.dimnames.2[-col.to.rm])
		print (paste ("variables dropped: ", ncol(curr.mat)+1, sep=""))
	  }
	  curr.var <- curr.var + 1
	}
	## here, ncol(curr.mat) >= maxNvar OR curr.var > ncol(x)
	nextVar <- curr.var
   }

   list(curr.mat=curr.mat, stopVar=stopVar, nextVar=nextVar)
}



## parameters: nbest: input to bicglm
##             sortedA: matrix of independent variables, sorted by columns
##                  expected G * C, all of our data
##             y: vector of response
##             thresProbne0: threshold to rm probne0, default = 1%
##	      maxIter: max ## iterations in repeating bic.glm
##             maxNvar: max ## variables to feed into bic.glm, default=30
## returns: an object of class bic.glm if the algorithm finishes
## 	   output from applying bic.glm on selected variables
##	   otherwise, returns -1

iterateBMAglm.wrapper <- function (sortedA, y, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1) {
  ## get the top "maxNvar" variables
  ret.bma.init <- iterateBMAinit (x=sortedA, maxNvar)

  ## call bic.glm repeatedly
  ret.bma <- iterateBMAglm (x=sortedA, y=y, curr.mat=ret.bma.init$curr.mat, stopVar=ret.bma.init$stopVar, nextVar=ret.bma.init$nextVar, nbest, thresProbne0, maxNvar, maxIter)

  if (ret.bma$stopVar == 1) {
     ## apply bic.glm again using selected genes
     ret.bic.glm <- bic.glm (x=ret.bma$curr.mat, y=y, maxCol=(maxNvar+1), glm.family="binomial", nbest=nbest)
     return (ret.bic.glm)
  } else {
     return (-1)
  }
}


########################################################


iterateBMAglm.train <- function (train.expr.set, train.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1) {

  train.dat <- t(exprs(train.expr.set))
  ## 1. order all the columns (variables) by BSS/WSS ratios
  sorted.vec <- BssWssFast (train.dat, train.class, numClass = 2)
  sorted.train.dat <- train.dat[, sorted.vec$ix[1:p]]
 
  ## 2. run iterative bic.glm
  ## run iterative BMA using logistic regression
  ret.bma <- iterateBMAglm.wrapper (sorted.train.dat, y=train.class, nbest=nbest, maxNvar=maxNvar, maxIter=maxIter, thresProbne0=thresProbne0)

  # convert gene names in "namesx" and "label"
  ret.bma$namesx <- as.character(sapply (ret.bma$namesx, convertSingleName, orig.expr.set=train.expr.set))
  ret.bma$label <- as.character(sapply (ret.bma$label, convertModelName, orig.expr.set=train.expr.set))
  ret.bma
}


iterateBMAglm.train.predict <- function (train.expr.set, test.expr.set, train.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1) {

  train.dat <- t(exprs(train.expr.set))
  test.dat <- t(exprs(test.expr.set))
  ## use the training data to select relevant genes
  ret.bma <- iterateBMAglm.train (train.expr.set, train.class, p, nbest, maxNvar, maxIter, thresProbne0)

  ## predict the classes of the test samples using the selected genes
  if (nrow(test.dat) > 1) {
     curr.test.dat <- test.dat [, ret.bma$namesx]
  } else {    # for LOOCV
     curr.test.dat <- matrix(test.dat [, ret.bma$namesx], ncol=length(ret.bma$namesx), dimnames=list(NULL, ret.bma$namesx))
  }
  y.pred.test <- apply (curr.test.dat, 1, bma.predict, postprobArr=ret.bma$postprob, mleArr=ret.bma$mle)
}


iterateBMAglm.train.predict.test <- function (train.expr.set, test.expr.set, train.class, test.class, p=100, nbest=10, maxNvar=30, maxIter=20000, thresProbne0=1) {

  train.dat <- t(exprs(train.expr.set))
  test.dat <- t(exprs(test.expr.set))
  ## use the training data to select relevant genes
  ret.bma <- iterateBMAglm.train (train.expr.set, train.class, p, nbest, maxNvar, maxIter, thresProbne0)

  ## predict the classes of the test samples using the selected genes
  if (nrow(test.dat) > 1) {
     curr.test.dat <- test.dat [, ret.bma$namesx]
  } else {    # for LOOCV
     curr.test.dat <- matrix(test.dat [, ret.bma$namesx], ncol=length(ret.bma$namesx), dimnames=list(NULL, ret.bma$namesx))
  }
  y.pred.test <- apply (curr.test.dat, 1, bma.predict, postprobArr=ret.bma$postprob, mleArr=ret.bma$mle)

  ## compute the number of errors and brier score
  curr.table <- table (round(y.pred.test, 0), test.class)
  curr.error <- sum (curr.table) - sum(diag(curr.table))
  ## print (curr.table)
  ## compute the brier score
  ## print (round(y.pred.test, 3))
  curr.score <- brier.score (y.pred.test, test.class)
  ## print (curr.score)

  selected.genes <- which (ret.bma$probne0 > 0)
  ## print (paste("# selected genes = ", length(selected.genes), sep=""))
  ## print (paste("# selected models = ", length(ret.bma$postprob), sep=""))
  list (num.genes=length(selected.genes), num.model=length(ret.bma$postprob), num.err=curr.error, brierScore=curr.score)

}

















