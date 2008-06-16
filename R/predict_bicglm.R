## goal: compute the posterior probability of Pr[Y=1 | D]
## Pr[Y=1|D] = sum_k postprob_k * (e^sum / 1+ e^sum)
## where sum = sum_j beta_k^j * Xj 
## where beta_k^j is the mle est of variable j in model k
## arguments: postprobArr = postprob from bic.glm
##	     mleArr = mle from bic.glm
##	     newdataArr = vector of test data (1 observation)
## value = NA if newdataArr contains NA

maxExpValueKY <- 700
minExpValueKY <- -700

bma.predict <- function (newdataArr, postprobArr, mleArr) {
       nModel <- length (postprobArr)
       ## sum over each model k
       retprob <- 0
       if (any(is.na(newdataArr))) {
       	  retprob <- NA
       } else {
       	 for (k in 1: nModel) {
           ## sum over all variables, assume the first mle is the intercept
	   myexp <- mleArr[k,1] + sum (mleArr[k,-1] * newdataArr)
           if (myexp > 0) {
              myexp <- min (myexp, maxExpValueKY)
           } else {
              myexp <- max (myexp, minExpValueKY)
           } 
	   myexpval <- exp (myexp)
	   retprob <- retprob + (postprobArr[k] * (myexpval / (1+ myexpval)))
       	 }
       }
       retprob
}


