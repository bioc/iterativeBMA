
## BssWssFast: compute the ratio of between-group to within group sum of squares
## arguments: X is a matrix where columns are variables, rows are observations
## return values: a vector of sorted (in descending order) ratio's and their ranks
## to get the top N variables: X[, sortedvec$ix[1:N]]

BssWssFast <- function (X, givenClassArr, numClass = 2) {
       
       ## create indicator variables for each class
       ## assume class numbers start from 0
       classVec <- matrix (0, numClass, length (givenClassArr))
       for (k in 1: numClass) {
           temp <- rep(0, length(givenClassArr))
       	   temp [givenClassArr == (k-1)] <- 1
	   classVec[k, ] <- temp
       }

       ## compute BSS/WSS for each variable
       classMeanArr <- rep(0, numClass)
       ratio <- rep (0, ncol(X))
       for (j in 1: ncol(X)) {
       	   ## compute mean across all samples
	   overallMean <- sum (X[, j])/length(X[, j])
	   ## compute class mean for variable j
	   for (k in 1: numClass) {
	       classMeanArr[k] <- sum(classVec[k, ] * X[, j])/sum(classVec[k, ])
       	   }
	   ## to re-write classMeanArr as a vector containing class mean for each entry
	   classMeanVec <- classMeanArr[givenClassArr+1]

	   ## compute BSS and WSS
	   bss <- sum( (classMeanVec - overallMean)^2)
	   wss <- sum( (X[, j] - classMeanVec)^2)
	   ratio[j] <- bss/wss
       } 
       sort (ratio, decreasing=TRUE, index = TRUE)
}
