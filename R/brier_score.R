## goal: to summarize pi^ (predicted P[Y=1]) using different scores

## for 2 classes only
## Brier score = sum_i (y_i - p_i^)^2

brier.score <- function (predictedArr, truthArr) {
	## make sure the 2 arrays are the same length
	if (length(predictedArr) != length(truthArr)) {
	   print ("ERROR: length NOT equal!!")
        }
 	## sum over all the samples
	temp.vec <- (truthArr - predictedArr)^2
	sum (temp.vec)
}





