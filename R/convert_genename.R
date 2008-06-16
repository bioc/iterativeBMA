# define the punctunation character string
bma.punct.string <- "-|/|\\.|#|&|\\*|'|~|\\{|\\}|\\(|\\)|:|;|\\+|=|%|$"

# to convert gene names from "namesx" of an bic.glm object to original gene names
convertSingleName <- function (curr.name, orig.expr.set) {
  curr.string <- unlist(strsplit (curr.name, bma.punct.string))
  ret.string <- curr.name

  # only need to convert if the string contains the period
  # map the R processed gene name back to what is in the original ExpressionSet
  if (length (curr.string) > 1) {
    match.vec <- sapply (dimnames(exprs(orig.expr.set))[[1]], function (x) {
      temp.string <- unlist(strsplit (x, bma.punct.string))
      if (length(temp.string) == length(curr.string)) {
        ret.val <- all (temp.string == curr.string)
      } else {
        ret.val <- FALSE
      }
      ret.val})
    temp.ind <- which (match.vec == T)
    if (length(temp.ind) == 1) {
      # exactly 1 match
      ret.string <- names(match.vec[temp.ind])
    } 
  }
  ret.string
}


# to convert model namess from "label" of an bic.glm object to original gene names
convertModelName <- function (curr.model, orig.expr.set) {
  model.arr <- unlist (strsplit (curr.model, ","))
  ret.names <- sapply (model.arr, convertSingleName, orig.expr.set=orig.expr.set)
  paste (ret.names, collapse=",")
}

