#' Auxiliary function to form matrices equal to  x everywhere except from the upper/lower k diagonal, which values are vec
#'
#'@param  x initial matriw
#'@param vec vector with specified new values for the upper/lower k diagonal elements of x
#'@param k index of the diagonal
#'
#'@return x
#'
#'@examples
#'K=3
#'u <- sqrt(1/(4-1/seq(1,(K-1))^2))
#'n = length(u)+1
#'trans = myDiag(matrix(0,n, n),u,1) + myDiag(matrix(0,n, n),u,-1)
#'

myDiag <- function(x, vec, k) {
  x[row(x) == col(x) - k] <- vec
  x
}
