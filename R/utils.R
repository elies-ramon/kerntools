#######################
### General HELPERS ###
#######################


#' Trace of a matrix (should be squared)
#' @keywords internal
#' @noRd
tr <- function(X) sum(diag(X))


#' Matrix or data frame to factor
#' @keywords internal
#' @noRd
toFactor <- function(X) {
  if(!methods::is(X,"data.frame"))  {
    X <- as.data.frame(X,stringsAsFactors = T)
  } else {
    X[] <- lapply(X,as.factor)
  }
  return(X)
}


#' Kernel matrix class, test if squared and numeric
#' @keywords internal
#' @noRd
kprecondition_helper <- function(K) {
  if(nrow(K) != ncol(K))  stop("K should be squared (and symmetric and PSD)")
  if(!methods::is(K,"matrix")) stop("K should be class = matrix")
  if(!methods::is(diag(K),"numeric")) stop("K should contain only numbers")
}

