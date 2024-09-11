#################################
### Data fusion using kernels ###
#################################

### Combining data via kernel matrices


## Combining kernel matrices

#' Multiple Kernel (Matrices) Combination
#'
#' Combination of kernel matrices coming from different datasets / feature types
#' into a single kernel matrix.
#'
#' @param K A three-dimensional \emph{NxDxM} array containing \emph{M} kernel matrices.
#' @param coeff A vector of length \emph{M} with the weight of each kernel matrix.
#' If NULL, all kernel matrices have the same weight. (Defaults: NULL)
#'
#' @return A kernel matrix.
#'
#' @examples
#'
#' # For illustrating a possible use of this function, we work with a dataset
#' # that contains numeric and categorical features.
#'
#' summary(mtcars)
#' cat_feat_idx <- which(colnames(mtcars) %in% c("vs", "am"))
#'
#' # vs and am are categorical variables. We make a list, with the numeric features
#' # in the first element and the categorical features in the second:
#' DATA <- list(num=mtcars[,-cat_feat_idx], cat=mtcars[,cat_feat_idx])
#' # Our N, D and M dimensions are:
#' N <- nrow(mtcars); D <- ncol(mtcars); M <- length(DATA)
#'
#' # Now we prepare a kernel matrix:
#' K <- array(dim=c(N,N,M))
#' K[,,1] <- Linear(DATA[[1]],cos.norm = TRUE) ## Kernel for numeric data
#' K[,,2] <- Dirac(DATA[[2]]) ## Kernel for categorical data
#'
#' # Here, K1 has the same weight than K2 when computing the final kernel, although
#' # K1 has 9 variables and K2 has only 2.
#' Kconsensus <- MKC(K)
#' Kconsensus[1:5,1:5]
#'
#' # If we want to weight equally each one of the 11 variables in the final
#' # kernel, K1 will weight 9/11 and K2 2/11.
#' coeff <- sapply(DATA,ncol)
#' coeff
#' Kweighted <- MKC(K,coeff=coeff)
#' Kweighted[1:5,1:5]
#'
#' @importFrom methods is hasArg
#' @export

MKC <- function(K,coeff=NULL) {
  K <- aperm(K,c(3,1,2))
  if(!methods::hasArg(coeff)) return(colMeans(K))
  if(dim(K)[1] != length(coeff)) stop("Length of the coefficients vector different to the number of matrices")
  if(sum(coeff) != 1) coeff <- coeff / sum(coeff)
  return(colSums(coeff*K))
}

