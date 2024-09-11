#############################
### KERNEL MATRIX STUDY ###
#############################

### Functions for studying kernel matrix(ces) properties


## Single matrix

#' Kernel matrix histogram
#'
#' `histK()` plots the histogram of a kernel matrix.
#'
#' @details Information about the von Neumann entropy can be found at '?vonNeumann()'.
#'
#' @param K Kernel matrix (class "matrix").
#' @param main Plot title.
#' @param vn If TRUE, the value of the von Neumann entropy is shown in the plot.
#' (Defaults: FALSE).
#' @param ...  further arguments and graphical parameters passed to `plot.histogram`.
#'
#' @return An object of class "histogram".
#' @export
#' @importFrom graphics hist mtext
#'
#' @examples
#' data <- matrix(rnorm(150),ncol=50,nrow=30)
#' K <- RBF(data,g=0.01)
#' histK(K)

histK <- function(K,main="Histogram of K",vn = FALSE,...) {
  ## Errors
  kprecondition_helper(K)

  hist(get_lower_tri(K),main=main,xlab="Value",...)
  if(vn)  mtext(paste("von Neumann entropy = ",round(vonNeumann(K),digits=3)))
}


#' Kernel matrix heatmap
#'
#' `heatK()` plots the heatmap of a kernel matrix.
#'
#' @param K Kernel matrix (class "matrix").
#' @param cos.norm If TRUE, the cosine normalization is applied to the kernel matrix
#' so its elements have a maximum value of 1. (Defaults: FALSE).
#' @param title Heatmap title (optional).
#' @param color A vector of length 2 containing two colors. The first color will be
#' used to represent the minimum value and the second the maximum value of the kernel matrix.
#' @param raster In large kernel matrices, raster = TRUE will draw quicker and
#' better-looking heatmaps. (Defaults=FALSE).
#'
#' @return A `ggplot2` heatmap.
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @examples
#' data <- matrix(rnorm(150),ncol=50,nrow=30)
#' K <- Linear(data)
#' heatK(K)

heatK <- function(K,cos.norm=FALSE,title=NULL,color=c("red","yellow"),raster=FALSE) {
  ## Errors
  kprecondition_helper(K)

  if(cos.norm) K <- cosNorm(K)
  melted_cormat <- reshape2::melt(get_lower_tri(K), na.rm = TRUE)
  Var1 <- Var2 <- value <- NULL
  colnames(melted_cormat) < c("Var1","Var2","value")
  q <- ggplot2::ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value))
  if(raster) {
    q <- q+ geom_raster()
  } else {
    q <- q+ geom_tile(color = "white")
  }
    q <- q + scale_fill_gradient(low = color[1], high =color[2], limit = c(min(0,min(melted_cormat$value)),max(melted_cormat$value))) +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_blank())+  coord_fixed() +ggtitle(title)
    return(q)
}


#' Von Neumann entropy
#'
#' `vonNeumann()` computes the von Neumann entropy of a kernel matrix.
#' Entropy values close to 0 indicate that all its elements are very similar,
#' which may result in underfitting when training a prediction model. Instead,
#' values close to 1 indicate a high variability which may produce overfitting.
#'
#' @references Belanche-Muñoz, L.A. and Wiejacha, M. (2023)
#' Analysis of Kernel Matrices via the von Neumann Entropy and Its Relation to RVM Performances.
#' Entropy, 25, 154. doi:10.3390/e25010154. \href{https://www.mdpi.com/1099-4300/25/1/154}{Link}
#'
#' @param K Kernel matrix (class "matrix").
#'
#' @return Von Neumann entropy (a single value).
#' @export
#'
#' @examples
#' data <- matrix(rnorm(150),ncol=50,nrow=30)
#' K <- Linear(data)
#' vonNeumann(K)

vonNeumann <- function(K) {
  ## Errors
  kprecondition_helper(K)

  autov <- eigen(K)$values
  autov <- autov[autov>0]
  autov <- autov/sum(autov)
  return(-sum(autov*log2(autov))/log2(nrow(K)))
}


#' Gamma hyperparameter estimation (RBF kernel)
#'
#' @description
#' This function returns an estimation of the optimum value for the gamma hyperparameter
#' (required by the RBF kernel function) using different heuristics:
#'
#' \describe{
#'   \item{\emph{D} criterion}{It returns the inverse of the number of features in X.}
#'   \item{Scale criterion}{It returns the inverse of the number of features,
#'   normalized by the total variance of X.}
#'   \item{Quantiles criterion}{A range of values, computed with the function
#'   `kernlab::sigest()`.}
#' }
#'
#' @param X Matrix or data.frame that contains real numbers ("integer", "float" or "double").
#'
#' @return A list with the gamma value estimation according to different criteria.
#' @importFrom kernlab sigest
#' @export
#'
#' @examples
#' data <- matrix(rnorm(150),ncol=50,nrow=30)
#' gamma <- estimate_gamma(data)
#' gamma
#' K <- RBF(data, g = gamma$scale_criterion)
#' K[1:5,1:5]

estimate_gamma <- function(X) {
  X <- as.matrix(X)
  qu <- kernlab::sigest(X,frac=1,scaled=FALSE)
  return(list(d_criterion=1/ncol(X),scale_criterion=1/(ncol(X)*stats::var(as.vector(X))),
              quantiles_criterion=qu))
}


## Two matrices, or matrix vs target

#' Kernel matrix similarity
#'
#' `simK()` computes the similarity between kernel matrices.
#'
#' @details It is a wrapper of `Frobenius()`.
#'
#' @param Klist A list of \emph{M} kernel matrices with identical \emph{NxN} dimension.
#' @return Kernel matrix (dimension: \emph{MxM}).
#'
#' @export
#'
#' @examples
#' K1 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))
#' K2 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))
#' K3 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))
#'
#' simK(list(K1,K2,K3))

simK <- function(Klist) {
  message("Remember that Klist should contain only kernel matrices (i.e. squared, symmetric and PSD).
  This function does NOT verify the symmetry and PSD criteria.")

  if(any(sapply(Klist,nrow) != sapply(Klist,ncol)))  stop("At least one matrix is not squared")
  return(Frobenius(DATA=Klist,cos.norm = TRUE))

  # Forma alternativa:
  # similarities <- outer(1:length(Klist), 1:length(Klist),
  #                       FUN = Vectorize(function(i, j) {
  #                         out <- tr(Klist[[i]] %*% Klist[[j]])
  #                        # out <- out / (norm(Klist[[i]], type="F") * norm(Klist[[j]], type="F"))
  #                         return(out)
  #                       }))
  # rownames(similarities) <- colnames(similarities) <- names(Klist)
  # return(similarities)
}


#' Kernel-target alignment
#'
#' `KTA()` computes the alignment between a kernel matrix and a target variable.
#'
#' @param K A kernel matrix (class: "matrix").
#' @param y The target variable. A numeric vector or a factor with two levels.
#' @return Alignment value.
#'
#' @export
#' @importFrom methods is
#'
#' @examples
#' K1 <- RBF(iris[1:100,1:4],g=0.1)
#' y <- factor(iris[1:100,5])
#' KTA(K1,y)

# # De moment només funciona amb SVM per classificació. y ha de ser binària.
# crec que per regressió hauria de servir també
# # Tenc la sensació que sempre dóna molt baix.

KTA <- function(K,y) {
  if(methods::is(y,"factor")) {
    if(nlevels(y)>2) stop("y should have 2 levels")
    y2 <- 1*(y==levels(y)[1])
    y2[y2==0] <- -1
  }
  ## Errors
  kprecondition_helper(K)

  return(simK(list(K,Linear(y2)))[1,2])
}


## Helpers

#' Helper for the heatmap of kernel matrices
#' @keywords internal
#' @noRd
get_lower_tri <- function(K){
  K[upper.tri(K)]<- NA
  return(K)
}
