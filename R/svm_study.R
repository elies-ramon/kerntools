###################
### SVM HELPERS ###
###################

### Extracting and studying alphas and SVM importance from ksvm() (kernlab)


#' SVM feature importance
#'
#' Recovering the features importances from a SVM model.
#' @details This function may be not valid for all kernels. Do not use it with
#' the RBF, Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels unless
#' you know exactly what you are doing.
#'
#' Usually the sign of the importances is irrelevant, thus justifying working with the
#' absolute or squared values; see for instance Guyon et al. (2002). Some classification
#' tasks are an exception to this, when it can be demonstrated that the feature space
#' is strictly nonnegative. In that case, a positive importance implies that a feature
#'  contributes to the "positive" class, and the same with a negative importance
#'  and the "negative" class.
#'
#' @references Guyon, I., Weston, J., Barnhill, S., and Vapnik, V. (2002) Gene selection
#' for cancer classification using support vector machines. Machine learning, 46, 389-422.
#' \href{https://link.springer.com/content/pdf/10.1023/a:1012487302797.pdf}{Link}
#'
#' @param X Matrix or data.frame that contains real numbers ("integer", "float" or "double").
#' X is NOT the kernel matrix, but the original dataset used to compute the kernel matrix.
#' @param svindx Indices of the support vectors.
#' @param coeff target * alpha.
#' @param result A string. If "absolute", the absolute values of the importances
#' are returned. If "squared", the squared values are returned. Any other input will
#' result in the original (positive and/or negative) importance values (see Details). (Defaults: "absolute").
#' @param cos.norm  Boolean. Was the data cosine normalized prior to training the model? (Defaults: FALSE).
#' @param center Boolean. Was the data centered prior to training the model? (Defaults: FALSE).
#' @param scale Boolean. Was the data scaled prior to training the model? (Defaults: FALSE).
#'
#' @return The importance of each feature (a vector).
#'
#' @export
#'
#' @examples
#' data1 <- iris[1:100,]
#' sv_index <- c( 24, 42, 58, 99)
#' coefficients <- c(-0.2670988, -0.3582848,  0.2129282,  0.4124554)
#' # This SV and coefficients were obtained from a model generated with kernlab:
#' # model <- kernlab::ksvm(Species ~ .,data=data1, kernel="vanilladot",scaled = TRUE)
#' # sv_index <- unlist(kernlab::alphaindex(model))
#' # coefficients <- kernlab::unlist(coef(model))
#' # Now we compute the importances:
#' svm_imp(X=data1[,-5],svindx=sv_index,coeff=coefficients,center=TRUE,scale=TRUE)

svm_imp <- function(X,svindx,coeff,result="absolute",cos.norm=FALSE, center=FALSE,scale=FALSE) {
  message("Do not use this function if the SVM model was created with the RBF,
          Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels")
  # y <- as.numeric(target[aindx])
  # y[y==1] <- -1
  # y[y==2] <- 1
  if(scale|center) X <- scale(X, center=center, scale=scale)
  if(cos.norm) X <- cosnormX(X)

  # return(colSums(coeff*(dat[aindx,])))
  feats <- crossprod(coeff, as.matrix(X[svindx,]))[,,drop=TRUE]
  if(result=="absolute") {
    return(abs(feats))
  } else if(result=="squared") {
    return(feats^2)
  } else {
    return(feats)
  }
}


