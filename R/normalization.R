############################################
### NORMALIZATION & KERNEL NORMALIZATION ###
############################################

### Several widespread normalization techniques


## Kernel normalization

#' Cosine normalization of a kernel matrix
#'
#' It is equivalent to compute K using the normalization `X/sqrt(sum(X^2))` in Feature Space.
#'
#' @references Ah-Pine, J. (2010). Normalized kernels as similarity indices.
#' In Advances in Knowledge Discovery and Data Mining: 14th Pacific-Asia Conference,
#' PAKDD 2010, Hyderabad, India, June 21-24, 2010. Proceedings. Part II 14 (pp. 362-373).
#' Springer Berlin Heidelberg. \href{https://hal.science/hal-01504523/document}{Link}
#'
#' @param K Kernel matrix (class "matrix").
#'
#' @return Cosine-normalized K (class "matrix").
#' @export
#' @importFrom methods is
#'
#' @examples
#' dat <- matrix(rnorm(250),ncol=50,nrow=5)
#' K <- Linear(dat)
#' cosNorm(K)

cosNorm <- function(K) {
  ## Errors
  kprecondition_helper(K)

  names <- rownames(K)
  D <- diag(1/sqrt(diag(K)))
  K <- D %*% K %*% D
  colnames(K) <- rownames(K) <- names
  return(K)
}


#' Centering a kernel matrix
#'
#' It is equivalent to compute `K` over centered data (i.e. the mean of each
#' column is subtracted) in Feature Space.
#'
#' @param K Kernel matrix (class "matrix").
#'
#' @return Centered `K` (class "matrix").
#' @export
#' @importFrom methods is
#' @examples
#' dat <- matrix(rnorm(250),ncol=50,nrow=5)
#' K <- Linear(dat)
#' centerK(K)

centerK <- function(K) {
  ## Errors
  kprecondition_helper(K)

  names <- rownames(K)
  n <- nrow(K)
  Ip <- matrix(rep(1/n,n^2),nrow=n,ncol=n)
  K <- (K - (Ip %*% K) - (K %*% Ip) +  (Ip %*% K %*% Ip))
  colnames(K) <- rownames(K) <- names
  return(K)
}



## General matrix and data normalization

#' Frobenius normalization
#'
#' This function computes the Frobenius normalization of a matrix.
#'
#' @param X Numeric matrix of any size. It may be a kernel matrix.
#'
#' @return Frobenius-normalized X (class: "matrix").
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(50),ncol=5,nrow=10)
#' frobNorm(dat)

# Podria ser útil en MKL. També es pot usar sobre matrius no de kernel i no quadrades.
frobNorm <- function(X)  return(X/(norm(X, type="F")))


#'
#' This function deletes those columns and/or rows in a matrix/data.frame that
#' only contain 0s.
#'
#' @param X Numeric matrix or data.frame of any size.
#' @param dim A numeric vector. 1 indicates that the function should be applied
#' to rows, 2 to columns, c(1, 2) indicates rows and columns. (Defaults: 2).

#' @return X with less rows or columns. (Class: the same than X).
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(150),ncol=50,nrow=30)
#' dat[c(2,6,12),] <- 0
#' dat[,c(30,40,50)] <- 0
#' dim(desparsify(dat))
#' dim(desparsify(dat,dim=c(1,2)))

desparsify <- function(X,dim=2)  {
  X0 <- X!=0
  if(1 %in% dim) X <- X[rowSums(X0)>0,,drop=FALSE]
  if(2 %in% dim) X <- X[,colSums(X0)>0,drop=FALSE]
  return(X)
}


#' Total Sum Scaling
#'
#' This function transforms a dataset from absolute to relative frequencies
#' (by row or column).
#'
#' @param X Numeric matrix or data.frame of any size containing absolute frequencies.
#' @inheritParams centerX
#' @return A relative frequency matrix or data.frame with the same dimension than X.
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(50),ncol=5,nrow=10)
#' TSS(dat) #It can be checked that, after scaling, the sum of each row is equal to 1.

TSS <- function(X,rows=TRUE) {
  if(rows) {
    return(X/rowSums(X))
  } else {
    return(X/colSums(X))
  }
}


#' Cosine normalization of a matrix
#'
#' Normalizes a numeric matrix dividing each row (if rows=TRUE) or column (if rows=FALSE)
#' by their L2 norm. Thus, each row (or column) has unit norm.
#'
#' @param X Numeric matrix or data.frame of any size.
#' @inheritParams centerX
#' @return Cosine-normalized X.
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(50),ncol=5,nrow=10)
#' cosnormX(dat)

cosnormX <- function(X,rows=TRUE) {
  if(rows) {
    return(X/sqrt(rowSums(X^2)))
  } else {
    return(X/sqrt(colSums(X^2)))
  }
}


#' Centering a squared matrix by row or column
#'
#' It centers a numeric matrix with dimension \emph{N x N} by row (rows=TRUE) or column
#'  (rows=FALSE).
#'
#' @inheritParams cosnormX
#' @param rows If TRUE, the operation is done by row; otherwise, it is done by
#' column. (Defaults: TRUE).
#' @return Centered X (class "matrix").
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(25),ncol=5,nrow=5)
#' centerX(dat)

centerX <- function(X,rows=TRUE) {
  if(!rows) {
    m <- nrow(X)
    H <- diag(m) - matrix(rep(1/m,length(X)),nrow=nrow(X),ncol=ncol(X))
    return(H%*%X)
  } else {
    m <- ncol(X)
    H <- diag(m) - matrix(rep(1/m,length(X)),nrow=nrow(X),ncol=ncol(X))
    return(X%*%H)
  }
}


#' Minmax normalization
#'
#' Minmax normalization. Custom min/max values may be passed to the function.
#'
#' @inheritParams cosnormX
#' @param rows If TRUE, the minmax normalization is done by row; otherwise, it
#' is done by column. (Defaults: FALSE)
#' @param values (optional) A list containing two elements, the "max"
#' values and the "min" values. If no value is passed, the typical minmax normalization
#' (which normalizes the dataset between 0 and 1) is computed with the observed
#' maximum and minimum value in each column (or row) of X.
#' @return Minmax-normalized X.
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(100),ncol=10,nrow=10)
#' dat_minmax <- minmax(dat)
#' apply(dat_minmax,2,min) ## Min values = 0
#' apply(dat_minmax,2,max) ## Max values = 1
#' # We can also explicitly state the max and min values:
#' values <- list(min=apply(dat,2,min),max=apply(dat,2,max))
#' dat_minmax <- minmax(dat,values=values)

minmax <- function(X,rows=FALSE,values=NULL) {
  if(rows) {
    if(is.null(values)) {
      values <- list(min=apply(X,1,min), max=apply(X,1,max))
    } else {
      minmax_errors(values=values, X=X, rows=rows) # Check Errors
    }
  } else { # per column
    if(is.null(values)) {
      values <- list(min=apply(X,2,min), max=apply(X,2,max))
      } else {
        minmax_errors(values=values, X=X, rows=rows)  # Check Errors
      }
    values$min <- matrix(values$min,ncol=ncol(X),nrow=nrow(X),byrow = TRUE)
    values$max <- matrix(values$max,ncol=ncol(X),nrow=nrow(X),byrow = TRUE)
  }
  return((X-values$min)/(values$max-values$min))
}


#' Levels per factor variable
#'
#' This function gives the categories ("levels") per categorical variable ("factor").
#'
#' @param X A matrix, or a data.frame containing factors. (If the columns are of
#' any other class, they will be coerced into factors anyway).
#' @return A list with the levels.
#' @export
#' @importFrom methods is
#'
#' @examples
#' summary(showdata)
#' dummy_var(showdata)

dummy_var <- function(X)  {
  X <- toFactor(X)
  Xlev <-  lapply(X,levels)
  # if(length(Xlev)==1) Xlev <- Xlev$X
  return(Xlev)
}


#' Convert categorical data to dummies.
#'
#' Given a matrix or data.frame containing character/factors, this function
#' performs one-hot-encoding.
#'
#' @inheritParams dummy_var
#' @param lev (optional) A vector with the categories ("levels") of each factor.
#' @return X (class: "matrix") after performing one-hot-encoding.
#' @export
#' @importFrom methods hasArg
#'
#' @examples
#' summary(CO2)
#' CO2_dummy <- dummy_data(CO2[,1:3],lev=dummy_var(CO2[,1:3]))
#' CO2_dummy[1:10,1:5]

dummy_data <- function(X,lev=NULL) { # Convert X to dummies
  rowsnames <- rownames(X)
  X <- toFactor(X)
  if(!methods::hasArg(lev)) lev <- dummy_var(X)
  position <- 1:length(lev)
  if(length(position)!=ncol(X)) stop("Different number of variables in X and in lev")
  colsnames <-   unlist(lev)
  nlev <- sapply(lev,length)

  # colsnames <- paste(substr(names(colsnames),start = 1,stop=3),colsnames,sep="_")
  colsnames <- paste( rep(names(nlev),nlev),colsnames,sep="_")
  # lapply(colsnames,function(x)paste(names(x),x))

  X <- sapply(X,as.integer) #Level names to integers
  ddata <- matrix(0,ncol=sum(nlev),nrow=nrow(X))
  indices <- vector(mode="numeric",length=ncol(X)) ## Where each variable starts after binary expansion
  indices[1] <- 0
  for (i in position[-1]) {
    indices[i] <- indices[i-1] + nlev[i-1]
  }
  for (i in position) {
    for(j in 1:nlev[i]) {
      pos <- indices[i] + j
      ddata[,pos] <- X[,i]==j
    }
  }
  colnames(ddata) <- colsnames
  rownames(ddata) <- rowsnames
  return(ddata)
}


## Helpers

#' Helper minmax
#' @keywords internal
#' @noRd
minmax_errors <- function(values, X, rows) {
  if(length(values$min) != length(values$max)) stop("min and max values should have the same length")
  if(rows) {
    l <- nrow(X)
    dimn <- "nrow(X)"
  } else {
    l <- ncol(X)
    dimn <- "ncol(X)"
  }
  if(length(values$min) != l) stop(paste("min and max values should have length",dimn))
}
