########################
### Kernel functions ###
########################


### Several standard and non-standard kernel functions that can be coupled
### to any kernel method.

# - Real vectors: Linear, RBF, Laplacian, Polynomial [not yet], Sigmoid [not yet]
# - Real matrices: Frobenius
# - Counts (absolute and relative frequencies): Ruzicka, Bray-Curtis
# - Categorical data: Overlap/Dirac
# - Sets: Intersect, Jaccard
# - Ordinal data: Kendall's tau
# ??? Chi-squared kernel? (histograms?)


## Kernel functions for real numbers

#' Linear kernel
#'
#' `Linear()` computes the inner product between all possible pairs of rows of a
#' matrix or data.frame with dimension \emph{NxD}.
#'
#' @param X Matrix or data.frame that contains real numbers ("integer", "float" or "double").
#' @param cos.norm  Should the resulting kernel matrix be cosine normalized? (Defaults: FALSE).
#' @param coeff (optional) A vector of length \emph{D} that weights each one of the
#' features (columns). When cos.norm=TRUE, `Linear()` first does the weighting and
#' then the cosine-normalization.
#'
#' @return Kernel matrix (dimension: \emph{NxN}).
#'
#' @export
#' @importFrom methods hasArg
#'
#' @examples
#' dat <- matrix(rnorm(250),ncol=50,nrow=5)
#' Linear(dat)

Linear <- function(X,cos.norm=FALSE,coeff=NULL) {
  X <- as.matrix(X)
  if(methods::hasArg(coeff)) X <- weight_helper(X,coeff)
  K <- tcrossprod(X)
  if(cos.norm) K <- cosNorm(K) # es pot fer més ràpid?
  colnames(K) <- rownames(K) <- rownames(X)
  return(K)
}


#' Gaussian RBF (Radial Basis Function) kernel
#'
#' `RBF()` computes the RBF kernel between all possible pairs of rows of a
#' matrix or data.frame with dimension \emph{NxD}.
#' @details Let \eqn{x_i,x_j} be two real vectors. Then, the RBF kernel is defined as:
#' \deqn{K_{RBF}(x_i,x_j)=\exp(-\gamma \|x_i - x_j \|^2)}
#'
#' Sometimes the RBF kernel is given a hyperparameter called sigma. In that case:
#' \eqn{\gamma = 1/\sigma^2}.
#'
#' @param X Matrix or data.frame that contains real numbers ("integer", "float" or "double").
#' @param g Gamma hyperparameter. If g=0 or NULL, `RBF()` returns the matrix of squared Euclidean
#' distances instead of the RBF kernel matrix.
#'
#' @return Kernel matrix (dimension: \emph{NxN}).
#'
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(250),ncol=50,nrow=5)
#' RBF(dat,g=0.1)

RBF <- function(X,g=NULL) { ## g = 1/sigma^2. g = NULL retorna la distància euclidiana al quadrat
  X <- as.matrix(X)
  N <- nrow(X)
  kk <- tcrossprod(X)
  dd <- diag(kk)
  D <- 2*kk-matrix(dd,N,N)-t(matrix(dd,N,N))
  colnames(D) <- rownames(D) <- rownames(X)
  if(is.null(g) ||g == 0 ) return(-D)
  return(exp(g*D))
}


#' Laplacian kernel
#'
#' `Laplace()` computes the laplacian kernel between all possible pairs of rows of a
#' matrix or data.frame with dimension \emph{NxD}.
#' @details Let \eqn{x_i,x_j} be two real vectors. Then, the laplacian kernel is defined as:
#' \deqn{K_{Lapl}(x_i,x_j)=\exp(-\gamma \|x_i - x_j \|_1)}
#'
#' @param X Matrix or data.frame that contains real numbers ("integer", "float" or "double").
#' @param g Gamma hyperparameter. If g=0 or NULL, `Laplace()` returns the Manhattan distance
#' (L1 norm between two vectors).
#'
#' @return Kernel matrix (dimension: \emph{NxN}).
#'
#' @export
#' @importFrom stats dist
#'
#' @examples
#' dat <- matrix(rnorm(250),ncol=50,nrow=5)
#' Laplace(dat,g=0.1)

Laplace <- function(X,g=NULL) { ## g = NULL retorna la distància de Manhattan
  X <- as.matrix(X)
  D <- as.matrix(dist(X,method="manhattan"))
  if(is.null(g) ||g == 0 ) return(D)
  return(exp(-g*D))
}



## Kernel functions for real matrices

#' Frobenius kernel
#'
#' `Frobenius()` computes the Frobenius kernel between numeric matrices.
#'
#' @details The Frobenius kernel is the same than the Frobenius inner product between
#' matrices.
#' @param DATA A list of \emph{M} matrices or data.frames containing only real
#' numbers (class "integer", "float" or "double").
#' All matrices or data.frames should have the same number of rows and columns.
#' @inheritParams Dirac
#' @inheritParams Linear
#' @return Kernel matrix (dimension:\emph{NxN}), or a list with the kernel matrix and the
#' feature space.
#' @export
#'
#' @examples
#'
#' data1 <- matrix(rnorm(250000),ncol=500,nrow=500)
#' data2 <- matrix(rnorm(250000),ncol=500,nrow=500)
#' data3 <- matrix(rnorm(250000),ncol=500,nrow=500)
#'
#' Frobenius(list(data1,data2,data3))

Frobenius <- function(DATA,cos.norm=FALSE, feat_space=FALSE) {
  if( (length(unique(sapply(DATA,ncol))) > 1) || (length(unique(sapply(DATA,nrow))) > 1) ) {
    stop("All matrices or data.frames should have the same dimensions")
  }

  DATA <- lapply(DATA,function(X)as.matrix(X))
  if(cos.norm) DATA <- lapply(DATA,frobNorm)

  similarity <- outer(1:length(DATA), 1:length(DATA),
                        FUN = Vectorize(function(i, j) {
                          return(sum(DATA[[i]]* DATA[[j]]))
                        }))
  rownames(similarity) <- colnames(similarity) <- names(DATA)
  if(feat_space) {
    return(list(K=similarity,feat_space=t(sapply(DATA,as.vector))))
  } else {
    return(similarity)
  }
}


## Kernel functions for frequencies (Non-negative real numbers)

#' Kernels for count data
#'
#' Ruzicka and Bray-Curtis are kernel functions for absolute or relative
#' frequencies and count data. Both kernels have as input a matrix or data.frame
#' with dimension \emph{NxD} and \emph{N}>1, \emph{D}>1, containing strictly non-negative real numbers.
#' Samples should be in the rows. Thus, when working with relative frequencies,
#' `rowSums(X)` should be 1 (or 100, or another arbitrary number) for \emph{all} rows
#' (samples) of the dataset.
#'
#' @details For more info about these measures, please check Details in
#' ?vegan::vegdist(). Note that, in the vegan help page, "Ruzicka" corresponds to
#' "quantitative Jaccard". `BrayCurtis(X)` gives the same result than
#'  `1-vegan::vegdist(X,method="bray")`, and the same with `Ruzicka(data)` and
#'  `1-vegan::vegdist(data,method="jaccard")`.
#'
#' @param X Matrix or data.frame that contains absolute or relative frequencies.
#' @return Kernel matrix (dimension: \emph{NxN}).
#'
#' @export
#' @examples
#' data <- matrix(rpois(5000,lambda=3),ncol=100,nrow=50)
#' Kruz <- Ruzicka(data)
#' Kbray <- BrayCurtis(data)
#' Kruz[1:5,1:5]
#' Kbray[1:5,1:5]

BrayCurtis <- function(X) return(freqkerns(X, kernel="bray"))

#' @rdname BrayCurtis
#' @export
Ruzicka <- function(X)  return(freqkerns(X=X, kernel="ruzicka"))


## Kernel functions for categorical data (factors) and sets

#' Kernels for categorical variables
#'
#' @description
#' From a matrix or data.frame with dimension \emph{NxD}, where \emph{N}>1, \emph{D}>0,
#' `Dirac()` computes the simplest kernel for categorical data. Samples
#' should be in the rows and features in the columns. When there is a single feature,
#' `Dirac()` returns 1 if the category (or class, or level) is the same in
#' two given samples, and 0 otherwise. Instead, when \emph{D}>1, the results for the
#' \emph{D} features are combined doing a sum, a mean, or a weighted mean.
#'
#' @references
#' Belanche, L. A., and Villegas, M. A. (2013).
#' Kernel functions for categorical variables with application to problems in the life sciences.
#' Artificial Intelligence Research and Development (pp. 171-180). IOS Press.
#' \href{https://upcommons.upc.edu/bitstream/handle/2117/23347/KernelCATEG_CCIA2013.pdf}{Link}
#'
#' @param X Matrix (class "character") or data.frame (class "character", or columns = "factor").
#' The elements in X are assumed to be categorical in nature.
# @param na_value (optional) Treatment for NAs during the comparison.
# For instance, na_action = FALSE will set all NAs to FALSE (that is: 0).
#' @param comp When \emph{D}>1, this argument indicates how the variables
#' of the dataset are combined. Options are: "mean", "sum" and "weighted". (Defaults: "mean")
#' \itemize{
#'   \item "sum" gives the same importance to all variables, and returns an
#'   unnormalized kernel matrix.
#'   \item "mean" gives the same importance to all variables, and returns a
#'   normalized kernel matrix (all its elements range between 0 and 1).
#'   \item "weighted" weights each variable according to the `coeff` parameter, and returns a
#'   normalized kernel matrix.
#' }
#' @param coeff (optional) A vector of weights with length \emph{D}.
#' @param feat_space If FALSE, only the kernel matrix is returned. Otherwise,
#' the feature space is also returned. (Defaults: FALSE).
#' @return Kernel matrix (dimension: \emph{NxN}), or a list with the kernel matrix and the
#' feature space.
#'
#' @export
#' @importFrom methods is
#' @examples
#' # Categorical data
#' summary(CO2)
#' Kdirac <- Dirac(CO2[,1:3])
#' ## Display a subset of the kernel matrix:
#' Kdirac[c(1,15,50,65),c(1,15,50,65)]

Dirac <- function(X, comp="mean", coeff=NULL,feat_space=FALSE) {
  return(catkerns(X=X, kernel="dirac", comp=comp, coeff=coeff,feat_space=feat_space))
}


#' Kernels for sets
#'
#' @description
#' `Intersect()` or `Jaccard()` compute the kernel functions of the same name,
#' which are useful for set data. Their input is a matrix or data.frame with
#' dimension \emph{NxD}, where \emph{N}>1, \emph{D}>0. Samples should be in the
#' rows and features in the columns. When there is a single feature,
#' `Jaccard()` returns 1 if the elements of the set are exactly the same in
#' two given samples, and 0 if they are completely different (see Details). Instead,
#' in the multivariate case (\emph{D}>1), the results (for both `Intersect()` and
#' `Jaccard()`) of the \emph{D} features are combined with a sum, a mean, or a
#' weighted mean.
#'
#' @details Let \eqn{A,B} be two sets. Then, the Intersect
#' kernel is defined as:
#'
#' \deqn{K_{Intersect}(A,B)=|A \cap B| }
#'
#' And the Jaccard kernel is defined as:
#'
#' \deqn{K_{Jaccard}(A,B)=|A \cap B| / |A \cup B|}
#'
#' This specific implementation of the Intersect and Jaccard kernels expects
#' that the set members (elements) are character symbols (length=1). In case the
#' set data is multivariate (\emph{D}>1 columns, and each one contains a set feature),
#' elements for the \emph{D} sets should come from the same domain (universe).
#' For instance, a dataset with two variables, so the elements
#' in the first one are colors c("green","black","white","red") and the second are names
#' c("Anna","Elsa","Maria") is not allowed. In that case, set factors should be recoded
#' to colors c("g","b","w","r") and names c("A","E","M") and, if necessary, 'Intersect()'
#' (or `Jaccard()`) should be called twice.
#'
#' @references
#' Bouchard, M., Jousselme, A. L., and Doré, P. E. (2013).
#' A proof for the positive definiteness of the Jaccard index matrix.
#' International Journal of Approximate Reasoning, 54(5), 615-626.
#'
#' Ruiz, F., Angulo, C., and Agell, N. (2008).
#' Intersection and Signed-Intersection Kernels for Intervals.
#' Frontiers in Artificial Intelligence and Applications. 184. 262-270.
#' doi: 10.3233/978-1-58603-925-7-262.
#'
#' @param elements All potential elements (symbols) that can appear in the sets. If there are
#' some elements that are not of interest, they can be excluded so they are not
#' taken into account by these kernels. (Defaults: LETTERS).
# @param na_value (optional) Treatment for NAs during the comparison.
# For instance, na_action = FALSE will set all NAs to FALSE (that is: 0).
#' @param feat_space (not available for the Jaccard kernel). If FALSE, only the
#' kernel matrix is returned. Otherwise, the feature space is returned too. (Defaults: FALSE).
#' @inheritParams Dirac
#' @return Kernel matrix (dimension: \emph{NxN}), or a list with the kernel matrix and the
#' feature space.
#'
#' @export
#' @importFrom methods is
#' @examples
#' # Sets data
#' ## Generating a dataset with sets containing uppercase letters
#' random_set <- function(x)paste(sort(sample(LETTERS,x,FALSE)),sep="",collapse = "")
#' max_setsize <- 4
#' setsdata <- matrix(replicate(20,random_set(sample(2:max_setsize,1))),nrow=4,ncol=5)
#'
#' ## Computing the Intersect kernel:
#' Intersect(setsdata,elements=LETTERS,comp="sum")
#'
#' ## Computing the Jaccard kernel weighting the variables:
#' coeffs <- c(0.1,0.15,0.15,0.4,0.20)
#' Jaccard(setsdata,elements=LETTERS,comp="weighted",coeff=coeffs)

Jaccard <- function(X, elements=LETTERS, comp="mean", coeff=NULL) {
  return(catkerns(X=X, elements=elements, kernel="jaccard", comp=comp, coeff=coeff,feat_space=FALSE))
}

#' @rdname Jaccard
#' @export
Intersect <- function(X, elements=LETTERS,  comp="mean", coeff=NULL,feat_space=FALSE) {
  return(catkerns(X=X, elements=elements, kernel="intersect", comp=comp, coeff=coeff,feat_space=feat_space))
}


## Kernel functions for strings (character sequences)

#' Spectrum kernel
#'
#' `Spectrum()` computes the basic Spectrum kernel between strings. This kernel
#' computes the similarity of two strings by counting how many matching substrings
#' of length \emph{l} are present in each one.
#'
#' @details In large datasets this function may be slow. In that case, you may use the `stringdot()`
#' function of the `kernlab` package, or the `spectrumKernel()` function of the `kebabs` package.
#'
#' @references Leslie, C., Eskin, E., and Noble, W.S.
#' The spectrum kernel: a string kernel for SVM protein classification.
#' Pac Symp Biocomput. 2002:564-75. PMID: 11928508.
#' \href{http://psb.stanford.edu/psb-online/proceedings/psb02/abstracts/p564.html}{Link}
#'
#' @param x Vector of strings (length \emph{N}).
#' @param alphabet  Alphabet of reference.
#' @param l Length of the substrings.
#' @param group.ids (optional) A vector with ids. It allows to compute the kernel
#' over groups of strings within x, instead of the individual strings.
#' @param weights (optional) A numeric vector as long as x. It allows to weight differently
#' each one of the strings.
#' @param feat_space If FALSE, only the kernel matrix is returned. Otherwise,
#' the feature space (i.e. a table with the number of times that a substring of
#' length \emph{l} appears in each string) is also returned (Defaults: FALSE).
#' @inheritParams Linear
#'
#' @return Kernel matrix (dimension: \emph{NxN}), or a list with the kernel matrix and the
#' feature space.
#'
#' @export
#' @importFrom dplyr %>% group_by id summarise_all
#' @importFrom stringi stri_count
#'
#' @examples
#' ## Examples of alphabets. _ stands for a blank space, a gap, or the
#' ## start or the end of sequence)
#' NT <- c("A","C","G","T","_") # DNA nucleotides
#' AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T",
#' "V","W","Y","_") ##canonical aminoacids
#' letters_ <- c(letters,"_")
#' ## Example of data
#' strings <- c("hello_world","hello_word","hola_mon","kaixo_mundua",
#' "saluton_mondo","ola_mundo", "bonjour_le_monde")
#' names(strings) <- c("english1","english_typo","catalan","basque",
#' "esperanto","galician","french")
#' ## Computing the kernel:
#' Spectrum(strings,alphabet=letters_,l=2)

Spectrum <- function(x, alphabet, l=1, group.ids=NULL, weights=NULL, feat_space=FALSE, cos.norm = FALSE) {
  if(l>1) alphabet <- permute_rep(alphabet=alphabet,l=l)
  count_table <- searchSubs(x, fixed = alphabet, overlap=TRUE)
  if(!is.null(weights)) {
    if(length(weights) != length(x)) stop("weights length should be the same than x, the string vector provided")
    count_table <- weights * count_table
  }
  if(!is.null(group.ids)) {
    if(length(group.ids) != length(x)) stop("Ids length should be the same than x, the string vector provided")
    group_id <- data.frame(id = as.factor(group.ids), count_table)
    colnames(group_id)[-1] <- colnames(count_table)
    group_id <- group_id %>% group_by(id) %>% summarise_all(sum)
    count_table <- as.matrix(group_id[,-1])
    # count_table <- apply(group_id,2,as.numeric)
    rownames(count_table) <- group_id$id
  } else {
    rownames(count_table) <- names(x)
  }
  K <- Linear(X=count_table,cos.norm=cos.norm)
  if(feat_space) {
    if(cos.norm) count_table <- cosnormX(count_table)
    return(list(K=K,feat_space=count_table))
  } else {
    return(K)
  }
}


## Kernel functions for ordinal data (rankings)

#' Kendall's tau kernel
#'
#' `Kendall()` computes the Kendall's tau, which happens to be a kernel function
#' for ordinal variables, ranks or permutations.
#'
#' @references Jiao, Y. and Vert, J.P.
#' The Kendall and Mallows kernels for permutations. International Conference on Machine Learning.
#' PMLR, 2015. \href{https://proceedings.mlr.press/v37/jiao15.html}{Link}
#'
#' @param X When evaluating a single ordinal feature, X should be a numeric matrix
#' or data.frame. If data is multivariate, X should be a list, and each ordinal/ranking
#' feature should be placed in a different element of the list (see Examples).
#' @param NA.as.0 Should NAs be converted to 0s? (Defaults: TRUE).
#' @param samples.in.rows If TRUE, the samples are considered to be in the rows.
#' Otherwise, it is assumed that they are in the columns. (Defaults: FALSE).
#' @param comp If X is a list, this argument indicates how the ordinal/ranking variables
#' are combined. Options are: "mean" and "sum". (Defaults: "mean").
#'
#' @return Kernel matrix (dimension: \emph{NxN}).
#'
#' @export
#' @importFrom stats cor

#' @examples
#' # 3 people are given a list of 10 colors. They rank them from most (1) to least
#' # (10) favorite
#' color_list <-  c("black","blue","green","grey","lightblue","orange","purple",
#' "red","white","yellow")
#' survey1 <- 1:10
#' survey2 <- 10:1
#' survey3 <- sample(10)
#' color <- cbind(survey1,survey2,survey3) # Samples in columns
#' rownames(color) <- color_list
#' Kendall(color)
#'
#' # The same 3 people are asked the number of times they ate 5 different kinds of
#' # food during the last month:
#' food <- matrix(c(10, 1,18, 25,30, 7, 5,20, 5, 12, 7,20, 20, 3,22),ncol=5,nrow=3)
#' rownames(food) <- colnames(color)
#' colnames(food) <- c("spinach", "chicken", "beef" , "salad","lentils")
#' # (we can observe that, for person 2, vegetables << meat, while for person 3
#' # is the other way around)
#' Kendall(food,samples.in.rows=TRUE)
#'
#' # We can combine this results:
#' dataset <- list(color=color,food=t(food)) #All samples in columns
#' Kendall(dataset)

Kendall <-  function(X, NA.as.0=TRUE,samples.in.rows=FALSE,comp="mean")  UseMethod("Kendall",X)

#'@exportS3Method
Kendall.default  <- function(X, NA.as.0=TRUE,samples.in.rows=FALSE,comp="mean") {
  if(samples.in.rows) X <- t(X)
  if(NA.as.0) X[is.na(X)] <- 0
  stats::cor(X,method="kendall",use="all.obs")
}

#'@exportS3Method
Kendall.list <- function(X,NA.as.0=TRUE,samples.in.rows=FALSE,comp="mean")  {
  if(samples.in.rows) {
    if( length(unique(sapply(X,nrow))) > 1 ) {
      stop("All list's elements should have the same number of rows")
    }
  } else {
    if( length(unique(sapply(X,ncol))) > 1 ) {
      stop("All list's elements should have the same number of columns")
    }
  }

  K <- lapply(X,Kendall)
  K <- array(unlist(K), c(dim(K[[1]]), length(K)))
  if(comp =="mean") {
    cat("Composition: Mean",sep="\n")
    return(rowMeans(K, dims = 2))
  }
  else if(comp == "sum") {
    cat("Composition: Sum",sep="\n")
    return(rowSums(K, dims = 2))
  }
  else {
    stop(paste("Option not available."))
  }
}



## Helpers

#' Helper for computing weights in some kernels
#' @keywords internal
#' @noRd
weight_helper <- function(X,coeff) {
  dnames <- dimnames(X)
  coeff <- sqrt(coeff)
  if(methods::is(X,"matrix") ) {
    if(length(coeff)!=ncol(X)) stop("length(coeff) should be equal to ncol(X)")
    X <- X %*% diag(coeff)
    dimnames(X) <- dnames
    return(X)
  } else {
    D <- dim(X)[3]
    for(i in 1:D ) { ## això igual es pot millorar
      if(length(coeff)!=D) stop("length(coeff) should be equal to ncol(X)")
      X[,,i] <- X[,,i]*coeff[i]
    }
    dimnames(X) <- dnames
    return(X)
  }
}


#' Helper for constructing some kernel functions. x is a vector
#' @keywords internal
#' @noRd
expand.grid.mod <- function(x, rep=FALSE) {
  g <- function(i) {
    z <- setdiff(x, x[seq_len(i-rep)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


#' Permutations with repetition
#' @keywords internal
#' @noRd
permute_rep <- function(alphabet,l) {
  topermute <- as.data.frame(matrix(rep(alphabet,l),nrow=length(alphabet),ncol=l))
  return(apply(expand.grid(topermute),1,function(x)paste(x,collapse = "")))
}


#' Helper for the Spectrum kernel
#' @keywords internal
#' @noRd
searchSubs <- Vectorize(stringi::stri_count,vectorize.args = "fixed")


#' "core" bray-curtis function
#' @keywords internal
#' @noRd
braycurtis <- function (DATA,i) {
  I <- rowSums(abs(DATA[i[,1], ] - DATA[i[,2], ])) #Comparison matrix, with dimension: (n^2-n)/2 * d)
  U <- rowSums(DATA[i[,1], ] + DATA[i[,2], ])
  return(1-(I/U))
}


#' "core" weighted jaccard/ruzicka function
#' @keywords internal
#' @noRd
ruzicka <- function (DATA,i) {
  I <- rowSums(pmin(DATA[i[,1], ], DATA[i[,2], ]))
  U <- rowSums(pmax(DATA[i[,1], ], DATA[i[,2], ]))
  return(I/U)
}


#' "core" overlap function
#' @keywords internal
#' @noRd
overlap <- function(DATA, i) {
  comparacio <- DATA[i[,1], ] == DATA[i[,2], ] #Comparison matrix, with dimension: (n^2-n)/2 * d)
  return(comparacio)
}


#' "core" intersect function
#' @keywords internal
#' @noRd
intersec <- function(array,i,col) { #array is a three-dimensional array, and i are the indexes
  D <- 1:col
  Int <- matrix(nrow=nrow(i),ncol=col)
  for(j in D)    Int[,j] <- rowSums(array[i[,1],,j] * array[i[,2],,j])
  return(Int)
}


#' "core" jaccard function
#' @keywords internal
#' @noRd
jaccard <- function(array,i,col) { #array is a three-dimensional array, and i are the indexes
  D <- 1:col
  Jac <- matrix(nrow=nrow(i),ncol=col)
  for(j in D) {
    I <-   array[i[,1],,j] & array[i[,2],,j] # Intersection of binary vectors
    U <-   array[i[,1],,j] | array[i[,2],,j] # Union of binary vectors
    Jac[,j] <- rowSums(I)/rowSums(U) ## JACCARD INDEX: size of the I / size of the U.
  }
  return(Jac)
}


#' Jaccard & Intersect kernel auxiliar functions
#'
#' Assigns to each amino acid a number between 1 and 20
#' @keywords internal
#' @noRd
Let.to.Num <- function(M,alphabet) { # M is a matrix
  col <- ncol(M)
  row <- nrow(M)
  return(matrix(match(M,alphabet),ncol=col,nrow=row))
}

#' Converts the integer i (between 1 and 20) into a vector of 0s with a 1 at the i-th position
#' @keywords internal
#' @noRd
To.Logical.Vector <- function(i,lengthCode) { # x is an integer
  V <- vector(mode="numeric",length=lengthCode)
  V[i] <- 1
  return(V)
}

#' Converts the M matrix (n rows x d cols) into a 3D array (n rows x d cols x alphabet)
#' @keywords internal
#' @noRd
matrix3D <- function(M,members) {  # M is a matrix
  z <- max(nchar(M))
  lengthCode <- length(members)
  result <- array(0,dim=c(lengthCode,nrow(M),ncol(M)), dimnames=list(members,NULL,colnames(M)))
  for(i in 1:z) {   # Split matrix in z submatrices (allele mixtures)
    DATAi <- substr(M, start=i, stop=i)
    DATAi.NUM <- Let.to.Num(DATAi,alphabet=members) # Letter code to Number code
    DATAi.LOG <- apply(DATAi.NUM,c(1,2),function(x)To.Logical.Vector(x,lengthCode)) #Number to binary vector
    result <- DATAi.LOG + result
  }

  ### Delete categories without occurrences
  # index <- rowSums(result)
  # index <- which(index != 0)
  # result <- result[index,,]
  return(result) #Sum of binary vectors
}


#' Kernel Jaccard/Ruzicka or Bray-Curtis Kernels wrapper
#' @keywords internal
#' @noRd
freqkerns <- function(X,kernel="ruzicka") {
  data <- as.matrix(X)
  if(nrow(data)<2) stop("X should be a matrix or data.frame with at least two rows")
  ## Number of rows, and comparison id's.
  n <- nrow(data)
  id <- expand.grid.mod(1:n)

  if(kernel =="ruzicka") {
    Composicio <- ruzicka(DATA=data,i=id)
  } else if(kernel=="bray"){
    Composicio <- braycurtis(DATA=data,i=id)
    } else {
    stop(paste("Kernel not available"))
  }

  ## Building the kernel matrix
  Ncomb <- 1:((n^2-n)/2)
  K <- matrix(0,ncol = n,nrow = n)
  colnames(K) <- 1:n
  rownames(K) <- colnames(K)
  for (i in Ncomb)  K[id[i,1],id[i,2]] <- Composicio[i] # Upper triangular matrix
  tK <- t(K)
  K <- tK + K # Upper triangular matrix to symmetric matrix
  diag(K) <- 1
  return(K)
}



#' Dirac/Jaccard/Intersect Kernels wrapper
#' @keywords internal
#' @noRd
catkerns <- function (X, kernel="dirac",  comp="mean", coeff=NULL,
                      elements=LETTERS,feat_space=FALSE) {
  data <- as.matrix(X)
  n <- nrow(data)
  d <- ncol(data)

  # Errors
  if(n<2) stop("X should be a matrix or data.frame with at least two rows")
  if(comp=="weighted") {
    if(is.null(coeff)) stop("Please provide weights")
    if(length(coeff) != d) stop(paste("Number of weights and number of variables do not coincide."))
    if(sum(coeff)>1) coeff <- coeff/sum(coeff) ## Transform an absolute importance value into a relative one
  } else {
    coeff <- rep(1/d,d)
  }

  ## Comparison id's.
  id <- expand.grid.mod(1:n,rep=TRUE)

  # 1. Basic kernels
  if(kernel == "dirac") {
    if(feat_space) {
      levs <- dummy_var(data)
      DATA.LOG <- dummy_data(data,levs)
      coeff2 <- rep(coeff,as.numeric(sapply(levs,length)))
    }
    Comparacio <- overlap(DATA=data,i=id)
  } else if(kernel =="intersect") {
    DATA.LOG <- matrix3D(data,members=elements)
    DATA.LOG <- aperm(DATA.LOG,c(2,1,3)) #Transpose to obtain: 1st index: individuals, 2nd index: alphabet and 3rd index: feature position
    if(comp!="sum") for(j in 1:d) DATA.LOG[,,j] <- cosnormX(DATA.LOG[,,j]) # cosnorm so comparisons are between 0 and 1, like overlap and jaccard. in sum, implicitly more weight is given to variables with more elements
    Comparacio <- intersec(array=DATA.LOG,i=id, col=d)
    coeff2 <- coeff
  }  else if(kernel =="jaccard") {
    DATA.LOG <- matrix3D(data,members=elements)
    DATA.LOG <- aperm(DATA.LOG,c(2,1,3)) #Transpose to obtain: 1st index: individuals, 2nd index: alphabet and 3rd index: feature position
    Comparacio <- jaccard(array=DATA.LOG,i=id, col=d) # Jaccard
    DATA.LOG <- NA
  }  else {
    stop(paste("Kernel not available"))
  }

  # 2. Composition step
  if(methods::is(Comparacio,"matrix")) {
    if(comp =="mean") {
      Composicio <- rowMeans(Comparacio)
      if(feat_space)  DATA.LOG <- weight_helper(DATA.LOG,coeff=coeff2)
    } else if(comp == "sum") {
      Composicio <- rowSums(Comparacio)
    } else if (comp == "weighted") {
      tComp <- t(Comparacio) * coeff
      Comparacio <- t(tComp)
      rm(tComp)
      Composicio <- rowSums(Comparacio)
      if(feat_space) DATA.LOG <- weight_helper(DATA.LOG,coeff2)
    } else {
      stop(paste("Option not available."))
    }
  } else {
    Composicio <- Comparacio
    rm(Comparacio)
  }

  ## 3. Building the kernel matrix
  Ncomb <- 1:((n^2+n)/2)
  K <- matrix(0,ncol = n,nrow = n)
  colnames(K) <- 1:n
  rownames(K) <- colnames(K)
  for (i in Ncomb)  K[id[i,1],id[i,2]] <- K[id[i,2],id[i,1]] <- Composicio[i]
  # if(hasArg(g)) K <- exp(g*K)/exp(g)
  if(feat_space) {
    return(list(K=K,feat_space=DATA.LOG))
  } else {
    return(K)
  }
}
