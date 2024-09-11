#############################
### PERFORMANCE MEASURES ###
#############################

### Several simple and widespread measures for computing a model's performance


## Regression

#' NMSE (Normalized Mean Squared Error)
#'
#' `nmse()` computes the Normalized Mean Squared Error between the output
#' of a regression model and the actual values of the target.
#'
#' @details The Normalized Mean Squared error is defined as:
#'
#' \deqn{NMSE=MSE/((N-1)*var(target))}
#'
#' where MSE is the Mean Squared Error.
#'
#' @param target Numeric vector containing the actual values.
#' @param pred Numeric vector containing the predicted values.
#' (The order should be the same than in the target)
#'
#' @return The normalized mean squared error (a single value).
#'
#' @importFrom stats var
#' @importFrom methods is
#' @export
#'
#' @examples
#' y <- 1:10
#' y_pred <- y+rnorm(10)
#' nmse(y,y_pred)

nmse <- function(target,pred) {
  N <- length(target)

  # Errors
  if(!(methods::is(target,"numeric")&methods::is(pred,"numeric"))) stop("target and pred should be numeric")
  if(N<2) stop("target should have length > 1")
  if(N != length(pred)) stop("target and pred should have the same length")
  if(var(target)==0) stop("target should have non-zero variance")

  error <- sum((target-pred)^2)/((N-1)*var(target))
  return(error)
}


## Classification

#' Accuracy
#'
#' `Acc()` computes the accuracy between the output
#' of a classification model and the actual values of the target.
#' It can also compute the weighted accuracy, which is useful in
#' imbalanced classification problems. The weighting is applied according
#' to the class frequencies in the target. In balanced problems, weighted Acc = Acc.
#' @param ct Confusion Matrix.
#' @param weighted If TRUE, the weighted accuracy is returned. (Defaults: FALSE).
#' @return Accuracy of the model (a single value).
#' @examples
#' y <- c(rep("a",3),rep("b",2))
#' y_pred <- c(rep("a",2),rep("b",3))
#' ct <- table(y,y_pred)
#' Acc(ct)
#' Acc(ct,weighted=TRUE)
#' @importFrom methods is
#' @export

Acc <- function(ct,weighted=FALSE) {
  # Errors
  if(!(methods::is(ct,"table")|methods::is(ct,"matrix"))) stop("ct should be class matrix or table")
  if(nrow(ct)!=ncol(ct)) stop("ct should be squared")

  if(weighted) {
    weights <- rowSums(ct)
    return(sum(diag(ct)/weights/length(weights)))
  } else   {
    return(sum(diag(ct))/sum(ct))
  }
}


#' Accuracy of a random model
#'
#' `Acc_rnd()` computes the expected accuracy of a random classifier based on the
#' class frequencies of the target. This measure can be used as a benchmark when contrasted
#' to the accuracy (in test) of a given prediction model.
#' @param target A character vector or a factor. Alternatively, a numeric vector
#' (see below).
#' @param freq TRUE if `target` contains the frequencies of the classes (in this
#' case, `target` should be numeric), FALSE otherwise. (Defaults: FALSE).
#' @return Expected accuracy of a random classification model (a single value).
#' @examples
#' # Expected accuracy of a random model:
#' target <- c(rep("a",5),rep("b",2))
#' Acc_rnd(target)
#' # This is the same than:
#' freqs <- c(5/7,2/7)
#' Acc_rnd(freqs,freq=TRUE)
#' @export

Acc_rnd <- function(target,freq=FALSE) {
  # Errors
  if(length(target)<2) stop("target should have length > 1")
  if(freq) {
    if(!(methods::is(target,"numeric"))) stop("target should be numeric")
    if(sum(target)!=1) stop("frequencies should sum to 1")
  } else {
    # Code
    target <- table(target)/sum(table(target))
  }
  return(sum(target^2))
}


#' Confidence Interval using Normal Approximation
#'
#' `Normal_CI()` computes the Confidence Interval (CI) of a performance measure
#' (for instance, accuracy) using normal approximation. Thus, it is advisable
#' that the test has a size of, at least, 30 instances.
#' @param value Performance value (a single value).
#' @param ntest Test set size (a single value).
#' @param confidence Confidence level; for instance, 95\% or 99\%. (Defaults: 95).
#' @return A vector containing the CI.
#' @examples
#' # Computing accuracy
#' y <- c(rep("a",30),rep("b",20))
#' y_pred <- c(rep("a",20),rep("b",30))
#' ct <- table(y,y_pred)
#' accuracy <- Acc(ct)
#' # Computing 95%CI
#' Normal_CI(accuracy, ntest=length(y), confidence=95)
#' @importFrom stats qnorm
#' @export

Normal_CI <- function(value, ntest, confidence=95) {
  #Errors
  if(length(value)>1||length(ntest)>1) stop("value and ntest should have length==1")
  if(!(methods::is(value,"numeric")&&methods::is(ntest,"numeric"))) stop("value and ntest should be numeric")
  if(value>1||value<0) stop("value cannot be lower than 0 or greater than 1")
  if(ntest<30) stop("ntest should be 30 or greater")
  if(confidence>=100||confidence<=0) stop("confidence should be greater than 0 and lower than 100")

  confidence <- (1-confidence/100)
  z <- qnorm(p=confidence/2, lower.tail=FALSE)
  dev <- z*sqrt((value*(1-value))/ntest)
  return(c(value-dev,value+dev))
}


#'Confidence Interval using Bootstrap
#'
#' `Boots_CI()` computes the Confidence Interval (CI) of a performance measure
#' (for instance, accuracy) via bootstrapping.
#' @param target Numeric vector containing the actual values.
#' @param pred Numeric vector containing the predicted values.
#' (The order should be the same than the target's).
#' @param index Performance measure name, in lowercase. (Defaults: "acc").
#' @param nboots Number of bootstrapping replicas.
#' @param confidence Confidence level; for instance, 95\% or 99\%. (Defaults: 95).
#' @param ... Further arguments to be passed to the performance measures functions;
#' notably, multi.class="macro" or multi.class="micro" for the macro or micro
#' performance measures. (Defaults: "macro").
#' @return A vector containing the bootstrap estimate of the performance and its CI.
#' @examples
#' y <- c(rep("a",30),rep("b",20))
#' y_pred <- c(rep("a",20),rep("b",30))
#' # Computing Accuracy with their 95%CI
#' Boots_CI(target=y, pred=y_pred, index="acc", nboots=1000, confidence=95)
#' @importFrom stats quantile
#' @export

Boots_CI <- function(target, pred, index="acc", nboots, confidence=95,...) {
  # Errors
  if(length(target)!=length(pred)) stop("target and pred should have the same length")
  if(nboots<399) stop("please set a higher nboots")
  if(confidence>=100||confidence<=0) stop("confidence should be greater than 0 and lower than 100")

  confidence <- confidence/100
  quant <- c((1-confidence)/2,1-(1-confidence)/2)
  score_bootstrap <- rep(0,nboots)

  for(i in 1:nboots) {
    ids <- sample(1:length(target),replace = TRUE)
    y1_b <-  target[ids]
    y2_b <- pred[ids]
    if(index=="acc") {
      score_bootstrap[i] <- Acc(table(y1_b,y2_b),...)
    } else if(index=="nmse") {
      score_bootstrap[i] <- nmse(y1_b,y2_b)
    } else if(index=="prec") {
      score_bootstrap[i] <- Prec(table(y1_b,y2_b),...)
    } else if(index=="rec") {
      score_bootstrap[i] <- Rec(table(y1_b,y2_b),...)
    } else if(index=="spe") {
      score_bootstrap[i] <- Spe(table(y1_b,y2_b),...)
    } else if(index=="f1") {
      score_bootstrap[i] <- F1(table(y1_b,y2_b),...) ## nota: fa les versions "macro" o "micro"
    } else{
      stop("Measure not available")
    }
  }
  return(c(mean(score_bootstrap),quantile(score_bootstrap,quant)))
}


#' Precision or PPV
#'
#' `Prec()` computes the Precision of PPV (Positive Predictive Value) between the output
#' of a classification model and the actual values of the target.
#' The precision of each class can be aggregated. Macro-precision is the average of the
#' precision of each classes. Micro-precision is the weighted average.
#' @param ct Confusion Matrix.
#' @param multi.class Should the results of each class be aggregated, and how?
#' Options: "none", "macro", "micro". (Defaults: "macro").
#' @return PPV (a single value).
#' @examples
#' y <- c(rep("a",3),rep("b",2))
#' y_pred <- c(rep("a",2),rep("b",3))
#' ct <- table(y,y_pred)
#' Prec(ct)
#' @importFrom methods is
#' @export

Prec <- function(ct,multi.class="macro") {
  # Errors
  imbfunct_helper(ct,multi.class)

  # if(multi.class=="micro")     pr <- sum(diag(ct))/sum(colSums(ct))
  if(multi.class=="micro") return(acc_equivalent(ct))

  pr <- prec(ct,1:ncol(ct))
  isnan <- is.nan(pr)
  if(any(isnan)) pr[which(isnan)] <- 0

  if(multi.class=="none") {
    names(pr) <- colnames(ct)
    return(pr)
  }  else {
    message("It is identical to weighted Accuracy")
    return(mean(pr))
  }
}


#' Recall or Sensitivity or TPR
#'
#' `Rec()` computes the Recall, also known as Sensitivity or TPR (True Positive Rate),
#'  between the output of a classification model and the actual values of the target.
#' @inheritParams Prec
#' @return TPR (a single value).
#' @examples
#' y <- c(rep("a",3),rep("b",2))
#' y_pred <- c(rep("a",2),rep("b",3))
#' ct <- table(y,y_pred)
#' Rec(ct)
#' @importFrom methods is
#' @export

Rec <-  function(ct,multi.class="macro") {
  # Errors
  imbfunct_helper(ct,multi.class)

  # if(multi.class=="micro")     rc <- sum(diag(ct))/sum(rowSums(ct))
  if(multi.class=="micro") return(acc_equivalent(ct))

  rc <- rec(ct,1:ncol(ct))
  isnan <- is.nan(rc)
  if(any(isnan)) rc[which(isnan)] <- 0

  if(multi.class=="none") {
    names(rc) <- colnames(ct)
    return(rc)
  }  else {
    return(mean(rc))
  }
}


#' Specificity or TNR
#'
#' `Spe()` computes the Specificity or TNR (True Negative Rate) between the output
#' of a classification prediction model and the actual values of the target.
#' @inheritParams Prec
#' @return TNR (a single value).
#' @examples
#' y <- c(rep("a",3),rep("b",2))
#' y_pred <- c(rep("a",2),rep("b",3))
#' ct <- table(y,y_pred)
#' Spe(ct)
#' @importFrom methods is
#' @export

Spe <-  function(ct,multi.class="macro") {
  # Errors
  imbfunct_helper(ct,multi.class)

  if(multi.class=="micro") {
    tp <- sum(spe_mic(ct,1:ncol(ct)))
    sp <- tp/(sum(ct)*nrow(ct) - sum(rowSums(ct)))
  } else {
    sp <- spe(ct,1:ncol(ct))
  }

  isnan <- is.nan(sp)
  if(any(isnan)) sp[which(isnan)] <- 0

  if(multi.class=="none") {
    names(sp) <- colnames(ct)
    return(sp)
  }  else {
    return(mean(sp))
  }
}


#' F1 score
#'
#' `F1()` computes the F1 score between the output of a classification prediction model
#' and the actual values of the target.
#' @details F1 corresponds to the harmonic mean of Precision and Recall.
#'
#' @inheritParams Prec
#' @return F1 (a single value).
#' @examples
#' y <- c(rep("a",3),rep("b",2))
#' y_pred <- c(rep("a",2),rep("b",3))
#' ct <- table(y,y_pred)
#' F1(ct)
#' @importFrom methods is
#' @export

F1 <-  function(ct,multi.class="macro") {
  # Errors
  imbfunct_helper(ct,multi.class)

  if(multi.class=="micro") return(acc_equivalent(ct))

  REC <- Rec(ct,multi.class="none")
  PREC <- Prec(ct,multi.class="none")

  f1_score <- (2*PREC*REC)/(PREC+REC)
  isnan <- is.nan(f1_score)
  if(any(isnan)) f1_score[which(isnan)] <- 0
  if(multi.class=="macro") return(mean(f1_score))
  return(f1_score)
}


## Helpers

#' Check errors in Prec, Rec, Spe and F1
#' @keywords internal
#' @noRd
imbfunct_helper <- function(ct,multi.class) {
  if(!(methods::is(ct,"table")|methods::is(ct,"matrix"))) stop("ct should be class matrix or table")
  if(nrow(ct)!=ncol(ct)) stop("ct should be squared")
  if(!(multi.class %in% c("none","macro","micro","weighted"))) stop("Option not available")
}

#' A certain performance measure is equivalent to accuracy
#' @keywords internal
#' @noRd
acc_equivalent <- function(ct) {
    message("It is equivalent to usual Accuracy")
    class(ct)
    return(Acc(ct))
}

# Helpers (basic performance scores)
#' @keywords internal
#' @noRd
prec <- Vectorize(function(ct,min.class)ct[min.class,min.class]/sum(ct[,min.class]),vectorize.args = "min.class")
rec <-  Vectorize(function(ct,min.class)ct[min.class,min.class]/sum(ct[min.class,]),vectorize.args = "min.class")
spe <-  Vectorize(function(ct,min.class)sum(ct[-min.class,,drop=FALSE][-min.class])/sum(ct[-min.class,]),
                  vectorize.args = "min.class")
spe_mic <-  Vectorize(function(ct,min.class)sum(ct[-min.class,][-min.class]),  vectorize.args = "min.class")
