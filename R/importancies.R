#######################################
###### Feature Importances Study ######
#######################################

### Useful to analyze feature importance


#' Importance barplot
#'
#' `plotImp()` displays the barplot of a numeric vector, which is assumed to contain the
#' features importance (from a prediction model) or the contribution of each
#' original variable to a Principal Component (PCA). In the barplot, features/PCs
#' are sorted by decreasing importance.
#'
#' @param x Numeric vector containing the importances.
#' @param y (optional) Numeric vector containing a different, independent variable to
#' be contrasted with the feature importances. Should have the same length and order
#' than `x`.
#' @param ylegend (optional) It allows to add a text explaining what is `y` (only
#' if `y` is not NULL).
#' @param leg_pos If `ylegend` is TRUE, the position of the legend. (Defaults: "right").
#' @param relative If TRUE, the barplot will display relative importances. (Defaults: TRUE).
#' @param absolute If FALSE, the bars may be positive or negative, which will affect
#' the order of the features Otherwise, the absolute value of `x` will be taken (Defaults: TRUE).
#' @param nfeat (optional) The number of top (most important) features displayed in the plot.
#' @param names (optional) The names of the features, in the same order than `x`.
#' @param main (optional) Plot title.
#' @param xlim (optional) A numeric vector. If absent, the minimum and maximum
#' value of `x` will be used to establish the axis' range.
#' @param color Color(s) chosen for the bars. A single value or a vector. (Defaults: "grey").
# @param labcolor (optional) Numeric vector containing the actual values.
#' @param leftmargin (optional) Left margin space for the plot.
#' @param ... (optional) Additional arguments (such as `axes`, `asp`,...) and graphical
#' parameters (such as `par`). See `?graphics::barplot()`.
#' @return A list containing:
#'
#' * The vector of importances in decreasing order. When `nfeat` is not NULL, only
#' the top `nfeat` are returned.
#'
#' * The cumulative sum of (absolute) importances.
#'
#' * A numeric vector giving the coordinates of all the drawn bars' midpoints.
#'
#' @importFrom methods hasArg
#' @importFrom graphics axis barplot legend lines par points
#' @export
#'
#' @examples
#' importances <- rnorm(30)
#' names_imp <- paste0("Feat",1:length(importances))
#'
#' plot1 <- plotImp(x=importances,names=names_imp,main="Barplot")
#' plot2 <- plotImp(x=importances,names=names_imp,relative=FALSE,
#' main="Barplot",nfeat=10)
#' plot3 <- plotImp(x=importances,names=names_imp,absolute=FALSE,
#' main="Barplot",color="coral2")

plotImp <- function(x,y=NULL, relative=TRUE, absolute=TRUE, nfeat=NULL,
                    names=NULL, main=NULL, xlim=NULL, color="grey",
                    leftmargin=NULL, ylegend=NULL, leg_pos="right",...) {

  original_mar <-  graphics::par()$mar
  if(!is.null(leftmargin)) {
    # if(is.null(y)) {
    graphics::par(mar = c(5, leftmargin, 1, 1))
    # } else {
    #   par(mar = c(5, 5, 1, 1))
    # }
  }

  x_abs <- abs(x)

  if(!is.null(names))   {
    names(x) <- names(x_abs) <- names
  } else {
    if(is.null(names(x))) {
      names(x) <- names(x_abs) <- names <- 1:length(x)
    } else {
      names(x) -> names
    }
  }

  x_den <- sum(x_abs)
  xord <- order(x_abs,decreasing = TRUE)
  x <- x[xord]
  x_abs <- x_abs[xord]

  if(!is.null(nfeat)) {
    if(nfeat>length(x)) stop("nfeat cannot be larger than x")
    x <- x[1:nfeat]
    x_abs <- x_abs[1:nfeat]
  }


  xord <- order(x_abs,decreasing = FALSE)


  if(relative) {
    x_abs <- x_abs/x_den
    x  <- x/x_den
  }

  if(methods::hasArg(y)) {
    if(methods::hasArg(names)) {
      names(y) <- names
    }
    if(relative) y <- y/sum(abs(y))
    if(absolute) y <- abs(y)
    y <- y[names(x_abs)]
    y <- y[xord]
  }

  if(length(color)==length(names))  {
    names(color) <- names
    color <- color[names(x_abs)]
    color <- color[xord]
  }

  if(absolute) {
    if(is.null(xlim)) xlim <- c(0,max(c(x_abs[xord],y))+0.1*max(c(x_abs[xord],y)))
    p <- graphics::barplot(x_abs[xord],las=2,horiz=TRUE, main=main,xlim=xlim, border=color, col=color,...)
    names <- names(x_abs[xord])
  } else {
    if(is.null(xlim)) xlim <- c(min(c(y,x[xord]))-0.1*abs(min(c(y,x[xord]))),max(c(x[xord],y))+0.1*max(c(x[xord],y)))
    signs <- c(sign(xlim[1]),sign(xlim[2]))
    if(all(signs==c(-1,-1))) xlim[2] <- 0
    if(all(signs==c(1,1))) xlim[1] <- 0
    p <- graphics::barplot(x[xord],las=2,horiz=TRUE, main=main,xlim=xlim, border=color, col=color,...)
    names <- names(x[xord])
  }
  if(methods::hasArg(y)) {
    graphics::lines(y = p, x= y,col="grey15")
    graphics::points(y =p,x=y,pch=20,col="grey15")
    if(methods::hasArg(ylegend)) graphics::legend(leg_pos, legend=ylegend, col="grey15", lty=6, cex=0.8,box.lty=0,inset=0.1)
  }

  # if(methods::hasArg(labcolor)) {
  #   for(j in 1:length(labcolor)) {
  #     ids <- labcolor[[j]]
  #     col <- names(labcolor)[j]
  #     ids <- names %in% ids
  #     graphics::axis(2, at = p[ids], labels = names[ids], col.axis = col, las=1)
  #   }
  # }

  toreturn <- list(first_features=names(x), cumsum=sum(x_abs),barplot = p)
  graphics::par(mar = original_mar)
  return(toreturn)
}
# labcolor <- paste0("<i style='color:",labcolor,"'>",names,"</i>")



#' Aggregate importances
#'
#' `aggregate_imp()` sums the importances present in a matrix or data.frame according
#' to some user-specified grouping criterion.
#'
#' @param X Matrix or data.frame containing the importances (in rows or in columns).
#' @param samples Samples are in rows or in columns? (Defaults: "rows").
#' @param lev (optional) The grouping elements. `lev` should be as long as the dimension (cols or rows)
#' that one wants to aggregate. If this parameter is absent, the colnames (if samples="rows") or rownames will
#' be used to that effect. In that case, it is expected that the col/rownames follow
#' this pattern: "V_Y", and the variables with the same "V" will be summed. (Check the colnames of
#' a typical output of `dummy_data()` for more info).
#' @return X, a matrix or data.frame containing the aggregated importances.
#'
#' @importFrom methods hasArg
#' @importFrom dplyr %>% group_by summarise_all
#' @export
#'
#' @examples
#' importances <- matrix(rnorm(90),nrow=3,ncol=30)
#' rownames(importances) <- c("sample1","sample2","sample3")
#' colnames(importances) <- paste0("Feat",
#' rep(1:5,times=2*(1:5)), "_", unlist(lapply(2*(1:5),function(x)LETTERS[1:x])))
#'
#' ## The grouping criterion is:
#' groups <- paste0("Feat",1:5)
#' aggregate_imp(X=importances,samples="rows",lev=groups)

#' ## We can also use the colnames:
#' colnames(importances)
#' aggregate_imp(X=importances,samples="rows")

aggregate_imp <- function(X,lev=NULL,samples="rows") {
  if(samples=="rows") X <- t(X)
  if(!hasArg(lev)) {
    lev <- sub('_.*', '', rownames(X))
  } else {
    ids <- sapply(lapply(paste0(lev,"_"),function(x)grep(x,rownames(X))),length)
    lev <- rep(lev,times= ids )
  }
  levs <- NULL
  group_id <- data.frame(X,levs=lev) %>% group_by(levs) %>% summarise_all(sum)
  group_id <- as.data.frame(group_id)
  colsnames <- group_id$levs
  group_id <- t(group_id[,-1])
  colnames(group_id) <- colsnames
  if(samples=="rows") {
      return(group_id[,unique(lev)])
  } else {
    return(t(group_id[,unique(lev)]))
  }
}
