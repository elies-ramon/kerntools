##################
### Kernel PCA ###
##################

### Kernel PCA and related functions


#' Kernel PCA
#'
#' @description
#'`kPCA()` computes the kernel PCA from a kernel matrix and, if desired, produces
#' a plot. The contribution of the original variables to the Principal Components (PCs),
#' sometimes referred as "loadings", is NOT returned (to do so, go to `kPCA_imp()`).
#'
#' @details As the ordinary PCA, kernel PCA can be used to summarize, visualize and/or
#' create new features of a dataset. Data can be projected in a linear or nonlinear
#' way, depending on the kernel used. When the kernel is `Linear()`, kernel PCA
#' is equivalent to ordinary PCA.
#'
#' @param K Kernel matrix (class "matrix").
#' @param center A logical value. If TRUE, the variables are zero-centered before
#' the PCA. (Defaults: TRUE).
#' @param Ktest (optional) An additional kernel matrix corresponding to test samples,
#' with dimension \emph{Ntest x Ntraining}. These new samples are projected
#' (using the color defined by `na_col`) over the kernel PCA computed from K.
#' Remember than the data that generated `Ktest` should be centered beforehand, using
#' the same values used for centering `K`.
#' @param plot (optional) A `ggplot2` is displayed. The input should be a vector of
#' integers with length 2, corresponding to the two Principal Components to be displayed in the plot.
#' @param y (optional) A factor, or a numeric vector, with length equal to `nrow(K)`
#' (number of samples). This parameter allows to paint the points with different colors.
#' @param colors A single color, or a vector of colors. If `y` is numeric, a gradient of colors
#' between the first and the second entry will be used to paint the points. (Defaults: "black").
#' @param na_col Color of the entries that have a NA in the parameter `y`, or the entries
#' corresponding to `Ktest` (when `Ktest` is not NULL). Otherwise, this parameter is ignored.
#' @param title Plot title.
#' @param pos_leg Position of the legend.
#' @param name_leg Title of the legend. (Defaults: blank)
#' @param labels (optional) A vector of the same length than nrow(K). A name will be
#' displayed next to each point.
#' @param ellipse (optional) A float between 0 and 1. An ellipse will be drawn for
#' each group of points defined by `y`. Here `y` should be of class "factor." This
#' parameter will indicate the spread of the ellipse.
#' @return A list with two objects:
#'
#' * The PCA projection (class "matrix"). Please note that if K was computed from a \emph{NxD}
#' table with \emph{N > D}, only the first \emph{N-D} PCs may be useful.
#'
#' * (optional) A `ggplot2` plot of the selected PCs.
#' @export
#' @import ggplot2
#' @importFrom methods hasArg
#' @importFrom stats qchisq
#' @importFrom methods is
#' @examples
#' dat <- matrix(rnorm(150),ncol=50,nrow=30)
#' K <- Linear(dat)
#'
#' ## Projection's coordinates only:
#' pca <- kPCA(K)
#'
#' ## Coordinates + plot of the two first principal components (PC1 and PC2):
#' pca <- kPCA(K,plot=1:2, colors = "coral2")
#' pca$plot


kPCA <- function(K, center=TRUE, Ktest=NULL, plot=NULL, y=NULL, colors="black", na_col="grey70",
                 title="Kernel PCA", pos_leg="right",name_leg="",labels=NULL,ellipse=NULL) {

  #kPCA
  if(center) K <- centerK(K)
  spectral_K <- svd(K)
  # La projecció és UD
  spectral_K$d <- sqrt(spectral_K$d)
  projected <- as.data.frame(spectral_K$u %*% diag(spectral_K$d))
  rownames(projected) <- rownames(K)
  colnames(projected) <- paste0("PC",1:ncol(projected))

  #Projecció d'un test (si n'hi ha)
  if(hasArg(Ktest)) {
    # La predicció d'un nou punt és kUD^(-1) on k=k(x,x_i)
    projectest <- Ktest %*% spectral_K$u %*% diag(1/spectral_K$d)
    colnames(projectest) <- colnames(projected)
    projected <- rbind(projected,projectest)
    rownames(projected) <- c(rownames(K),rownames(Ktest))
    if(is.null(y)) {
      y <- rep(1,nrow(K))
      pos_leg <- "none"
    }
    y2 <- c(y,rep(NA,nrow(Ktest)))
    names(y2) <- rownames(projected)
    if(methods::is(y,"factor"))   y <- factor(y2,levels=levels(factor(y2)),labels=levels(y))
  }

  #Plot
  if(hasArg(plot)) {
    i <- plot[1]
    j <- plot[2]
    spectral_K$d <- spectral_K$d^2 ##calculat com a prcomp
    xpercent <-  spectral_K$d[i]/sum(spectral_K$d)*100
    ypercent <-  spectral_K$d[j]/sum(spectral_K$d)*100
    q <-  ggplot2::ggplot(projected, aes(projected[,i], projected[,j])) +  theme_bw() +
      ggtitle(title) + xlab(paste0("PC ", i,": ", format(xpercent,digits=3), "%")) +
      ylab(paste0("PC ", j,": ", format(ypercent,digits=3), "%")) +
      # coord_fixed() +
      ggplot2::theme(legend.position = pos_leg,legend.box = "horizontal",plot.title = element_text(hjust = 0.5,face = "bold"))
    if(hasArg(y)) {
      q <- q + ggplot2::geom_point(ggplot2::aes(colour = y))
      if(hasArg(colors)) {
        if(methods::is(y,"factor")) {
          q <- q + ggplot2::scale_color_manual(name=name_leg,values =  colors,na.value = na_col)
        } else {
          q <- q + ggplot2::scale_colour_gradientn(name=name_leg, colours=colors, na.value = na_col)

        }
      } else {
        q$labels$colour <-  name_leg
      }
    } else if(hasArg(colors)) {
      q <-  q + ggplot2::geom_point(colour = colors)
    } else {
      q <- q +  ggplot2::geom_point()
    }
    if(hasArg(labels)) q <- q + ggplot2::geom_text(ggplot2::aes(label=labels,hjust=0.5, vjust=2))

    if((is.factor(y)|hasArg(Ktest)) && hasArg(ellipse)) {
      if(hasArg(Ktest)) {
        y[!is.na(y)] <- "Training"
        y[is.na(y)] <- "Test"
        y <- as.factor(y)
      }
      ell <- tapply(data.frame(xvar=projected[,i],yvar=projected[,j],groups=y),y,  function(X)ellipse_helper(X,ellipse=ellipse))
      ell <- data.frame(do.call(rbind, ell))
      colnames(ell)[1:3] <- c("xvar", "yvar","groups")
      xvar <- yvar <- groups <- NULL
      q <- q + geom_path(data = ell, aes(xvar,yvar,group=groups,color=groups))
      # +   coord_fixed()
    }
    return(list(projection=projected,plot=q))
  }
  return(projected)
}


#' Contributions of the variables to the Principal Components ("loadings")
#'
#' @description
#'`kPCA_imp()` performs a PCA and a kernel PCA simultaneously and returns
#' the contributions of the variables to the Principal Components (sometimes, these
#' contributions are called "loadings") in Feature Space. Optionally, it can also
#' return the samples' projection (cropped to the relevant PCs) and the values used
#' to centering the variables in Feature Space.
#' It does not return any plot, nor it projects test data. To do so, please use `kPCA()`.
#' @details This function may be not valid for all kernels. Do not use it with
#' the RBF, Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels unless
#' you know exactly what you are doing.
#' @param DATA A matrix or data.frame (NOT a kernel matrix) containing the data in
#' feature space. Please note that nrow(DATA) should be higher than ncol(DATA).
#' If the Linear kernel is used, this feature space is simply the original space.
#' @param center	A logical value. If TRUE, the variables are zero-centered. (Defaults: TRUE).
#' @param projected (optional) If desired, the PCA projection (generated, for example, by `kPCA()`)
#' can be included. If DATA is big (especially in the number of rows) this may save
#' some computation time.
#' @param secure (optional) If TRUE, it tests the quality of the loadings
#' This may be slow. (Defaults: FALSE).
#' @return A list with three objects:
#'
#' * The PCA projection (class "matrix") using only the relevant Principal Components.
#'
#' * The loadings.
#'
#' * The values used to center each variable in Feature Space.
#' @export
#' @import ggplot2
#' @examples
#' dat <- matrix(rnorm(150),ncol=30,nrow=50)
#' contributions <- kPCA_imp(dat)
#' contributions$loadings[c("PC1","PC2"),1:5]

kPCA_imp <- function(DATA, center=TRUE, projected=NULL,secure=FALSE) {

  message("Do not use this function if the PCA was created with the RBF,
          Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels")

  n <- nrow(DATA)
  d <- ncol(DATA)

  if(n < d) stop("nrow should be greater than ncol")

  # Centrar les dades
  if(center) {
    X <- scale(DATA,center = TRUE,scale=FALSE)
  } else {
    X <- DATA
  }

  # Importàncies
  spectral_X <-  svd(X)
  loading <- spectral_X$v
  rownames(loading) <- colnames(DATA) ##equivalent a rownames(loadings( pca))
  colnames(loading) <- paste0("PC",1:d) ##equivalent a colnames(loadings( pca))
  loading <- t(loading)

  #kPCA
  if(is.null(projected))   {
    projected <- kPCA(K=Linear(X,cos.norm = FALSE),center=center)
  } else {
    if(ncol(projected)<ncol(X)) message("Ignoring projected data (ncol of projected < ncol of data)")
  }

  # Comprovació de signes
  relneg <- as.data.frame(projected[,1:d] / spectral_X$u %*% diag(spectral_X$d))

  ## Alguns PC tenen el signe canviat (estan a l'inrevés). Així que canviam es signe dels loadings:
  relneg <- which(round(colSums(relneg)) ==-n)
  loading[as.numeric(relneg),] <- -loading[as.numeric(relneg),]

  if(secure) { ## Mirem que la reconstrucció sigui correcta. Aquest pas pot ser lent
    oldvars <- matrix(0,nrow=nrow(DATA),ncol=ncol(loading))
    rownames(oldvars) <- rownames(DATA)
    colnames(oldvars) <- colnames(loading)
    for(i in 1:n)   oldvars[i,] <- solve(loading,t(projected)[1:d,i])
    oldvars <- oldvars + matrix(attr(X,which = "scaled:center"),nrow=n,ncol=ncol(oldvars),byrow = TRUE)
    if(mean(matrix(abs(DATA-oldvars)))>1e-7) warning("Reconstruction may be inacurate")
  }
  return(list(rotated=projected[,1:d], loadings=loading, center=attr(X,which = "scaled:center")))
}


#' Plot the original variables' contribution to a PCA plot
#'
#' @description
#'`kPCA_arrows()` draws arrows on a (kernel) PCA plot to represent the contribution
#' of the original variables to the two displayed Principal Components (PCs).
#'
#' @details It is important to note that the arrows are scaled to match the samples' projection
#' plot. Thus, arrows' directions are correct, but do not expect that their magnitudes
#' match the output of `kPCA_imp()` or other functions(`prcomp`, `princomp...`).
#' (Nevertheless, they should at least be proportional to the real magnitudes.)
#' @param plot A kernel PCA plot generated by `kPCA()`.
#' @param contributions The variables contributions, for instance obtained via `kPCA_imp()`.
#' It is not mandatory to draw all the original variables; a subset of interest
#' can be passed on to this argument.
#' @param colour Color of arrows and labels. (Defaults: "steelblue").
#' @param size Size of the labels. (Defaults: 4).
#' @param ... Additional parameters passed on to geom_segments() and geom_text().
#' @return The PCA plot with the arrows (`ggplot2` object).
#' @export
#' @import ggplot2
#' @examples
#' dat <- matrix(rnorm(500),ncol=10,nrow=50)
#' K <- Linear(dat)
#'
#' ## Computing the kernel PCA. The plot represents PC1 and PC2:
#' kpca <- kPCA(K,plot=1:2)
#'
#' ## Computing the contributions to all the PCS:
#' pcs <- kPCA_imp(dat,secure=FALSE)
#'
#' ## We will draw the arrows for PC1 and PC2.
#' contributions <- t(pcs$loadings[1:2,])
#' rownames(contributions) <- 1:10
#' kPCA_arrows(plot=kpca$plot,contributions=contributions)

kPCA_arrows <- function(plot, contributions, colour="steelblue", size=4, ...) {

  pcs <- plot$data[,as.numeric(sub("PC","",colnames(contributions)))]
  ratio <- max(apply(contributions,2,range)/apply(pcs,2,range))
  X <- Y <- label <- NULL
  colnames(contributions)[1:2] <- c("X","Y")
  if(!("label" %in% colnames(contributions))) {
    contributions  <- data.frame(contributions/ratio,label=rownames(contributions))
  }
  plot <- plot  + ggplot2::geom_segment(data = contributions,  aes(x=0, xend=X, y=0, yend=Y),
                                        arrow = arrow(length = unit(0.3, "cm")),colour=colour)

  # Add text
  plot <- plot + ggplot2::geom_text(data=contributions, aes( x=X, y=Y, label=label),                 ,
                                    size=size,hjust=0.5, vjust=2,colour=colour )
  return(plot)
}


#' Procrustes Analysis
#'
#' Procrustes Analysis compares two PCA/PCoA/MDS/other ordination methods'
#' projections after "removing" the translation, scaling and rotation effects.
#' Thus, they are compared in their configuration of "maximum similarity".
#' Samples in the two projections should be related. The similarity of the projections
#' X1 and X2 is quantified using a correlation-like statistic derived from the
#' symmetric Procrustes sum of squared differences between X1 and X2.
#' @details `Procrustes()` performs a Procrustes Analysis equivalent to
#' `vegan::procrustes(X,Y,scale=FALSE,symmetric=TRUE)`.
#'
#' @param X1 A matrix or data.frame containing a PCA/PCoA/MDS projection.
#' @param X2 A second matrix or data.frame containing a different PCA/PCoA/MDS projection,
#' with the same number of rows than X1.
#' @param labels (optional) A vector of the same length than nrow(X1), or instead,
#' nrow(X1)+nrow(X2). A name will be displayed next to each point.
#' @inheritParams kPCA
#' @return A list containing:
#'
#' * X1 (zero-centered and scaled).
#'
#' * X2 superimposed over X1 (after translating, scaling and rotating X2).
#'
#' * Procrustes correlation between X1 and X2.
#'
#' * (optional) A `ggplot2` plot.
#'
#' @export
#' @import ggplot2
#' @importFrom methods hasArg
#' @examples
#' data1 <- matrix(rnorm(900),ncol=30,nrow=30)
#' data2 <- matrix(rnorm(900),ncol=30,nrow=30)
#' pca1 <- kPCA(Linear(data1),center=TRUE)
#' pca2 <- kPCA(Linear(data2),center=TRUE)
#' procr <- Procrustes(pca1,pca2)
#' # Procrustean correlation between pca1 and pca2:
#' procr$pro.cor
#' # With plot (first two axes):
#' procr <- Procrustes(pca1,pca2,plot=1:2,labels=1:30)
#' procr$plot

Procrustes <- function(X1,X2,plot=NULL,labels=NULL) { ## simètric
  if (nrow(X1) != nrow(X2))
    stop(gettextf("matrices have different number of rows: %d and %d",
                  nrow(X1), nrow(X2)))
  if (ncol(X1) < ncol(X2)) {
    warning("X1 has fewer axes than X2: X1 adjusted to comform X2\n")
    addcols <- ncol(X2) - ncol(X1)
    for (i in 1:addcols) X1 <- cbind(X1, 0)
  } else if(ncol(X1) > ncol(X2)) {
    warning("X2 has fewer axes than X1\n")
  }
  colnms <- colnames(X1)

  ## Center
  X1 <- scale(X1,scale=FALSE)
  X2 <- scale(X2,scale=FALSE)

  ## Scale
  X1 <- frobNorm(as.matrix(X1))
  X2 <- frobNorm(as.matrix(X2))

  # Rotate
  desc <- svd(crossprod(X1,X2))
  Q <- tcrossprod(desc$v, desc$u)
  X2_rot <- X2%*%Q
  colnames(X2_rot) <- colnms
  rownames(X2_rot) <- rownames(X1)

  if(hasArg(plot))   {
    i <- plot[1]
    j <- plot[2]
    colors <- factor(paste("Table",c(rep(1,nrow(X1)),rep(2,nrow(X2)))))
    X <-  data.frame(rbind(X1,X2_rot))
    q <-  ggplot(X, aes(X[,i], X[,j])) +  theme_bw() +
      ggtitle("Procrustes") + xlab(paste0("Dimension", i)) +
      ylab(paste0("Dimension", j)) + coord_fixed() + ggplot2::geom_point(ggplot2::aes(colour = colors))
    if(hasArg(labels)) {
      if(length(labels)==nrow(X1)) {
        labels <- c(labels,labels)
      } else if(length(labels)!=(2*nrow(X1))) {
        stop("label length does not match the dimension of X,Y")
      }
      q <- q +   geom_text(ggplot2::aes(label=labels,hjust=0, vjust=0))
      }
    return(list(X1=X1, X2_rot=X2_rot, pro.cor= sqrt(1-tr( crossprod( X1-(X2_rot),X1-(X2_rot) ) )) ,plot=q))
  } else {
    return(list(X1=X1, X2_rot=X2_rot, pro.cor= sqrt(1-(tr(crossprod( X1-(X2_rot),X1-(X2_rot)))) )))
  }
}


## Helpers

# Ellipse helper kPCA()
#' @keywords internal
#' @noRd
ellipse_helper <-  function(X,ellipse) {
  if (nrow(X) <= 2)       return(NULL)
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  sigma <- var(cbind(X$xvar, X$yvar))
  mu <- c(mean(X$xvar), mean(X$yvar))
  ed <- sqrt(stats::qchisq(ellipse, df = 2))
  data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                   mu, FUN = "+"), groups =X$groups[1])
}
