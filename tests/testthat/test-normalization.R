X <- matrix( c(0.3,-0.1,1.1,1.8,-0.5,-0.2,-0.2,0.1,1.5,3.0,-0.3,-2.3,-0.8,-0.9, 1.1),ncol=5,nrow=3)
K <- Linear(X)



# Cosine normalization (kernel matrix)
test_that("cosNorm computes the right value", {
  K0 <- matrix(c(rep(c(1,0,0,0),2),1),nrow=3,ncol=3)
  expect_equal(cosNorm(K0),K0)
  cosK <- matrix(c(1,-0.3,-0.7,-0.3,1,0,-0.7,0,1), nrow=3,ncol=3)
  expect_equal(round(cosNorm(K),digits=1),cosK)
  expect_equal(round(cosNorm(K),digits=1),round(Linear(cosnormX(X)),digits=1))

})
#> Test passed

# Centering a kernel matrix
test_that("cosNorm computes the right value", {
  centerK <- matrix(c(11.1,-1.7,-9.3,-1.7,1.8,-0.1,-9.3,-0.1,9.4), nrow=3,ncol=3)
  expect_equal(round(centerK(K),digits=1),centerK)
  centeredX <- scale(X,center=TRUE,scale=FALSE)
  expect_equal(centerK,round(Linear(centeredX),digits=1))
})
#> Test passed


# Frobenius norm
test_that("frobNorm computes the right value", {
  frobK <- matrix(c(0.6,-0.1,-0.4, -0.1, 0.1, 0.0, -0.4, 0.0, 0.5),nrow=3,ncol=3)
  expect_equal(round(frobNorm(K),digits=1),frobK )
})
#> Test passed


# Desparsify
test_that("Desparsify computes the right value", {
  N <- 5; D <- 9
  dat <- matrix(rnorm(45),ncol=D,nrow=N)
  row0 <- c(2,5); col0 <- c(1,6,8)
  dat[row0,] <- 0
  dat[,col0] <- 0
  expect_equal(dim(desparsify(dat)), c(N,D-3))
  expect_equal(dim(desparsify(dat,dim=1)), c(N-2,D))
  dat2 <- desparsify(dat,dim=c(1,2))
  expect_equal(dim(dat2), c(N-2,D-3))
  expect_equal(dat2, dat[-row0,-col0])
})


# TSS norm
test_that("TSS computes the right value", {
  expect_equal( round(rowSums( TSS(X))), rep(1,3))
})
#> Test passed


# Euclidean norm (any matrix)
test_that("Euclidean norm  computes the right value", {
  expect_equal(apply(cosnormX(X) ,1,function(X)norm(X,type="2")),rep(1,3))
})
#> Test passed


# Centering a squared matrix by row or column
test_that("centerX computes the right value", {
  Xcrop <- X[,1:3]
  Xrows <- t(scale(t(Xcrop),center=TRUE,scale=FALSE))
  Xcols <- scale(Xcrop,center=TRUE,scale=FALSE)
  attr(Xrows,"scaled:center") <- attr(Xcols,"scaled:center") <- NULL
  expect_equal( centerX(Xcrop,rows=TRUE) , Xrows)
  expect_equal( centerX(Xcrop,rows=FALSE) , Xcols)
})
#> Test passed


# Minmax normalization
test_that("minmax throws errors", {
  valX <- list(min=c(-0.1 ,-0.5, -0.2, -2.3 ),max=apply(X,2,max))
  expect_error(minmax(X,values=valX) , "min and max values should have the same length")
  valX$max <- c( 1.1, 1.8, 1.5, 3.0)
  expect_error(minmax(X,values=valX) , "min and max values should have length")
})
#> Test passed


test_that("minmax computes the right value", {
  minmaxX <- minmax(X,rows=TRUE)
  expect_equal( apply(minmaxX,1,min), rep(0,length=nrow(X)))
  expect_equal( apply(minmaxX,1,max), rep(1,length=nrow(X)))

  minmaxX <- minmax(X)
  expect_equal( apply(minmaxX,2,min), rep(0,length=ncol(X)))
  expect_equal( apply(minmaxX,2,max), rep(1,length=ncol(X)))

  valX <- list(min=apply(X,2,min),max=apply(X,2,max))
  expect_equal(  minmaxX,minmax(X,values=valX))

})
#> Test passed


# Levels per factor variable
test_that("Correct levels per factor variable", {
  levels_list <- list()
  levels_list$Favorite.color <- c("black", "blue", "green", "grey", "lightblue",
                                  "orange", "purple","red","white", "yellow"  )
  levels_list$Liked.new.show <- c("No" , "Yes")
  expect_equal(dummy_var(showdata$Favorite.color)$X,levels_list$Favorite.color)
  expect_equal(dummy_var(showdata[,c(1,5)]),levels_list)

})
#> Test passed


# One-hot-encoding (dummy variables)
test_that("One-hot-encoding is correct", {
  CO2dummies <- dummy_data(CO2[,1:3])
  CO2dummies2 <- dummy_data(CO2[,1:3],lev=dummy_var(CO2[,1:3]))
  expect_equal( CO2dummies, CO2dummies2)
  expect_equal( colnames(CO2dummies),
                c("Plant_Qn1"   ,"Plant_Qn2" , "Plant_Qn3"  ,"Plant_Qc1" ,"Plant_Qc3" ,"Plant_Qc2",
                  "Plant_Mn3"   ,"Plant_Mn2" , "Plant_Mn1"  ,"Plant_Mc2" ,"Plant_Mc3" ,"Plant_Mc1",
                  "Type_Quebec" ,"Type_Mississippi", "Treatment_nonchilled","Treatment_chilled"   ))
  colsums <- colSums(CO2dummies)
  names(colsums) <-NULL
  expect_equal( colsums,  c(rep(7,12),rep(42,4)))
  expect_equal( rowSums(CO2dummies), rep(3,84),ignore_attr=TRUE)
})
#> Test passed


