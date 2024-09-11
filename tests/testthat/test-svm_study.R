
test_that("svm_imp works", {
  X <- matrix( c(0.3,-0.1,1.1,1.8,-0.5,-0.2,-0.2,0.1,1.5,3.0,-0.3,-2.3,-0.8,-0.9, 1.1),ncol=5,nrow=3)
  sv_index <- c(1,3)
  coeff <- c(-0.5,0.5)

  importances <- svm_imp(X = X,coeff = coeff,svindx = sv_index)
  Xidx <- X[c(1,3),]
  Xidx[1,] <- -0.5*Xidx[1,]
  Xidx[2,] <- 0.5*Xidx[2,]
  expect_equal(importances,abs(colSums(Xidx)))

  importances <- svm_imp(X = X,coeff = coeff,svindx = sv_index, center=TRUE,scale=TRUE,
                         result="original")
  Xidx <- scale(X,center=TRUE,scale=TRUE)[c(1,3),]
  Xidx[1,] <- -0.5*Xidx[1,]
  Xidx[2,] <- 0.5*Xidx[2,]
  expect_equal(importances,colSums(Xidx))

  importances <- svm_imp(X = X,coeff = coeff,svindx = sv_index, cos.norm = TRUE,
                         result="squared")
  Xidx <- cosnormX(X[c(1,3),])
  Xidx[1,] <- -0.5*Xidx[1,]
  Xidx[2,] <- 0.5*Xidx[2,]
  expect_equal(importances,colSums(Xidx)^2)
})
