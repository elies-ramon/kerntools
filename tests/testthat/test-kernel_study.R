X <- matrix( c(0.3,-0.1,1.1,1.8,-0.5,-0.2,-0.2,0.1,1.5,3.0,-0.3,-2.3,-0.8,-0.9, 1.1),ncol=5,nrow=3)
K <- RBF(X,g=0.01)


# Kernel matrix as input
test_that("Kernel matrices throw errors", {
  pseudoK <- matrix(1:8,nrow=2,ncol=4)
  expect_error(kprecondition_helper(pseudoK), "should be squared")
  pseudoK <- data.frame(A=c(1,0,0),B=c(0,1,0),C=c(0,0,1))
  expect_error(kprecondition_helper(pseudoK), "should be class = matrix")
  pseudoK <- matrix(LETTERS[-1],nrow=5,ncol=5)
  expect_error(kprecondition_helper(pseudoK), "should contain only numbers")
})


test_that("vonNeumann works", {
  expect_equal(round(vonNeumann(K),digits=2),0.40)
})


test_that("Estimate gamma works", {
  result <- estimate_gamma(X)
  expect_equal( result$d_criterion, 1/5)
  expect_equal( round(result$scale_criterion,digits=2), 0.12)
})

test_that("simK works", {
  Klist <- list(rbf=K,cosnorm = Linear(X,cos.norm = T),linear=Linear(X))
  result <- round(simK(Klist),digits=1)
  similarity <- matrix( c(1,0.3,0.2,0.3,1,0.9,0.2,0.9,1),ncol=3,nrow=3)
  colnames(similarity) <- rownames(similarity) <- c("rbf","cosnorm","linear")
  expect_equal( result,similarity)
})

test_that("simK throws errors", {
  Klist <- list(rbf=K,rbf2 = K[-1,])
  expect_error( simK(Klist),"not squared")
})

test_that("KTA (?) works", {
  K1 <- RBF(iris[1:100,1:4],g=0.1)
  y <- factor(iris[1:100,5])
  expect_equal(round(KTA(K1,y),digits=2),0.41)

  y <-  (c(-4,1,2,1,3,-2,-1,2,-6,3,1,-1,-1,-2,4,-3,1,-1,-1,5))
  K <- Linear(y,cos.norm = T)
  y[y>0] <- 1
  y[y<0] <- -1
  expect_equal( KTA(K=K,y=factor(y)),1)
})

test_that("KTA throws errors", {
  K1 <- RBF(iris[ ,1:4],g=0.1)
  y <- factor(iris[ ,5])
  expect_error(KTA(K1,y), "should have 2 levels")
})
