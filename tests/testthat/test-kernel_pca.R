X <- scale(USArrests,center = TRUE,scale = TRUE)
prcomp_pca <- stats::prcomp(X,center = FALSE,scale. = FALSE)
prcomp_proj <- round(prcomp_pca$x,digits=4)

K <- Linear(X,cos.norm = FALSE)
our_pca <-  kPCA(K=K,center = FALSE,plot=1:2)
our_proj <- as.matrix(round(our_pca$projection[,1:4],digits=4))
our_rot <- kPCA_imp(X,center=FALSE,secure=FALSE)

quocient <- apply( our_proj/prcomp_proj,2,unique)
signs <- rep(1,4)
names(signs) <- c("PC1", "PC2", "PC3", "PC4" )

our_pca3 <-  kPCA(K=RBF(X,g=0.05),center = TRUE)


# Test kernel PCA
test_that("Normal (kernel = linear) PCA works", {
  expect_equal( abs(our_proj),abs(prcomp_proj))
  expect_equal(abs(quocient),signs)
})


test_that("Projecting a test over a PCA", {
  N <- nrow(X)
  K2 <- K
  rownames(K2)<-colnames(K2)<-1:N
  our_pca2 <- kPCA(K=K,center = FALSE,Ktest = K2)
  our_proj2 <- as.matrix(round(our_pca2[,1:4],digits=4))
  our_test <- our_proj2[(N+1):nrow(our_proj2),]
  rownames(our_test) <- rownames(X)
  expect_equal(our_proj2[1:N,],our_test)
  expect_equal(our_proj2[1:N,],our_proj)

})


test_that("Kernel PCA works", {
  pc1 <- c(-0.29709337,-0.42643848,-0.44461601, 0.02126754,-0.52059390,-0.37872522, 0.35751141,-0.01414932,-0.63131157 ,
  -0.42227493, 0.23906290, 0.43746703,-0.36549892, 0.15458532, 0.55067872, 0.23950352, 0.19316097,-0.42467210 ,
  0.55907924,-0.47577930, 0.14151164,-0.53579895, 0.45116121,-0.27750203,-0.20154968, 0.32997911, 0.36238738 ,
  -0.57741030, 0.56677343,-0.03433462,-0.51905781,-0.42635637,-0.29604985, 0.60354160, 0.07657383, 0.09976540 ,
  -0.01037478, 0.25833660, 0.20770536,-0.36063240, 0.48711536,-0.29508661,-0.36461582, 0.16211123, 0.53497866 ,
  0.02726138, 0.07318451, 0.46847962, 0.51514787, 0.18159152)
  expect_equal(our_pca3$PC1  ,pc1)
})



# Test contribution to PCs
test_that("kPCA_imp works", {
  prcomp_rot <- prcomp_pca$rotation
  expect_equal(prcomp_rot  ,t(quocient*our_rot$loadings))

  our_rot2 <- kPCA_imp(X, center=TRUE,secure=FALSE)
  expect_lt(max(our_rot2$center),1e-15)

  our_rot3 <- kPCA_imp(X,projected = our_pca$projection,center=FALSE,secure=FALSE)
  expect_identical(our_rot,our_rot3)
})


# Test Kernel arrows
test_that("kernel arrows' plotis different than kPCA's", {
  our_arrows <- kPCA_arrows(plot = our_pca$plot,contributions = t(our_rot$loadings[1:2,]))
  expect_false( identical(our_pca$plot,our_arrows))
})


# Test Procrustes
test_that("Procrustes throws errors", {
  expect_error(Procrustes(our_pca$projection[1:20,],our_pca3), "different number of rows")
  expect_warning(Procrustes(our_pca3,our_pca$projection[,1:2]), "X2 has fewer axes than X1")
  expect_warning(diffcols <- Procrustes(our_pca$projection[,1:2],our_pca3), "X1 has fewer axes than X2")
  expect_equal(ncol(diffcols$X1),ncol(our_pca3))
  expect_equal(sum(diffcols$X1[,-c(1:2)]),0)
})

test_that("Procrustes works", {
  proc1 <- Procrustes(our_pca3,our_pca$projection)
  proc2 <- Procrustes(our_pca$projection,our_pca3)
  expect_equal(proc1$pro.cor,proc2$pro.cor) ## symmetry
  expect_equal(round(proc1$pro.cor,digits=2),round(sqrt(1-0.1849),digits=2) ) #nombre calculat amb vegan::procrustes
  #x1f1 i x2f1 calculats amb vegan::procrustes()
  x1f1 <- c(-0.08,-0.11,-0.12, 0.01,-0.14,-0.10, 0.09, 0.00,-0.17,-0.11, 0.06, 0.11,-0.10, 0.04, 0.14, 0.06,  0.05,-0.11,0.15 ,
  -0.12, 0.04,-0.14, 0.12,-0.07,-0.05, 0.09, 0.09,-0.15, 0.15,-0.01,-0.14,-0.11,-0.08, 0.16, 0.02,  0.03, 0.00,0.07 ,
  0.05,-0.09, 0.13,-0.08,-0.10, 0.04, 0.14, 0.01, 0.02, 0.12, 0.13, 0.05)
  x2f1 <- c(-0.07,-0.14,-0.12, 0.01,-0.17,-0.10, 0.10, 0.00,-0.21,
             -0.12, 0.07, 0.12,-0.10, 0.04, 0.16, 0.06, 0.05,-0.11,
             0.17,-0.13, 0.04,-0.15, 0.12,-0.08,-0.05, 0.08, 0.09,
             -0.20, 0.17,-0.01,-0.14,-0.12,-0.09, 0.21, 0.02, 0.02,
             0.00, 0.06, 0.07,-0.10, 0.14,-0.07,-0.09, 0.04, 0.19,
             0.01, 0.02, 0.14, 0.15, 0.04)
  names(x1f1) <- names(x2f1) <-rownames(our_pca3)
  expect_equal(round(proc1$X1[,1L],digits=2) ,x1f1 )## symmetry
  expect_equal(round(proc1$X2[,1L],digits=2) ,x2f1 )## symmetry

})
