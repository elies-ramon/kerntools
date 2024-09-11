X <- matrix( c(1,5,2,4,2,8,1,4,-2,6,-3,4 ),ncol=4,nrow=3,byrow = TRUE)

# Linear kernel
test_that("Linear kernel works", {
  K <- Linear(X)
  Kmanual <- matrix( c(46,60,38,60,85,57,38,57,65),ncol=3,nrow=3) ## calculat a mà
  expect_equal(K,Kmanual)

  ## amb cos-normalització
  K <- Linear(X,cos.norm = TRUE)
  Kmanual <- matrix( c(1.00, 0.96, 0.69, 0.96, 1.00, 0.77, 0.69, 0.77,1.00),
                     ncol=3,nrow=3)
  expect_equal(round(K,digits=2),Kmanual)

  ## amb coeficients (a mà)
  Xcoef <-   matrix( c(0.6708204,2.738613, 0.7745967,1.264911,
                       1.3416408,4.381780, 0.3872983,1.264911,
                       -1.3416408, 3.286335, -1.1618950, 1.264911),ncol=4,nrow=3,byrow = TRUE)

  K <- Linear(X,coeff  = c(.45,.3,.15,.1))
  Kmanual <- Linear(Xcoef)
  expect_equal(round(K,digits=4),round(Kmanual,digits=4))
})

test_that("Linear kernel throws errors", {
  expect_error( Linear(X,coeff  = c(.45,.3,.15,.1,0.2)),"should be equal to")
})


# Gaussian RBF kernel
test_that("Gaussian RBF kernel works", {
  D <- RBF(X)
  D2 <- as.matrix(stats::dist(X,method="euclidean"))
  dimnames(D2) <- NULL
  expect_equal(sqrt(D),D2)
  expect_equal(D,RBF(X,g=0))

  K <- RBF(X,g=0.05)
  expect_equal(K,exp(-0.05*D))

  Kmanual <- matrix(c(1.0000000,0.5769498,0.1737739, 0.5769498,1.0000000,0.1652989,
  0.1737739,0.1652989,1.0000000),nrow=3,ncol=3)
  expect_equal(round(K,digits=4),round(Kmanual,digits=4))
})


# Laplacian kernel
test_that("Laplace kernel works", {
  D <- Laplace(X)
  D2 <- matrix(c(0,5,9, 5,0,10,9,10,0),nrow=3,ncol=3)
  colnames(D2) <- rownames(D2) <- 1:3
  expect_equal(D,D2)
  expect_equal(D,Laplace(X,g=0))

  K <- Laplace(X,g=0.1)
  expect_equal(K,exp(-0.1*D))

  Kmanual <- matrix(c( 1.0000000,0.6065307,0.4065697 ,
                       0.6065307,1.0000000,0.3678794 ,
                       0.4065697,0.3678794,1.0000000),nrow=3,ncol=3)
  # colnames(Kmanual) <- rownames(Kmanual) <- 1:3
  expect_equal(round(K,digits=4),round(Kmanual,digits=4),ignore_attr = TRUE)
})


# Frobenius kernel
test_that("Frobenius kernel throws errors", {
  DATA <- list(X1=X,X2=matrix( c(1,5,2,4,2,8,1,4),ncol=4,nrow=2),
               X3=matrix( c(0,4,2,-1,5,-1,-2,1,7,7,8,3 ),ncol=4,nrow=3))
  expect_error(Frobenius(DATA),"should have the same dimensions")
})


test_that("Frobenius kernel works", {
  DATA <- list(X1=X,X2=matrix( c(1,5,2,4,2,8,1,4,-2,6,-3,4 ),ncol=4,nrow=3),
               X3=matrix( c(0,4,2,-1,5,-1,-2,1,7,7,8,3 ),ncol=4,nrow=3),X4=X)
  frob <- Frobenius(DATA,feat_space=TRUE)
  Kmanual <- matrix(c(196,131, 81,196, 131,196, 40,131,
               81, 40,223, 81,  196,131, 81,196),ncol=4,nrow=4)
  fs_manual <- matrix(c( 1, 2 ,-2,  5,  8, 6, 2, 1,-3 , 4 ,  4 ,4,
                  1, 5 , 2,  4,  2, 8, 1, 4,-2 , 6 , -3 ,4,
                  0, 4 , 2, -1,  5,-1,-2, 1, 7 , 7 ,  8 ,3, 1, 2 ,-2,  5,  8, 6,
                  2, 1,-3 , 4 ,  4 ,4),ncol=12,nrow=4,byrow=TRUE)
  rownames(fs_manual) <- rownames(Kmanual) <-  colnames(Kmanual) <- c("X1","X2","X3","X4")
  expect_equal(frob$K,Kmanual)
  expect_equal(frob$feat_space,fs_manual)
  Knorm <- Frobenius(DATA,cos.norm = TRUE)
  expect_equal( Knorm[1,1] , 1)
  expect_equal( Knorm[1,1] , Knorm[1,4] )
  expect_equal(cosNorm(frob$K),Knorm)
})


# Kernels for count data: ruzicka, bray-curtis

test_that("Kernels for count data work", {
  Xcount <- matrix(c(3,2,1, 1, 1, 4, 7, 3, 2 , 5,
                     4,1,5, 3, 4, 2, 3, 6, 1 , 2,
                     2,5,2, 2, 2, 4, 1, 2, 6 , 1,
                     2,1,1, 4, 3, 2, 1, 4, 4 , 5,
                     3,4,6, 1, 2, 2, 4, 0, 2 , 1), nrow=5,byrow = TRUE)

  # Computed with: 1-vegan::vegdist(X,method="jaccard",diag=TRUE,upper=TRUE)
  Ruzicka_vegan <- matrix(c(1.0000000,0.4285714,0.4358974,0.5135135,0.4594595,
                            0.4285714,1.0000000,0.3809524,0.5263158,0.5135135,
                            0.4358974,0.3809524,1.0000000,0.5000000,0.4857143,
                            0.5135135,0.5263158,0.5000000,1.0000000,0.3333333,
                            0.4594595,0.5135135,0.4857143,0.3333333,1.0000000),nrow=5)

  # Computed with: 1-vegan::vegdist(X,method="bray",diag=TRUE,upper=TRUE)
  BC_vegan <- matrix(c(1.0000000,0.6000000,0.6071429,0.6785714,0.6296296,
                       0.6000000,1.0000000,0.5517241,0.6896552,0.6785714,
                       0.6071429,0.5517241,1.0000000,0.6666667,0.6538462,
                       0.6785714,0.6896552,0.6666667,1.0000000,0.5000000,
                       0.6296296,0.6785714,0.6538462,0.5000000,1.0000000),nrow=5)

  expect_equal(round(BrayCurtis(Xcount),digits=4),round(BC_vegan,digits=4),ignore_attr=TRUE)
  expect_equal(round(Ruzicka(Xcount),digits=4),round(Ruzicka_vegan,digits=4),ignore_attr=TRUE)

})

test_that("Kernels for count data errors", {
  X <- matrix(sample(10),nrow=1)
  expect_error(BrayCurtis(X),"X should be a matrix or data.frame with at least two rows")
})


# Dirac
test_that("Dirac kernel works", {
  K1 <- Dirac(showdata[1:10,5]) ## 1 variable
  K1_ma <- matrix(c(1,0,1,0,1,0,0,1,1,0,
                    0,1,0,1,0,1,1,0,0,1,
                    1,0,1,0,1,0,0,1,1,0,
                    0,1,0,1,0,1,1,0,0,1,
                    1,0,1,0,1,0,0,1,1,0,
                    0,1,0,1,0,1,1,0,0,1,
                    0,1,0,1,0,1,1,0,0,1,
                    1,0,1,0,1,0,0,1,1,0,
                    1,0,1,0,1,0,0,1,1,0,
                    0,1,0,1,0,1,1,0,0,1),nrow=10,ncol=10)
  expect_equal(K1, K1_ma,ignore_attr=TRUE)
  expect_equal(K1, Dirac(showdata[1:10,5,drop=FALSE],comp="sum"),ignore_attr=TRUE)

  K1_sum <- Dirac(showdata,feat_space = TRUE,comp="sum")
  K1_mean <- Dirac(showdata,feat_space = TRUE,comp="mean")
  expect_equal(K1_sum$K, Linear(K1_sum$feat_space),ignore_attr=TRUE)
  expect_equal(K1_mean$K, Linear(K1_mean$feat_space),ignore_attr=TRUE)
  expect_equal(K1_mean$K, cosNorm(K1_sum$K),ignore_attr=TRUE)

  absw <- c(1,0.5,0.5,1,1)
  relw <- c(1,0.5,0.5,1,1)/sum(c(1,0.5,0.5,1,1))
  sqrtw <- sqrt(relw)
  K1_w <- Dirac(showdata,feat_space = TRUE,comp="weighted",coeff=absw)
  halfvalue_idx <- grep("Favorite.act",colnames(K1_w$feat_space))
  val1 <- unique(as.vector(K1_w$feat_space[,-halfvalue_idx]))
  val2 <- unique(as.vector(K1_w$feat_space[,halfvalue_idx]))
  expect_equal(val1,c(0,0.5))
  expect_equal(val2,c(0,sqrt(0.125)))
  expect_equal(K1_w$K, Linear(K1_w$feat_space),ignore_attr=TRUE)
  expect_equal(K1_mean$K,  Dirac(showdata,comp="weighted",coeff=rep(1,5)))
})

test_that("Dirac kernel throws errors", {
  expect_error(Dirac(showdata[1,]),"X should be a matrix or data.frame with at least two rows")
  expect_error(Dirac(showdata,comp="weighted"),"Please provide weights")
  expect_error(Dirac(showdata,comp="weighted",coeff=1:3),
               "Number of weights and number of variables do not coincide")
  expect_error(Dirac(showdata,comp="asda"),
               "Option not available")

})


# Intersect, jaccard
absw <- c(1,0.5,0.5,1,1)
setsdata <- matrix(c( "co" ,"mz" ,"ey" ,"nqw","su" ,
                      "ao" ,"kps","hky","hm" ,"eg" ,
                      "asz","ag" ,"mq" ,"hqr","di" ,
                      "be" ,"bnu","qsy","ilw","ch" ,
                      "bf" ,"klu","dkl","em" ,"fq" ,
                      "jt" ,"gy" ,"ery","ghj","fnt",
                      "hrw","ls" ,"es" ,"bjv","dkp",
                      "jm" ,"dv" ,"cg" ,"gv" ,"oq" ,
                      "em" ,"bnw","dmt","vx" ,"aco",
                      "abr","lr" ,"brv","bip","kw" ),nrow=10,byrow = TRUE)
Ksum_ma <- matrix(c( 11, 2, 1, 2, 0, 2, 1, 0, 0, 0,
                     2,12, 2, 1, 3, 2, 1, 0, 0, 1,
                     1, 2,12, 1, 0, 2, 1, 0, 1, 1,
                     2, 1, 1,13, 2, 1, 1, 0, 4, 2,
                     0, 3, 0, 2,12, 1, 1, 1, 1, 2,
                     2, 2, 2, 1, 1,13, 2, 2, 0, 1,
                     1, 1, 1, 1, 1, 2,13, 1, 1, 4,
                     0, 0, 0, 0, 1, 2, 1,10, 3, 0,
                     0, 0, 1, 4, 1, 0, 1, 3,13, 0,
                     0, 1, 1, 2, 2, 1, 4, 0, 0,13),nrow=10)

test_that("Intersect kernel works", {

  Ksum <- Intersect(setsdata,elements=letters,comp="sum",feat_space = TRUE)
  expect_equal(Ksum$K, Ksum_ma,ignore_attr=TRUE)

  Ksum_fs <- cbind(Ksum$feat_space[,,1],Ksum$feat_space[,,2],Ksum$feat_space[,,3],
                   Ksum$feat_space[,,4],Ksum$feat_space[,,5])
  expect_equal(Ksum$K,Linear(Ksum_fs),ignore_attr=TRUE)

  Kmean <- Intersect(setsdata,elements=letters,feat_space = TRUE)
  Kmean_fs <- array(0,dim=c(10,10))
  for(i in 1:5)  Kmean_fs <- Linear(Kmean$feat_space[,,i]) + Kmean_fs
  expect_equal(Kmean$K, Kmean_fs,ignore_attr=TRUE)

  expect_equal(Kmean$K,  Intersect(setsdata,elements=letters,comp="weighted",coeff=rep(1,5)))

  relw <- c(1,0.5,0.5,1,1)/sum(c(1,0.5,0.5,1,1))
  sqrtw <- sqrt(relw)
  Kw <- Intersect(setsdata,elements=letters,feat_space = TRUE,comp="weighted",coeff=absw)

  Kw_fs <- array(0,dim=c(10,10))
  for(i in 1:5)  Kw_fs <- Linear(Kw$feat_space[,,i]) + Kw_fs
  expect_equal(Kw$K, Kw_fs,ignore_attr=TRUE)
  expect_equal(Kmean$feat_space[,,1]/sqrt(1/5),   Kw$feat_space[,,1] /sqrtw[1])
  expect_equal(Kmean$feat_space[,,2]/sqrt(1/5),   Kw$feat_space[,,2] /sqrtw[2])
})


test_that("Jaccard kernel works", {

  Ksum <- Jaccard(setsdata,elements=letters,comp="sum")
  Kmean <- Jaccard(setsdata,elements=letters)
  Kw1 <- Jaccard(setsdata,elements=letters,comp="weighted",coeff=rep(1,5))
  Kw <- Jaccard(setsdata,elements=letters,comp="weighted",coeff=absw)

  expect_equal(rep(5,10), diag(Ksum),ignore_attr=TRUE)
  expect_equal(rep(1,10), diag(Kmean),ignore_attr=TRUE)
  expect_equal(Kw1, Kmean)

  Kw_manual <- matrix(c(1.0000,0.1146,0.0500,0.0813,0.0000,0.0833,0.0417,0.0000, 0.0000,0.0000,
  0.1146,1.0000,0.1250,0.0250,0.1333,0.0875,0.0312,0.0000, 0.0000,0.0625,
  0.0500,0.1250,1.0000,0.0312,0.0000,0.0917,0.0625,0.0000, 0.0312,0.0500,
  0.0813,0.0250,0.0312,1.0000,0.1083,0.0250,0.0312,0.0000, 0.2083,0.1125,
  0.0000,0.1333,0.0000,0.1083,1.0000,0.0625,0.0312,0.0833, 0.0250,0.0938,
  0.0833,0.0875,0.0917,0.0250,0.0625,1.0000,0.0813,0.1458, 0.0000,0.0250,
  0.0417,0.0312,0.0625,0.0312,0.0312,0.0813,1.0000,0.0625, 0.0625,0.2042,
  0.0000,0.0000,0.0000,0.0000,0.0833,0.1458,0.0625,1.0000, 0.2292,0.0000,
  0.0000,0.0000,0.0312,0.2083,0.0250,0.0000,0.0625,0.2292, 1.0000,0.0000,
  0.0000,0.0625,0.0500,0.1125,0.0938,0.0250,0.2042,0.0000, 0.0000,1.0000),nrow=10)

  expect_equal(round(Kw,digits=4), Kw_manual,ignore_attr=TRUE)

})



# Spectrum kernel
letters_ <- c(letters,"_")
strings <- c("hello_world","hello_word","hola_mon","kaixo_mundua",
             "bonjour_le_monde")
names(strings) <- c("english1","english_typo","catalan","basque","french")

test_that("Spectrum kernel works", {
  K_1 <- Spectrum(strings,alphabet=letters_,l=1,feat_space = TRUE)
  K_1$feat_space <- desparsify(K_1$feat_space)

  fs_1ma <- matrix(c(0,0,1,1,1,0,0,0,3,0,0,2,1,0,1,0, 1,
                             0,0,1,1,1,0,0,0,2,0,0,2,1,0,1,0, 1,
                             1,0,0,0,1,0,0,0,1,1,1,2,0,0,0,0, 1,
                             2,0,1,0,0,1,0,1,0,1,1,1,0,2,0,1, 1,
                             0,1,1,2,0,0,1,0,1,1,2,3,1,1,0,0, 2),byrow=TRUE,
                    nrow=5,ncol=17)

  colnames(fs_1ma) <- c("a","b","d","e","h","i","j","k","l","m","n","o","r","u",
                        "w","x","_")
  rownames(fs_1ma) <- names(strings)

  K_1ma <- matrix(c( 19,16 ,  9 ,  4, 15, 16,14 ,  8 ,  4, 14,
                     9, 8 , 10 ,  7, 12, 4, 4 ,  7 , 16, 11,
                     15,14 , 12 , 11, 28),nrow=5,ncol=5)

  expect_equal(rowSums(K_1$feat_space),nchar(strings))
  expect_equal(K_1$K,K_1ma,ignore_attr = TRUE)
  expect_equal(K_1$feat_space, fs_1ma)
  expect_equal(K_1$K,Linear(K_1$feat_space),ignore_attr = TRUE)

  K_2cos <- Spectrum(strings,alphabet=letters_,l=2,cos.norm = TRUE,feat_space=TRUE)
  expect_equal(Linear(K_2cos$feat_space),K_2cos$K)

  K_2 <- Spectrum(strings,alphabet=letters_,l=2,feat_space=FALSE)
  expect_equal(cosNorm(K_2),K_2cos$K)

  expect_contains(class(K_2),"matrix")
  expect_type(K_2cos, "list")

  contains_pattern <-  K_2cos$feat_space[,colSums(K_2cos$feat_space)>0]!=0
  expect_identical(contains_pattern["french","mo"],contains_pattern["catalan","mo"])
  expect_identical(contains_pattern["french","nd"],contains_pattern["basque","nd"])
  expect_identical(contains_pattern["english1","o_"],contains_pattern["basque","o_"])

  ### weights / group.ids
  K_1w <- Spectrum(strings,alphabet=letters_,l=1,group.ids = c(1,1,2,3,4),
                   weights = c(0.5,0.5,1,1,1),feat_space = TRUE)
  K_1w$feat_space <- desparsify(K_1w$feat_space)
  expect_equal( K_1w$feat_space[2:4,],K_1$feat_space[3:5,],ignore_attr = TRUE)
  expect_equal( K_1w$feat_space[1,],colMeans(K_1$feat_space[1:2,]),ignore_attr = TRUE)
})

test_that("Spectrum kernel throws errors", {
  expect_error(Spectrum(strings,alphabet=letters_,l=1,group.ids = c(1,1,2),
                   weights = c(0.5,0.5,1,1,1)),"Ids length should be the same than x")
  expect_error(Spectrum(strings,alphabet=letters_,l=1,weights = c(0.5,0.5,2:5)),
               "weights length should be the same than x")

})


# Kendall's tau
color_list <-  c("black","blue","green","grey","lightblue","orange","purple",
                 "red","white","yellow")
survey1 <- 1:10
survey2 <- 10:1
survey3 <- c(10,3,4,7 , 8 , 1,  6,  2,  5 , 9)
color <- cbind(survey1,survey2,survey3) # samples is columns
rownames(color) <- color_list
food <- matrix(c(10, 1,18, 25,30, 7, 5,20, 5, 12, 7,20, 20, 3,22),ncol=5,nrow=3)
rownames(food) <- colnames(color)
colnames(food) <- c("spinach", "chicken", "beef" , "salad","lentils")

test_that("Kendall's tau kernel works", {

  K1 <- Kendall(color)
  expect_equal(nrow(K1),3)
  expect_equal(K1[1,2],-1)
  expect_equal(round(K1[1,3],digits = 1),0)
  expect_equal( K1[1,3],-K1[2,3])

  expect_equal( nrow(Kendall(food)),5)
  K2 <- Kendall(food,samples.in.rows=TRUE)
  Kmanual <- matrix(c(1.0, 0.2, 0.4, 0.2, 1.0 ,-0.4, 0.4,-0.4, 1.0), nrow=3,ncol=3)
  expect_equal(K2,Kmanual,ignore_attr = TRUE)

  X <- list(color=color,food=t(food)) #All samples in columns
  K <- array(dim=c(3,3,2))
  K[,,1] <- K1
  K[,,2] <- K2

  expect_equal(Kendall(X),MKC(K))
})

test_that("Kendall's tau kernel throws errors", {
  expect_error(Kendall(list(color=color,food= food)),
               "All list's elements should have the same number of columns")
  expect_error(Kendall(list(color=color,food= food),samples.in.rows = TRUE),
               "All list's elements should have the same number of rows")
  expect_error(Kendall(list(color=color,food= t(food)),comp="hola"),
               "Option not available")
})

