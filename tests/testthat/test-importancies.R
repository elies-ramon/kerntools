x <- c(0.42,1.29,-0.81,-0.22,0,0.88,1.83,-0.91,-0.38,1.03,0.21,-0.71,-0.46,-0.49,0.86,-1.77,1.01,0.52,-0.06,-0.49,0.18,-0.73,0.76,0.73,1.92,0.55,-0.12,0.02,-1.39,0.07)
names(x) <- paste0("Feat",1:30)

test_that("plotImp works", {
  pdf(NULL)
  plotimp <- plotImp(x,absolute=TRUE,relative=FALSE)
  expect_equal(plotimp$first_features,names(sort(abs(x),decreasing = TRUE)))
  plotimp2 <- plotImp(x,absolute=TRUE,relative=TRUE)
  expect_equal(plotimp$first_features,plotimp2$first_features)
  plotimp3 <- plotImp(x,absolute=FALSE,relative=TRUE)
  expect_equal(plotimp$first_features,plotimp3$first_features)
  plotimp4 <- plotImp(x,absolute=FALSE,relative=FALSE)
  expect_equal(plotimp$first_features,plotimp4$first_features)

  expect_equal(plotimp2$cumsum,1)
  expect_equal(plotimp3$cumsum,1)

  plotimp2 <- plotImp(x,absolute=TRUE,relative=TRUE,nfeat=10)
  expect_equal(plotimp$first_features[1:10],plotimp2$first_features)
  expect_equal(plotimp2$cumsum, sum(sort(abs(x)/sum(abs(x)),decreasing = TRUE)[1:10]))

  plotimp <- plotImp(x,absolute=TRUE,relative=FALSE,names=1:30)
  expect_equal(plotimp$first_features, sub("Feat","",names(sort(abs(x),decreasing = TRUE))))
})
