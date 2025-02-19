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

  plotimp <- plotImp(x,absolute=TRUE,relative=FALSE,names=1:30)
  expect_error( plotImp(x,nfeat=50), "nfeat cannot be larger than x")
})


test_that("aggregate_imp works", {
  vect_grep <- Vectorize(grepl)
  importances <- matrix(rnorm(120),nrow=4,ncol=30)
  rownames(importances) <- paste0("sample",1:4)
  colnames(importances) <- paste0("Feat",rep(1:5,times=2*(1:5)), "_", unlist(lapply(2*(1:5),function(x)LETTERS[1:x])))

  aggrimp <- aggregate_imp(X=importances,samples="rows")
  expect_equal(aggrimp[,"Feat3"],rowSums(importances[,vect_grep(sub('_.*', '',colnames(importances)),"Feat3")]))
  expect_equal(aggrimp[,"Feat5"],rowSums(importances[,vect_grep(sub('_.*', '',colnames(importances)),"Feat5")]))
  expect_equal(t(aggrimp), aggregate_imp(X=t(importances),samples="cols"))

  groups <- paste0("Feat",1:5)
  expect_equal(aggrimp, aggregate_imp(X=importances,lev=groups))
  groups <- paste0("Feat",2:4)
  expect_error(aggregate_imp(X=importances,lev=groups),"arguments imply differing number of rows")
  expect_error(aggregate_imp(X=importances,lev=31:35),"arguments imply differing number of rows")

})
