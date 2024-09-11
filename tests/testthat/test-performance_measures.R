# FALTA
# Revisar el resultat numèric dels tests a mà




# NMSE
test_that("NMSE throws errors", {
  expect_error(nmse(letters, LETTERS), "should be numeric")
  expect_error(nmse(1,1), "should have length > 1")
  expect_error(nmse(rep(1,3), 1:3), "should have non-zero variance")
  expect_error(nmse(1:3, 1:2),"should have the same length")
})
#> Test passed

test_that("NMSE computes the right error value", {
  expect_equal(nmse(1:3, 1:3), 0)
  expect_equal(nmse(1:3, c(1,1,3)), 0.5)
})
#> Test passed


# Accuracy
test_that("Accuracy throws errors", {
  expect_error(Acc(data.frame(a=c(1,2),b=c(3,4))), "should be class matrix or table")
  expect_error(Acc(1:4), "should be class matrix or table")
})
#> Test passed

test_that("Accuracy computes the right value", {
  expect_equal(Acc(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2)))),1)
  expect_equal(Acc(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1)))),0.5)
  expect_equal(Acc(table(factor( rep("b",4),levels=c("a","b")), factor( rep("a",4),levels=c("a","b")))),0)
})
#> Test passed


test_that(" weighted Accuracy computes the right value", {
  expect_equal(Acc(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),weighted=TRUE),1)
  expect_equal(Acc(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))),weighted=TRUE),2/3)
  expect_equal(Acc(table(factor( rep("b",4),levels=c("a","b")), factor( rep("a",4),levels=c("a","b"))),weighted=TRUE),NaN)
})
#> Test passed

test_that("weighted Accuracy = Accuracy when the classes in the target are balanced", {
  expect_equal(Acc(table(c(rep("a",2),rep("b",2)), c(rep("a",1),rep("b",3)))),
               Acc(table(c(rep("a",2),rep("b",2)), c(rep("a",1),rep("b",3))),weighted=TRUE))
})
#> Test passed


#  Accuracy of random model
test_that("Acc_rnd throws errors", {
  expect_error(Acc_rnd(1), "should have length > 1")
  expect_error(Acc_rnd(letters,freq=T), "should be numeric")
  expect_error(Acc_rnd(c(0.5,0.1),freq=T), "frequencies should sum to 1")
})
#> Test passed

test_that("Acc_rnd computes the right value", {
  expect_equal(Acc_rnd(c("a","a","a","a"),freq=F),1)
  expect_equal(Acc_rnd(c("a","a","b","b"),freq=F),0.5)
  expect_equal(Acc_rnd(c("a","a","a","b"),freq=F),0.625)
  expect_equal(Acc_rnd(c("a","a","b","b","c","c"),freq=F),1/3) # > 2 classes
  expect_equal(Acc_rnd(c(0.5,0.5),freq=T),0.5) # Freq
  expect_equal(Acc_rnd(c(0.75,0.25),freq=T),0.625) # Freq
  expect_equal(Acc_rnd(c(0.3,0.43,0.27),freq=T),0.3478)  # Freq and > 2 classes
})
#> Test passed


#  CI (Confidence Interval) using Normal Approximation
test_that("Normal_CI throws errors", {
  expect_error(Normal_CI(value=1:4,ntest=40), "should have length==1")
  expect_error(Normal_CI(value=0.4,ntest=40:45), "should have length==1")
  expect_error(Normal_CI(value="A",ntest="B"), "should be numeric")
  expect_error(Normal_CI(value=1.5,ntest=400), "cannot be lower than 0 or greater than 1")
  expect_error(Normal_CI(value=0.5,ntest=3), "should be 30 or greater")
  expect_error(Normal_CI(value=0.5,ntest=400,confidence=130), "should be greater than 0 and lower than 100")
})
#> Test passed

test_that("Normal_CI computes the right CI values", {
  expect_equal(round(Normal_CI(0.8, ntest=50, confidence=95),4),c(0.6891,0.9109))
})
#> Test passed


#  CI (Confidence Interval) using Bootstrapping
test_that("Boots_CI throws errors", {
  expect_error(Boots_CI(target=1:3, pred=1:2), "should have the same length")
  expect_error(Boots_CI(target=letters, pred=letters,nboots = 1), "please set a higher nboots")
  y <- c(rep("a",20),rep("b",30))
  y_pred <-  c(rep("a",30),rep("b",20))
  expect_error(Boots_CI(y, y_pred,nboots = 500,confidence=-10),"should be greater than 0 and lower than 100")
  expect_error(Boots_CI(y, y_pred,nboots = 500,index="asdas",confidence=95),
               "Measure not available")

})
#> Test passed

test_that("Boots_CI computes the right CI values", {
  y <- c(rep("a",20),rep("b",30))
  y_pred <-  c(rep("a",30),rep("b",20))
  result <- round(Boots_CI(target=y,pred=y_pred, index="acc", nboots=1000, confidence=95),1) ## Accuracy
  expect_equal(as.numeric(result),c(0.8,  0.7,  0.9))
  expect_equal(names(result),c("", "2.5%", "97.5%"))

})
#> Test passed


# Helpers
test_that("Prec, Rec, Spe and F1 throw errors", {
  expect_error(imbfunct_helper(data.frame(a=c(1,2),b=c(3,4)),multi.class="none"), "should be class matrix or table")
  expect_error(imbfunct_helper(table( c(rep("a",3),rep("b",4)), rep("a",7)),multi.class="none"), "should be squared")
  expect_error(imbfunct_helper(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))), multi.class="asdas"),"Option not available")
})


# Precision or PPV
test_that("Prec computes the right value", {
  expect_equal(Prec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2)))),1)
  expect_equal(Prec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "micro"),1)
  result <- c(1,1)
  names(result) <- c("a","b")
  expect_equal(Prec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "none"),result)

  expect_equal(Prec(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))),multi.class="macro"),2/3)
  expect_equal(Prec(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="micro"),0.8)
  result[2] <- 0.75
  expect_equal(Prec(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="none"),result)

  ct <- table(c(rep("b",1),rep("a",3),rep("c",3)), c(rep("a",3),rep("b",2),rep("c",2)))
  expect_equal(round(Prec(ct,multi.class = "micro"),4),0.5714)
  expect_equal(round(Prec(ct,multi.class = "macro"),4),0.5556)
  result <- c(2/3 ,0 ,1 )
  names(result) <- c("a","b","c")
  expect_equal(Prec(ct,multi.class = "none"),result)

})
#> Test passed

# Recall or Sensitivity or TPR
test_that("Rec computes the right value", {
  expect_equal(Rec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2)))),1)
  expect_equal(Rec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "micro"),1)
  result <- c(1,1)
  names(result) <- c("a","b")
  expect_equal(Rec(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "none"),result)

  expect_equal(Rec(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))),multi.class="macro"),2/3)
  expect_equal(Rec(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="micro"),0.8)
  result[1] <- 0.5
  expect_equal(Rec(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="none"),result)

  ct <- table(c(rep("b",1),rep("a",3),rep("c",3)), c(rep("a",3),rep("b",2),rep("c",2)))
  expect_equal(round(Rec(ct,multi.class = "micro"),4),0.5714)
  expect_equal(round(Rec(ct,multi.class = "macro"),4),0.4444)
  result <- c(2/3 ,0 ,2/3 )
  names(result) <- c("a","b","c")
  expect_equal(Rec(ct,multi.class = "none"),result)
})
#> Test passed


# Specificity or TNR
test_that("Spe computes the right value", {
  expect_equal(Spe(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2)))),1)
  expect_equal(Spe(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "micro"),1)
  result <- c(1,1)
  names(result) <- c("a","b")
  expect_equal(Spe(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "none"),result)

  expect_equal(Spe(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))),multi.class="macro"),2/3)
  expect_equal(Spe(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="micro"),0.8)
  result[2] <- 0.5
  expect_equal(Spe(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="none"),result)

  ct <- table(c(rep("b",1),rep("a",3),rep("c",3)), c(rep("a",3),rep("b",2),rep("c",2)))
  expect_equal(round(Spe(ct,multi.class = "micro"),4), 0.8571)
  expect_equal(round(Spe(ct,multi.class = "macro"),4),0.8333)
  result <- c(0.75 ,1 ,0.75)
  names(result) <- c("a","b","c")
  expect_equal(Spe(ct,multi.class = "none"),result)
})
#> Test passed


# F1
test_that("F1 computes the right value", {
  expect_equal(F1(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2)))),1)
  expect_equal(F1(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "micro"),1)
  result <- c(1,1)
  names(result) <- c("a","b")
  expect_equal(F1(table(c(rep("a",2),rep("b",2)), c(rep("a",2),rep("b",2))),multi.class = "none"),result)

  expect_equal(F1(table(c(rep("a",1),rep("b",3)), c(rep("a",3),rep("b",1))),multi.class="macro"),0.5)
  expect_equal(F1(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="micro"),0.8)
  result[1] <- 0.6667
  result[2] <-  0.8571
  expect_equal(round(F1(table(c(rep("a",2),rep("b",3)), c(rep("a",1),rep("b",4))),multi.class="none"),4),result)

  ct <- table(c(rep("b",1),rep("a",3),rep("c",3)), c(rep("a",3),rep("b",2),rep("c",2)))
  expect_equal(round(F1(ct,multi.class = "micro"),4),0.5714)
  expect_equal(round(F1(ct,multi.class = "macro"),4),0.4889)
  result <- c(2/3 ,0 ,0.8 )
  names(result) <- c("a","b","c")
  expect_equal(F1(ct,multi.class = "none"),result)
})
#> Test passed

