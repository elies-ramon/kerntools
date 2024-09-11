test_that("MKC works", {
  elements <- c(1,0.3,0.6,0.3,1,0.25,0.6,0.25,1,
                c(1,-0.3,-0.7,-0.3,1,0,-0.7,0,1))
  K <- array(elements,dim=c(3,3,2))

  expect_equal(MKC(K), (K[,,1]+K[,,2])/2)
  expect_equal(MKC(K,coeff = c(1/2,1/2)), (K[,,1]+K[,,2])/2)
  expect_equal(MKC(K,coeff = c(0.1,0.9)), (0.1*K[,,1])+(0.9*K[,,2]))
  expect_equal(MKC(K,coeff = c(5,4)), ((5*K[,,1])+(4*K[,,2]))/9)
})
#> Test passed


test_that("MKC throw errors", {
  K <- array(0,dim=c(3,3,4))
  expect_error(MKC(K,coeff = c(0.1,0.9)), "Length of the coefficients vector different to the number of matrices")
})
#> Test passed

