test_that("PLSIM works", {
  hsmodel <- "visual  =~ x1 + x2
              textual =~ x4 + x5"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, ncol(sigma0))
  k <- rep(7, ncol(sigma0))
  set.seed(1)
  sample1 <- rPLSIM(10^6,sigma0, s, k, reps=1 )[[1]][[1]]

  cov1 <- lavaan::lav_matrix_vech(cov(sample1))
  skew1 <- psych::skew(sample1)
  kurt1 <- unname(psych::kurtosi(sample1))


  expect_equal(cov1, target0, tolerance=0.01)
  expect_equal(skew1, s, tolerance=0.1)
  expect_equal(kurt1, k, tolerance=0.1)

})


