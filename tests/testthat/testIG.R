test_that("IG works", {
  hsmodel <- "visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, 6)
  k <- rep(7, 6)
  set.seed(1)
  sample1 <- rIG(10^6,sigma0, s, k, reps=1, typeA="triang")[[1]]
  sample2 <- rIG(10^6,sigma0, s, k, reps=1, typeA="symm")[[1]]

  cov1 <- lavaan::lav_matrix_vech(cov(sample1))
  cov2 <- lavaan::lav_matrix_vech(cov(sample2))
  skew1 <- psych::skew(sample1)
  skew2 <- psych::skew(sample2)
  kurt1 <- unname(psych::kurtosi(sample1))
  kurt2 <- unname(psych::kurtosi(sample2))

  expect_equal(cov1, target0, tolerance=0.01)
  expect_equal(cov2, target0, tolerance=0.01)
  expect_equal(skew1, s, tolerance=0.1)
  expect_equal(skew2, s, tolerance=0.1)
  expect_equal(kurt1, k, tolerance=0.1)
  expect_equal(kurt2, k, tolerance=0.1)

})


