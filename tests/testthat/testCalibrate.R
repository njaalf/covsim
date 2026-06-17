test_that("calibrate_ig + simulate matches target moments", {
  hsmodel <- "visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, 6)
  k <- rep(7, 6)

  set.seed(1)
  cal <- calibrate_ig(sigma0, s, k, typeA = "triang")
  expect_s3_class(cal, "covsim_calib_ig")
  expect_s3_class(cal, "covsim_calib")

  sample1 <- stats::simulate(cal, nsim = 1, N = 10^6)[[1]]
  expect_equal(lavaan::lav_matrix_vech(cov(sample1)), target0, tolerance = 0.01)
  expect_equal(unname(psych::skew(sample1)), s, tolerance = 0.1)
  expect_equal(unname(psych::kurtosi(sample1)), k, tolerance = 0.1)
})


test_that("calibrate_ig calibration can be reused with different N and nsim", {
  set.seed(1)
  sigma0 <- cov(MASS::mvrnorm(20, rep(0, 3), diag(3)))
  cal <- calibrate_ig(sigma0, skewness = rep(1, 3), excesskurtosis = rep(3, 3))

  set.seed(2)
  draws <- stats::simulate(cal, nsim = 4, N = 100)
  expect_length(draws, 4)
  expect_true(all(sapply(draws, nrow) == 100))
  expect_true(all(sapply(draws, ncol) == 3))

  # second call with different shape works
  draws2 <- stats::simulate(cal, nsim = 2, N = 500)
  expect_length(draws2, 2)
  expect_true(all(sapply(draws2, nrow) == 500))
})


test_that("rIG wrapper produces same shape as before", {
  set.seed(1)
  sigma0 <- cov(MASS::mvrnorm(20, rep(0, 3), diag(3)))
  out <- rIG(N = 200, sigma0, skewness = rep(0.5, 3),
             excesskurtosis = rep(2, 3), reps = 3)
  expect_type(out, "list")
  expect_length(out, 3)
  expect_true(all(sapply(out, nrow) == 200))
  expect_true(all(sapply(out, ncol) == 3))
})


test_that("calibrate_plsim + simulate matches target moments", {
  hsmodel <- "visual  =~ x1 + x2
              textual =~ x4 + x5"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, ncol(sigma0))
  k <- rep(7, ncol(sigma0))

  set.seed(1)
  cal <- calibrate_plsim(sigma0, s, k, verbose = FALSE)
  expect_s3_class(cal, "covsim_calib_plsim")

  sample1 <- stats::simulate(cal, nsim = 1, N = 10^6)[[1]]
  expect_equal(lavaan::lav_matrix_vech(cov(sample1)), target0, tolerance = 0.01)
  expect_equal(unname(psych::skew(sample1)), s, tolerance = 0.1)
  expect_equal(unname(psych::kurtosi(sample1)), k, tolerance = 0.1)
})


test_that("rPLSIM wrapper preserves legacy list(samples, model) format", {
  set.seed(1)
  sigma.target <- cov(MASS::mvrnorm(5, rep(0, 3), diag(3)))
  res <- rPLSIM(10^4, sigma.target, skewness = rep(1, 3),
                excesskurtosis = rep(4, 3), verbose = FALSE)
  expect_named(res, c("samples", "model"))
  expect_named(res$model, c("a", "b", "gamma", "z.corr"))
  expect_length(res$samples, 1)
})


test_that("calibrate_vita + simulate produces correct covariance (small Nmax)", {
  skip_on_cran()
  set.seed(1)
  sigma.target <- cov(MASS::mvrnorm(20, mu = rep(0, 3), Sigma = diag(1, 3)))
  marginsnorm <- lapply(X = sqrt(diag(sigma.target)),
                       function(X) list(distr = "norm", sd = X))

  cal <- calibrate_vita(marginsnorm, sigma.target = sigma.target,
                        Nmax = 10^4, cores = 1, verbose = FALSE)
  expect_s3_class(cal, "covsim_calib_vita")

  set.seed(2)
  draws <- stats::simulate(cal, nsim = 2, N = 10^4)
  expect_length(draws, 2)
  # covariance should match within reasonable tolerance for the (low) Nmax used
  cov_est <- cov(draws[[1]])
  expect_equal(as.numeric(cov_est), as.numeric(sigma.target), tolerance = 0.15)
})


test_that("vita calibration bias is small under large-sample averaging", {
  # Isolates calibration bias from sampling noise by averaging the empirical
  # covariance over several large simulated samples. Skipped on CRAN because
  # the large samples make this slow.
  skip_on_cran()

  set.seed(1)
  d <- 5
  sigma.target <- cov(MASS::mvrnorm(50, mu = rep(0, d), Sigma = diag(d)))
  margins <- lapply(sqrt(diag(sigma.target)),
                    function(s) list(distr = "norm", sd = s))

  set.seed(42)
  cal <- calibrate_vita(margins, sigma.target = sigma.target,
                        Nmax = 5e4, cores = 1, verbose = FALSE)

  cov_sum <- matrix(0, d, d)
  reps <- 3
  for (k in seq_len(reps)) {
    set.seed(100 + k)
    samp <- stats::simulate(cal, nsim = 1, N = 5e5, cores = 1)[[1]]
    cov_sum <- cov_sum + cov(samp)
  }
  cov_mean <- cov_sum / reps

  # Calibration bias (after averaging out MC noise) should be small.
  # Empirically ~0.005 at this Nmax/d; allow generous headroom for CI variance.
  expect_lt(max(abs(cov_mean - sigma.target)), 0.02)
})
