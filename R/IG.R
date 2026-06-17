#' Calibrate the IG (Independent Generator) model
#'
#' Performs the calibration step of the IG algorithm: computes the transformation
#' matrix \code{A} and fits a Pearson distribution to each independent generator
#' such that samples drawn from the calibrated model match the target covariance,
#' skewness, and excess kurtosis. Calibration is decoupled from simulation so the
#' same calibration object can be reused to draw many samples without recomputation.
#'
#' @param sigma.target Target population covariance matrix.
#' @param skewness Target skewness vector.
#' @param excesskurtosis Target excess kurtosis vector.
#' @param typeA Either \code{"triang"} (default, triangular Cholesky factor) or
#'   \code{"symm"} (symmetric square root).
#' @return An object of class \code{c("covsim_calib_ig", "covsim_calib")} that can
#'   be passed to \code{\link[stats]{simulate}} to draw samples.
#' @seealso \code{\link{rIG}}, \code{\link[stats]{simulate}}
#' @author Njål Foldnes  (\email{njal.foldnes@gmail.com})
#' @references Foldnes, N. and Olson, U. H. (2016). A simple simulation technique for nonnormal data with prespecified skewness,
#'             kurtosis, and covariance matrix. Multivariate behavioral research, 51(2-3), 207-219
#' @examples
#' set.seed(1)
#' sigma.target <- cov(MASS::mvrnorm(10, mu = rep(0, 3), Sigma = diag(3)))
#' cal <- calibrate_ig(sigma.target, skewness = rep(1, 3), excesskurtosis = rep(3, 3))
#' samples <- stats::simulate(cal, nsim = 5, N = 1000)
#' @export
calibrate_ig <- function(sigma.target, skewness, excesskurtosis, typeA = "triang") {
  if (is.null(sigma.target))
    stop("Please specify sigma.target")

  nvar <- dim(sigma.target)[2]

  if (typeA == "symm")
    A <- lavaan::lav_matrix_symmetric_sqrt(sigma.target)
  else
    A <- t(chol(sigma.target))

  function.skew <- function(IGvalues) {
    fval <- numeric(nvar)
    for (i in 1:nvar)
      fval[i] <- A[i, ]^3 %*% IGvalues / (sum(A[i, ]^2)^(3/2))
    fval - skewness
  }
  function.kurt <- function(IGvalues) {
    fval <- numeric(nvar)
    for (i in 1:nvar)
      fval[i] <- A[i, ]^4 %*% IGvalues / (sum(A[i, ]^2)^(2))
    fval - excesskurtosis
  }

  IGskew        <- nleqslv::nleqslv(x = skewness, function.skew)$x
  IGkurt.excess <- nleqslv::nleqslv(x = excesskurtosis, function.kurt)$x

  parlist <- list()
  for (i in 1:nvar) {
    parlist[[i]] <- tryCatch(
      PearsonDS::pearsonFitM(moments = c(mean = 0, variance = 1,
                                         skewness = IGskew[i],
                                         3 + IGkurt.excess[i])),
      error = function(c) {
        stop(paste("No valid solution for variable", i,
                   "- Try setting typeA==\"symm\"."))
      })
  }

  structure(
    list(
      A = A,
      parlist = parlist,
      sigma.target = sigma.target,
      skewness = skewness,
      excesskurtosis = excesskurtosis,
      typeA = typeA,
      nvar = nvar
    ),
    class = c("covsim_calib_ig", "covsim_calib")
  )
}


#' Simulate from a calibrated IG model
#'
#' S3 method for \code{\link[stats]{simulate}} that draws samples from a
#' \code{covsim_calib_ig} object produced by \code{\link{calibrate_ig}}.
#'
#' @param object A \code{covsim_calib_ig} object.
#' @param nsim Number of independent samples (datasets) to return.
#' @param seed Optional seed passed to \code{\link{set.seed}}.
#' @param N Number of observations per simulated dataset.
#' @param ... Unused.
#' @return A list of \code{nsim} matrices, each with \code{N} rows.
#' @export
simulate.covsim_calib_ig <- function(object, nsim = 1, seed = NULL, N, ...) {
  if (missing(N))
    stop("Please specify N (number of observations per dataset).")
  if (!is.null(seed))
    set.seed(seed)

  IGdata <- sapply(seq_len(object$nvar),
                   function(i) PearsonDS::rpearson(N * nsim, object$parlist[[i]]))

  simulated.samples <- IGdata %*% t(object$A)
  colnames(simulated.samples) <- colnames(object$sigma.target)
  idx <- rep(seq_len(nsim), each = N)
  lapply(seq_len(nsim), function(i) simulated.samples[i == idx, , drop = FALSE])
}


#' Simulation of non-normal data
#'
#' Backward-compatible wrapper that combines \code{\link{calibrate_ig}} and
#' \code{\link[stats]{simulate}} in one call. New code should prefer the
#' calibrate/simulate split so calibration can be reused.
#'
#' @inheritParams calibrate_ig
#' @param N Number of observations to simulate.
#' @param reps Number of simulated samples.
#' @return A list of simulated samples.
#' @author Njål Foldnes  (\email{njal.foldnes@gmail.com})
#' @references Foldnes, N. and Olson, U. H. (2016). A simple simulation technique for nonnormal data with prespecified skewness,
#'             kurtosis, and covariance matrix. Multivariate behavioral research, 51(2-3), 207-219
#' @examples
#' set.seed(1234)
#' model <- '
#'  # measurement model
#'    ind60 =~ x1 + x2 + x3
#'    dem60 =~ y1 + y2 + y3 + y4
#'    dem65 =~ y5 + y6 + y7 + y8
#'  # regressions
#'    dem60 ~ ind60
#'    dem65 ~ ind60 + dem60
#'  # residual correlations
#'    y1 ~~ y5
#'    y2 ~~ y4 + y6
#'    y3 ~~ y7
#'    y4 ~~ y8
#'    y6 ~~ y8'
#' fit  <- lavaan::sem(model, data=lavaan::PoliticalDemocracy)
#' population.sigma <- lavaan::lavInspect(fit, "sigma.hat")
#' population.skew  <- c(0, 0, 0, 0, 1, 1, 1, 1, 2,2,2 )
#' population.excesskurt <- c( 1 , 1, 1, 1, 3, 3, 3, 3, 15, 15, 15)
#' my.samples <- rIG(N=10^3, sigma=population.sigma,
#'         skewness=population.skew,
#'         excesskurt=population.excesskurt,
#'         reps=5)
#' @export
rIG <- function(N, sigma.target, skewness, excesskurtosis, reps = 1, typeA = "triang") {
  cal <- calibrate_ig(sigma.target, skewness, excesskurtosis, typeA)
  stats::simulate(cal, nsim = reps, N = N)
}
