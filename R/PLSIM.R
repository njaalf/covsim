#' Calibrate the PLSIM (Piecewise Linear) model
#'
#' Performs the calibration step of the PLSIM algorithm: fits piecewise linear
#' marginal transforms with the target skewness and kurtosis and searches for
#' intermediate normal correlations that yield the target covariance after
#' transformation. The (expensive) pairwise calibration is decoupled from
#' simulation so the resulting calibration object can be reused to draw many
#' samples without recomputation.
#'
#' @param sigma.target Target population covariance matrix.
#' @param skewness Target skewness vector.
#' @param excesskurtosis Target excess kurtosis vector.
#' @param numsegments The number of line segments in each marginal.
#' @param gammalist A list of breakpoints in each margin.
#' @param monot If \code{TRUE}, piecewise linear functions are forced to be monotonous.
#'   The implied copula is then normal.
#' @param verbose If \code{TRUE}, prints calibration progress.
#' @return An object of class \code{c("covsim_calib_plsim", "covsim_calib")} that
#'   can be passed to \code{\link[stats]{simulate}} to draw samples.
#' @seealso \code{\link{rPLSIM}}, \code{\link[stats]{simulate}}
#' @author Njål Foldnes  (\email{njal.foldnes@gmail.com})
#' @references Foldnes, N. and Grønneberg S. (2021). Non-normal data simulation using piecewise linear transforms.
#' @examples
#' set.seed(1)
#' sigma.target <- cov(MASS::mvrnorm(5, rep(0, 3), diag(3)))
#' cal <- calibrate_plsim(sigma.target, skewness = rep(1, 3), excesskurtosis = rep(4, 3))
#' samples <- stats::simulate(cal, nsim = 2, N = 1000)
#' @export
calibrate_plsim <- function(sigma.target, skewness, excesskurtosis,
                            numsegments = 4, gammalist = NULL,
                            monot = FALSE, verbose = TRUE) {
  if (is.null(sigma.target))
    stop("Please specify sigma.target")
  sigma.target <- as.matrix(sigma.target)
  nvar <- dim(sigma.target)[2]
  corr.target <- stats::cov2cor(sigma.target)
  sd_vec <- sqrt(diag(sigma.target))

  if (is.null(gammalist)) {
    gamma <- get_gamma(numsegments)
    gammalist <- rep(list(gamma), ncol(corr.target))
  }

  # univariate special case
  if (nvar == 1) {
    out <- stats::nlminb(start = rep(1, length(gammalist[[1]]) + 1),
                         objective = skewkurt_discrepancy,
                         gamma = gammalist[[1]], skew = skewness[1],
                         kurt = excesskurtosis[1],
                         lower = ifelse(monot, 0.1, -Inf))
    if (out$objective > 1e-4)
      return(NA)
    a <- out$par
    b <- get_bs(a, gammalist[[1]])
    mom <- get_pl_moments(a, b, gammalist[[1]])
    a <- a / sqrt(mom[2])
    b <- get_bs(a, gammalist[[1]])

    return(structure(
      list(
        alist = list(a),
        blist = list(b),
        gammalist = gammalist,
        z.corr = matrix(1, 1, 1),
        pre = matrix(1, 1, 1),
        sd = sd_vec,
        sigma.target = sigma.target,
        skewness = skewness,
        excesskurtosis = excesskurtosis,
        nvar = 1L,
        monot = monot,
        univariate = TRUE
      ),
      class = c("covsim_calib_plsim", "covsim_calib")
    ))
  }

  # general multidimensional case
  alist <- fit_univariate(gammalist, skewness, excesskurtosis,
                          scale = TRUE, monot = monot)

  if (sum(is.na(alist)) > 0)
    stop("Error: Can not fit univariate marginals for variable(s): ",
         which(is.na(alist)),
         " . Increase numsegments, set monot=FALSE, or provide gammalist.\n")

  blist <- lapply(seq_along(alist), function(i) get_bs(alist[[i]], gammalist[[i]]))

  # pairwise calibration
  z.corr <- diag(length(alist))
  solved.already <- list()
  solved.corrs <- NULL
  for (i in 1:(ncol(z.corr) - 1)) {
    for (j in ((i + 1):ncol(z.corr))) {
      if (verbose)
        message("Calibrating vars: ", i, "-", j)
      listelement <- list(corr.target[i, j], alist[[i]], blist[[i]], gammalist[[i]],
                          alist[[j]], blist[[j]], gammalist[[j]])
      if (length(solved.already) > 0)
        match <- which(sapply(solved.already, identical, listelement))
      else
        match <- vector()

      if (length(match) == 0) {
        out <- stats::nlminb(corr.target[i, j], function(rho) {
          (get_cov(alist[[i]], blist[[i]], gammalist[[i]],
                   alist[[j]], blist[[j]], gammalist[[j]], rho) -
             corr.target[i, j])^2
        }, lower = -0.998, upper = 0.998)
        if (out$objective < 1e-4) {
          z.corr[i, j] <- z.corr[j, i] <- out$par
          solved.already[[length(solved.already) + 1]] <- listelement
          solved.corrs <- c(solved.corrs, out$par)
        } else {
          stop("Error: No intermediate correlation exists for ", i, j, "\n")
        }
      } else {
        z.corr[i, j] <- z.corr[j, i] <- solved.corrs[match]
      }
    }
  }

  pre <- diag(ncol(z.corr))
  if (!all(eigen(z.corr)$values > 0)) {
    message("The intermediate matrix is not positive definite. It is corrected.
        Skewness and kurtosis values will not be exactly matched. \n")
    z.corr <- Matrix::nearPD(z.corr, corr = TRUE)$mat
    M <- z.corr
    for (i in 1:(ncol(z.corr) - 1)) {
      for (j in ((i + 1):ncol(z.corr))) {
        M[i, j] <- M[j, i] <- get_cov(alist[[i]], blist[[i]], gammalist[[i]],
                                      alist[[j]], blist[[j]], gammalist[[j]],
                                      z.corr[i, j])
      }
    }
    pre <- lavaan::lav_matrix_symmetric_sqrt(corr.target) %*%
      solve(lavaan::lav_matrix_symmetric_sqrt(M))
  }

  structure(
    list(
      alist = alist,
      blist = blist,
      gammalist = gammalist,
      z.corr = z.corr,
      pre = pre,
      sd = sd_vec,
      sigma.target = sigma.target,
      skewness = skewness,
      excesskurtosis = excesskurtosis,
      nvar = as.integer(nvar),
      monot = monot,
      univariate = FALSE
    ),
    class = c("covsim_calib_plsim", "covsim_calib")
  )
}


#' Simulate from a calibrated PLSIM model
#'
#' S3 method for \code{\link[stats]{simulate}} that draws samples from a
#' \code{covsim_calib_plsim} object produced by \code{\link{calibrate_plsim}}.
#'
#' @param object A \code{covsim_calib_plsim} object.
#' @param nsim Number of independent samples (datasets) to return.
#' @param seed Optional seed passed to \code{\link{set.seed}}.
#' @param N Number of observations per simulated dataset.
#' @param ... Unused.
#' @return A list of \code{nsim} matrices, each with \code{N} rows.
#' @export
simulate.covsim_calib_plsim <- function(object, nsim = 1, seed = NULL, N, ...) {
  if (missing(N))
    stop("Please specify N (number of observations per dataset).")
  if (!is.null(seed))
    set.seed(seed)

  if (isTRUE(object$univariate)) {
    return(lapply(seq_len(nsim), function(i) {
      x <- pl_fun(stats::rnorm(N), object$alist[[1]],
                  object$blist[[1]], object$gammalist[[1]])
      matrix(x * object$sd, ncol = 1)
    }))
  }

  samples <- lapply(seq_len(nsim), function(i) {
    z <- MASS::mvrnorm(N, rep(0, object$nvar), Sigma = object$z.corr)
    sapply(seq_len(object$nvar), function(j) {
      pl_fun(z[, j], object$alist[[j]], object$blist[[j]], object$gammalist[[j]])
    })
  })
  lapply(samples, function(x) x %*% t(object$pre) %*% diag(object$sd))
}


#' Simulation of non-normal data
#'
#' Backward-compatible wrapper that combines \code{\link{calibrate_plsim}} and
#' \code{\link[stats]{simulate}} in one call. New code should prefer the
#' calibrate/simulate split so the (expensive) pairwise calibration can be reused.
#'
#' @param N Number of observations to simulate.
#' @inheritParams calibrate_plsim
#' @param reps Number of simulated samples.
#' @return A list with two elements. First element: the list of simulated samples.
#'   Second element: the fitted piecewise linear functions and the intermediate
#'   correlations matrix.
#' @author Njål Foldnes  (\email{njal.foldnes@gmail.com})
#' @references Foldnes, N. and Grønneberg S. (2021). Non-normal data simulation using piecewise linear transforms.Under review.
#' @examples
#' set.seed(1)
#' sigma.target  <- cov(MASS::mvrnorm(5, rep(0,3), diag(3)))
#' res  <- covsim::rPLSIM(10^5, sigma.target, skewness=rep(1,3), excesskurtosis=rep(4,3))
#' my.sample  <- res[[1]][[1]]
#' @export
rPLSIM <- function(N, sigma.target, skewness, excesskurtosis, reps = 1,
                   numsegments = 4, gammalist = NULL,
                   monot = FALSE, verbose = TRUE) {
  cal <- calibrate_plsim(sigma.target, skewness, excesskurtosis,
                         numsegments = numsegments, gammalist = gammalist,
                         monot = monot, verbose = verbose)
  if (identical(cal, NA))
    return(NA)
  samples <- stats::simulate(cal, nsim = reps, N = N)
  list(samples = samples,
       model = list(a = cal$alist, b = cal$blist,
                    gamma = cal$gammalist, z.corr = cal$z.corr))
}
