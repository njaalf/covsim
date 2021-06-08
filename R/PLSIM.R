#' Simulation of non-normal data
#'
#' Using the piecewise linear PLSIM method to simulate non-normal data
#'
#'
#' @param N Number of observations to simulate.
#' @param sigma.target Target population covariance matrix
#' @param skewness Target skewness
#' @param excesskurtosis Target excess kurtosis
#' @param reps Number of simulated samples
#' @param numsegments The number of line segments in each marginal
#' @param gammalist A list of breakpoints in each margin
#' @param monot  True if piecewise linear functions are forced to be monotonous. The copula will then be normal.
#' @param verbose If true, progress details of the procedure are printed
#' @return A list with two elements. First element: the list of simulated samples. Second element: The fitted piecewise linear functions and the intermediate correlations matrix.
#' @author Njål Foldnes  (\email{njal.foldnes@gmail.com})
#' @references Foldnes, N. and Grønneberg S. (2021). Non-normal data simulation using piecewise linear transforms. Structural Equation Modeling, Volume (2-3), pp-pp
#' @examples
#' #'set.seed(1)
#'sigma.target  <- cov(MASS::mvrnorm(5, rep(0,3), diag(3)))
#'res  <- covsim::rPLSIM(10^5, sigma.target, skewness=rep(1,3), excesskurtosis=rep(4,3))
#'my.sample  <- res[[1]][[1]]
#' @export
rPLSIM <- function(N, sigma.target,  skewness, excesskurtosis, reps=1, numsegments=4, gammalist=NULL, monot=FALSE, verbose=TRUE) {
  if(  is.null(sigma.target))
    stop("Please specify sigma.target")
  sigma.target <-as.matrix(sigma.target)
  nvar <- dim(sigma.target)[2]
  corr.target <- cov2cor(sigma.target)
  if(is.null(gammalist)){
    gamma <- get_gamma(numsegments)
    gammalist <- rep(list(gamma), ncol(corr.target))
  }

  #univariate special case
  if(nvar==1){
    out <- nlminb(start=rep(1, length(gammalist[[1]])+1),
                  objective =skewkurt_discrepancy,
                  gamma=gammalist[[1]], skew=skewness[1], kurt=excesskurtosis[1], lower=ifelse(monot, 0.1,-Inf))
    if(out$objective > 1e-4)
      return(NA)
    a <- out$par; b  <-get_bs(a, gammalist[[1]])
    mom <- get_pl_moments(a, b, gammalist[[1]])
    a <- a/sqrt(mom[2]); b  <- get_bs(a, gammalist[[1]])

    samples <- lapply(1:reps, function(i) pl_fun(rnorm(N),a, b, gammalist[[1]]))

    return(list(samples=lapply(samples, function(x) x*sqrt(c(sigma.target))), model=list(a=list(a),b=list(b),gamma=gammalist)))
  }

  #general multidimensional case
  alist  <- fit_univariate(gammalist, skewness, excesskurtosis, scale=TRUE, monot=monot)

  if(sum(is.na(alist))>0)
    stop("Error: Can not fit univariate marginals for variable(s): ",
         which(is.na(alist)), " . Increase numsegments, set monot=FALSE, or provide gammalist.\n")
  #scale
  blist <- lapply(1:length(alist), function(i) get_bs(alist[[i]], gammalist[[i]]))

  #pairwise calibration
  z.corr <- diag(length(alist))
  for(i in 1:(ncol(z.corr)-1)){
    for(j in ((i+1):ncol(z.corr))){
      if(verbose)
        message("Calibrating vars: ",i,"-",j)
      out <- nlminb(corr.target[i, j], function(rho){
        (get_cov(alist[[i]],blist[[i]],gammalist[[i]],  alist[[j]],blist[[j]], gammalist[[j]], rho)-corr.target[i,j])^2
      }, lower = -0.998, upper=0.998)
      if(out$objective < 1e-4)
        z.corr[i, j] <- z.corr[j,i] <- out$par
      else {
        stop("Error: No intermediate correlation exists for ", i, j,"\n" )
      }
    }
  }

  correct <- FALSE
  pre <- diag(ncol(z.corr))
  if(!all(eigen(z.corr)$values > 0)){
    correct <- TRUE
    message("The intermediate matrix is not positive definite. It is corrected.
        Skewness and kurtosis values will not be exactly matched. \n")
    z.corr <- Matrix::nearPD(z.corr, corr=TRUE)$mat
    # calculate model-implied covariance matrix
    M <- z.corr
    for(i in 1:(ncol(z.corr)-1)){
      for(j in ((i+1):ncol(z.corr))){
        #cat(i, j, "\n")
        M[i, j] <- M[j,i] <- get_cov(alist[[i]],blist[[i]],gammalist[[i]],
                                     alist[[j]],blist[[j]], gammalist[[j]], z.corr[i,j])
      }
    }
    pre <- lavaan::lav_matrix_symmetric_sqrt(corr.target) %*% solve(lavaan::lav_matrix_symmetric_sqrt(M))
  }

  samples <- lapply(1:reps, function(i){
    z <- MASS::mvrnorm(N, rep(0, ncol(corr.target)), Sigma=z.corr)
    x <- sapply(1:ncol(z), function(i) pl_fun(z[, i], alist[[i]], blist[[i]], gammalist[[i]]) )
  })

  samples <- lapply(samples, function(x) x%*%t(pre) %*% diag(sqrt(diag(sigma.target))))

  return(list(samples=samples, model=list(a=alist,b=blist,gamma=gammalist, z.corr=z.corr)))
}




