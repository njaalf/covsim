#' Simulation of non-normal data
#'
#' Using the IG method to simulate non-normal data
#'
#'
#' @param N Number of observations to simulate.
#' @param sigma.target Target population covariance matrix
#' @param skewness Target skewness
#' @param excesskurtosis Target excess kurtosis
#' @param reps Number of simulated samples
#' @param typeA Symmetrical or triangular (default) A matrix
#' @return A list of simulated samples
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
rIG <- function(N, sigma.target,  skewness, excesskurtosis, reps=1, typeA="triang") {
  if(  is.null(sigma.target))
    stop("Please specify sigma.target")

  nvar          <- dim(sigma.target)[2]
  #define functions

  function.skew <- function(IGvalues){
    fval <- numeric(nvar)
    for (i in 1:nvar)
      fval[i] <-   A[i, ]^3 %*% IGvalues/(sum(A[i,]^2)^(3/2))

    fval-skewness
  }
  function.kurt <- function(IGvalues){
    fval <- numeric(nvar)
    for (i in 1:nvar)
      fval[i] <-   A[i, ]^4 %*% IGvalues/(sum(A[i,]^2)^(2))

    fval-excesskurtosis
  }

  #calculate A
  if(typeA=="symm")
    A  <- lavaan::lav_matrix_symmetric_sqrt(sigma.target)#symmetric
  else
    A  <- t(chol(sigma.target))

  IGskew        <- nleqslv::nleqslv(x=skewness, function.skew)$x
  IGkurt.excess <- nleqslv::nleqslv(x=excesskurtosis, function.kurt)$x
  parlist       <- list()
  for (i in 1:nvar){
    res <- tryCatch(PearsonDS::pearsonFitM(moments=c(mean=0, variance=1, skewness=IGskew[i], 3+IGkurt.excess[i])),
                    error = function(c) {
                      stop(paste("No valid solution for variable", i, "- Try setting typeA==\"symm\"."))})
    parlist[[i]] <- res
  }

  IGdata <- sapply(1:nvar, function(i) PearsonDS::rpearson(N*reps, parlist[[i]]) )

  simulated.samples           <- IGdata %*% t(A)
  colnames(simulated.samples) <- colnames(sigma.target)
  idx <- rep(1:reps, each=N)
  return(lapply(1:reps, function(i) simulated.samples[i==idx, ]))

}
##################################################
