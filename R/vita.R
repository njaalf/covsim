#' Calibrate a regular vine
#'
#' \code{vita} implements the VITA (VIne-To-Anything) algorithm.
#' Covariance matrix and margins are specified, and \code{vita} calibrates the
#' pair-copulas in each node of the tree to match the target covariance.

#' @param margins A list where each element corresponds to a margin. Each
#' margin element is a list containing the distribution family ("distr") and
#' additional parameters. Must be a distribution available in the stats package.
#' @param sigma.target The target covariance matrix that is to be matched. The diagonal
#' elements must contain the variances of marginal distributions.
#' @param vc A vine dist object as specified by the rvinecopulib package. This object
#' specifies the  vine that is to be calibrated. If not provided, a D-vine is assumed.
#' @param family_set A vector of one-parameter pair-copula families that is to
#' be calibrated at each node in the vine. Possible entries are "gauss", "clayton", "joe", "gumbel" and "frank".
#' Calibration of pair-copula families is attempted in the order provided.
#' @param Nmax The sample size used for calibration. Reduce for faster calibration,
#' at the cost of precision.
#' @param numrootpoints The number of estimated roots at the initial calibration stage, which
#' determines a search interval where Nmax samples are drawn
#' @param conflevel Confidence level for determining search interval
#' @param numpoints The number of samples drawn with size Nmax, to determine the root within search interval
#' To increase precision increase this number. To calibrate faster (but less precisely), may be reduced to a number no lower than 2
#'@param verbose If TRUE, outputs details of calibration of each bicopula
#'@param cores Number of cores to use. If larger than 1, computations are done in parallel. May be determined with parallel:detectCores()
#'@return If a feasible solution exists,
#' \code{vita} returns a vine object which may be used for simulation.
#'
#'@references Grønneberg, S and Foldnes, N. (2017). Covariance model simulation using regular vines.
#' Psychometrika, 82(4), 1035-1051
#'
#' @examples
#' set.seed(1)# define a target covariance. 3 dimensions.
#' sigma.target <- cov(MASS::mvrnorm(10, mu=rep(0,3), Sigma=diag(1, 3)))
#'
#' #normal margins that match the covariances:
#' marginsnorm <- lapply(X=sqrt(diag(sigma.target)),function(X) list(distr="norm", sd=X) )
#'
#' #calibrate with a default D-vine, with rather low precision (default Nmax is 10^6)
#' # if cores=1 is removed, all cores are used, with a speed gain
#' cv <- vita(marginsnorm, sigma.target =sigma.target, Nmax=10^5, cores=1)
#'
#' #check
#' #round(cov(rvinecopulib::rvine(10^5, cv))-sigma.target, 3)
#'
#' #margins are normal but dependence structure is not
#' #pairs(rvinecopulib::rvine(500, cv))
#'
#'
#'
#' @export
vita <- function(margins, sigma.target, vc = NULL,
                 family_set = c("clayton", "gauss",
                                "joe", "gumbel", "frank"),
                 Nmax = 10^6, numrootpoints=10, conflevel=0.995, numpoints=4, verbose=TRUE, cores=parallel::detectCores())
{
  # if vc is provided we try to calibrate with its pair-copula families,
  # if a pc is not feasible, we first rotate and then if not succesful,
  # switch to pc's provided in family_set if vc not provided, we generate
  # a d-vine vc and run through family_set in each node, looking for
  # feasibility:
  d <- ncol(sigma.target)
  if (is.null(vc))
  {
    # vc d-vine with pcs from first member of family_set
    dstructure <- rvinecopulib::dvine_structure(1:d)  #d-vine
    pcs <- unflatten(rep(list(rvinecopulib::bicop_dist(family = family_set[1])),
                         d * (d - 1)/2))
    vc <- rvinecopulib::vinecop_dist(pair_copulas = pcs, structure = dstructure)

  }

  vine_structure <- rvinecopulib::get_structure(vc)
  pcs <- rvinecopulib::get_all_pair_copulas(vc)
  pcs_list <- unlist(pcs, recursive = F)  #flattened
  d <- dim(vine_structure)[1]
  v_matrix <- rvinecopulib::as_rvine_matrix(vine_structure)[1:d, ]
  pair_idx <- get_pair_idx(v_matrix)
  Matrix <- v_matrix[rev(1:d), ]  # VineCopula package, for create.submatrix

  # candidate copulas allowed
  if (!all(family_set %in% c("gauss", "clayton", "gumbel", "frank", "joe")))
    stop(paste("ERR: 'family_set' allows only one-parameter families: gauss, clayton, gumbel, frank, or joe.\n"))

  # run algorithm
  counter <- 1
  for (i in seq_along(pcs))
  {
    if (verbose)
      cat("Tree", i, "\n")
    for (j in seq_along(pcs[[i]]))
    {
      var1 <- v_matrix[d + 1 - j, j]
      var2 <- v_matrix[i, j]
      if (i == 1)
      {
        conditional <- NULL
      } else
      {
        conditional <- v_matrix[1:(i - 1), j]
      }
      pair.index <- c(var1, var2)
      cond.index <- conditional
      if(verbose)
        cat("   ", var1, "-", var2, "(", counter, "of", d * (d - 1)/2,
          ")\n")

      res <- tryCatch(solve.param(sigma.target, pair.index, cond.index, Matrix,
                                  margins, pair_idx, pcs_list, family_set, Nmax=Nmax,
                                  numrootpoints=numrootpoints, conflevel=conflevel,
                                  numpoints=numpoints, cores=cores), error = function(err)
                                  {
                                    print(paste("\n Error message in solve.param: ", err))
                                    return(NA)
                                  })
      if (length(res) !=2 )
        stop("\n \n  The specified vine, marginal and covariances are not compatible. \n \n ")
      pcs_list[res[[1]]] <- res[[2]]
      counter <- counter + 1

    }
  }

  pcs_calibrated <- unflatten(pcs_list)

  rvinecopulib::vine_dist(margins, pair_copulas = pcs_calibrated, structure = v_matrix)
}



