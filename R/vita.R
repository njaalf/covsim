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
#' at the cost of precision. Since covsim 1.2.0 calibration uses cached uniform
#' inputs (via \code{\link[rvinecopulib]{inverse_rosenblatt}}) so the
#' root-finding objective is deterministic and \code{Nmax} only controls the
#' Monte Carlo precision of the matched covariance, not the stability of root
#' search.
#' @param numrootpoints Unused since covsim 1.2.0 (retained for backward
#' compatibility). The previous multi-stage stochastic root search was replaced
#' by a single deterministic \code{\link[stats]{uniroot}} over cached uniforms.
#' @param conflevel Unused since covsim 1.2.0 (retained for backward compatibility).
#' @param numpoints Unused since covsim 1.2.0 (retained for backward compatibility).
#'@param verbose If TRUE, outputs details of calibration of each bicopula
#'@param cores Number of cores to use. If larger than 1, computations are done in parallel. May be determined with parallel:detectCores()
#'@return If a feasible solution was found, a vine to be used for simulation
#'
#'
#'
#' @references Grønneberg, S., Foldnes, N., & Marcoulides, K. M. (2021). covsim: An r package for simulating non-normal data for structural equation models using copulas. Journal of Statistical Software. doi:10.18637/jss.v102.i03
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
#' calibrated.vine <- vita(marginsnorm, sigma.target =sigma.target, Nmax=10^5, cores=1)
#' #check
#' #round(cov(rvinecopulib::rvine(10^5, calibrated.vine))-sigma.target, 3)
#'
#' #margins are normal but dependence structure is not
#' #pairs(rvinecopulib::rvine(500, calibrated.vine))
#'
#'
#'
#' @export
vita <- function(margins, sigma.target, vc = NULL,
                 family_set = c("clayton", "gauss",
                                "joe", "gumbel", "frank"),
                 Nmax = 10^6, numrootpoints=10, conflevel=0.995,
                 numpoints=4, verbose=TRUE,
                 cores=parallel::detectCores())
{
  # candidate copulas allowed
  if (!all(family_set %in% c("gauss", "clayton", "gumbel", "frank", "joe")))
    stop(paste("ERR: 'family_set' allows only one-parameter families: gauss, clayton, gumbel, frank, or joe.\n"))
  #family_set <- match.arg(family_set, several.ok = TRUE)


  # if vc is provided we try to calibrate with its pair-copula families,
  # if a pc is not feasible, we first rotate and then if not succesful,
  # switch to pc's provided in family_set if vc not provided, we generate
  # a d-vine vc and run through family_set in each node, looking for
  # feasibility


  d <- ncol(sigma.target)

  #check correct margins wrt rvinecopulib
  dstructure <- rvinecopulib::dvine_structure(1:d)  #d-vine
  pcs <- unflatten(rep(list(rvinecopulib::bicop_dist(family = family_set[1])),
                       d * (d - 1)/2))
  res <- tryCatch(rvinecopulib::vine_dist(margins, pair_copulas = pcs,
                                          structure = dstructure), error = function(err)
                                          {
                                            stop(paste('\n Margins must be of type
  "beta", "cauchy", "chisq", "exp", "f", "gamma",
  "lnorm", "logis", "norm", "t", "unif", "weibull".'))
                                          })


  if (is.null(vc))
  {
    # vc d-vine with pcs from first member of family_set
    vc <- rvinecopulib::vinecop_dist(pair_copulas = pcs, structure = dstructure)
  }

  vine_structure <- rvinecopulib::get_structure(vc)
  pcs <- rvinecopulib::get_all_pair_copulas(vc)
  pcs_list <- unlist(pcs, recursive = F)  #flattened
  d <- dim(vine_structure)[1]
  v_matrix <- rvinecopulib::as_rvine_matrix(vine_structure)[1:d, ]
  pair_idx <- get_pair_idx(v_matrix)
  Matrix <- v_matrix[rev(1:d), ]  # VineCopula package, for create.submatrix


  # run algorithm
  counter <- 1; messages <- NULL
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
      res <- tryCatch(solve_param(sigma.target=sigma.target,pair.index=pair.index,
                                  cond.index=cond.index,Matrix=Matrix,
                                  margins=margins,
                                  pair_idx=pair_idx,
                                  pcs_list=pcs_list,
                                  family_set=family_set, Nmax=Nmax,
                                  numrootpoints=numrootpoints,
                                  conflevel=conflevel,
                                  numpoints=numpoints,
                                  cores=cores), error = function(err)
                                  {
                                    warning(paste("\n Error message in solve_param: ", err))
                                    return(NA)
                                  })
      if (is.na(res[[1]])){
        message("\n \n  The specified vine, marginal and covariances are not compatible. \n \n ")
        return(NULL)
      }
      pcs_list[res[[1]]] <- res[[2]]
      counter <- counter + 1
    }
  }

  pcs_calibrated <- unflatten(pcs_list)

  return(rvinecopulib::vine_dist(margins, pair_copulas = pcs_calibrated, structure = v_matrix))
}


#' Calibrate a regular vine (VITA)
#'
#' Wraps \code{\link{vita}} and returns a calibration object that can be passed
#' to \code{\link[stats]{simulate}}. The returned object embeds the calibrated
#' \code{vine_dist} so it can also be used directly with the \pkg{rvinecopulib}
#' simulator if desired.
#'
#' @inheritParams vita
#' @return If a feasible solution was found, an object of class
#'   \code{c("covsim_calib_vita", "covsim_calib")} that can be passed to
#'   \code{\link[stats]{simulate}}; otherwise \code{NULL}.
#' @seealso \code{\link{vita}}, \code{\link[stats]{simulate}}
#' @examples
#' set.seed(1)
#' sigma.target <- cov(MASS::mvrnorm(10, mu = rep(0, 3), Sigma = diag(1, 3)))
#' marginsnorm <- lapply(X = sqrt(diag(sigma.target)),
#'                      function(X) list(distr = "norm", sd = X))
#' cal <- calibrate_vita(marginsnorm, sigma.target = sigma.target,
#'                       Nmax = 10^5, cores = 1)
#' samples <- stats::simulate(cal, nsim = 1, N = 1000)
#' @export
calibrate_vita <- function(margins, sigma.target, vc = NULL,
                           family_set = c("clayton", "gauss",
                                          "joe", "gumbel", "frank"),
                           Nmax = 10^6, numrootpoints = 10, conflevel = 0.995,
                           numpoints = 4, verbose = TRUE,
                           cores = parallel::detectCores())
{
  vd <- vita(margins = margins, sigma.target = sigma.target, vc = vc,
             family_set = family_set, Nmax = Nmax,
             numrootpoints = numrootpoints, conflevel = conflevel,
             numpoints = numpoints, verbose = verbose, cores = cores)
  if (is.null(vd))
    return(NULL)
  structure(
    list(
      vine = vd,
      margins = margins,
      sigma.target = sigma.target
    ),
    class = c("covsim_calib_vita", "covsim_calib")
  )
}


#' Simulate from a calibrated VITA vine
#'
#' S3 method for \code{\link[stats]{simulate}} that draws samples from a
#' \code{covsim_calib_vita} object produced by \code{\link{calibrate_vita}}.
#'
#' @param object A \code{covsim_calib_vita} object.
#' @param nsim Number of independent samples (datasets) to return.
#' @param seed Optional seed passed to \code{\link{set.seed}}.
#' @param N Number of observations per simulated dataset.
#' @param cores Number of cores passed to \code{\link[rvinecopulib]{rvine}}.
#' @param ... Unused.
#' @return A list of \code{nsim} data frames (or matrices), each with \code{N} rows.
#' @export
simulate.covsim_calib_vita <- function(object, nsim = 1, seed = NULL, N,
                                       cores = 1, ...) {
  if (missing(N))
    stop("Please specify N (number of observations per dataset).")
  if (!is.null(seed))
    set.seed(seed)
  lapply(seq_len(nsim),
         function(i) rvinecopulib::rvine(N, object$vine, cores = cores))
}
