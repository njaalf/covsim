create.submatrix <- function(I, Matrix)
{
  # I <- unique(I)
  d <- dim(Matrix)[1]
  l <- length(I)

  for (i in (d:1))
  {
    if(sum(Matrix[, i] %in% I) == l)
      break
  }

  sub.matrix <- Matrix[(i:d), (i:d)]

  remove.ind <- which(!(diag(sub.matrix) %in% I))
  if (length(remove.ind) != 0)
  {
    sub.matrix <- sub.matrix[, -remove.ind]
  }

  new.sub.matrix <- NULL

  for (i in (1:(dim(sub.matrix)[2])))
  {
    indx <- which(sub.matrix[, i] %in% I)

    tmp.mat <- sub.matrix[indx, i]
    new.sub.matrix <- cbind(new.sub.matrix, c(rep(0, length(I) - length(tmp.mat)),
                                              tmp.mat))
  }

  sub.matrix <- new.sub.matrix
  return(list(sub.matrix = sub.matrix))
}

## functions that return the order of pcs in vines structure recursive!
get_uncond_pcs <- function(M)
{
  a <- cbind(diag(M[nrow(M):1, ]), M[1, ])
  return(a[-nrow(a), ])
}
get_pair_idx <- function(M)
{
  res <- matrix(get_uncond_pcs(M), ncol = 2)
  tree <- 1
  res <- cbind(res, tree)
  if (nrow(M) > 2)
  {
    for (tree in 2:(nrow(M) - 1))
    {
      tmpM <- M[tree:nrow(M), 1:length(tree:nrow(M))]
      tmpM <- matrix(get_uncond_pcs(tmpM), ncol = 2)
      tmpM <- cbind(tmpM, tree)
      res <- rbind(res, tmpM)
    }
  }
  res
}
# un-flatten
unflatten <- function(longlist)
{
  position <- length(longlist)
  ll <- 1
  mylist <- list()
  while (position > 0)
  {
    mylist <- c(mylist, list(longlist[(position - ll + 1):position]))
    position <- position - ll
    ll <- ll + 1
  }
  rev(mylist)
}

## adf test 2. We first estimate W Matrix from large sample, and use it
## repeatedly.
# test_covariance2 <- function(sigma.target, vc, n = 10^5, bign = 10^6, reps = 100)
# {
#   model <- NULL
#   for (i in 1:ncol(sigma.target)) model <- paste(model, paste0("x", i,
#                                                                "~~", diag(sigma.target)[i], "*x", i, "\n"))
#   for (i in 2:nrow(sigma.target))
#   {
#     for (j in 1:(i - 1)) model <- paste(model, paste0("x", i, "~~",
#                                                       sigma.target[i, j], "*x", j, "\n"))
#   }
#   ## large WLS weight matrix
#   bigsample <- data.frame(rvine(bign, vc, cores=parallel::detectCores()))
#   wls.V <- solve(lavaan:::lavGamma(bigsample)[[1]])
#   wls.0 <- lavaan::lav_matrix_vech(sigma.target)
#   pvalues <- numeric(reps)
#   set.seed(1)
#   pb <- txtProgressBar(min = 0, max = reps, style = 3)
#   for (i in 1:reps)
#   {
#     wls.obs <- lavaan::lav_matrix_vech(cov(rvinecopulib::rvine(n, vc,cores=parallel::detectCores())))
#
#     wls.res <- as.matrix(wls.obs - wls.0)
#     TestStat <- (n - 1) * as.numeric(t(wls.res) %*% wls.V %*% wls.res)
#     pvalues[i] <- 1 - pchisq(TestStat, df = length(wls.obs))
#     setTxtProgressBar(pb, i)
#   }
#   # pvalues should be uniformly distributed
#   res <- ks.test(pvalues, y = "punif")
#   res$p.value
# }



# one-parameter copulas
get_lowerupper <- function(name)
{
  names <- c("gaussian", "clayton", "gumbel", "frank", "joe", "indep")
  bound_list <- list(c(-1, 1), c(0.0001, 28), c(1, 50), c(-35, 34.8), c(1,
                                                                      30), c(0,0))
  bound_list[[grep(name, names, fixed = T)]]
}




