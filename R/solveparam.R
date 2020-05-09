solve.param <- function(sigma.target, pair.index, cond.index, Matrix, margins, pair_idx,
                        pcs_list, family_set, Nmax, numrootpoints, conflevel, numpoints, cores)
{
  # construct subvine
  I <- c(pair.index, cond.index)
  sub.matrix <- create.submatrix(I, Matrix)$sub.matrix  #the sub-vine we need for simulation
  sub.matrix <- sub.matrix[rev(1:ncol(sub.matrix)), ]
  sub.pcs <- get_pair_idx(sub.matrix)
  num_trees <- max(sub.pcs[, "tree"])
  # find idx in list of copulas
  copulalist <- vector(mode = "list", length = num_trees)
  for (tree in 1:num_trees)
  {
    tmp <- matrix(sub.pcs[sub.pcs[, 3] == tree, 1:2], ncol = 2)
    tmplist <- vector(mode = "list", length = nrow(tmp))
    for (i in 1:nrow(tmp))
    {
      varpair <- tmp[i, 1:2]
      idx <- which((varpair[1] == pair_idx[, 1] & varpair[2] == pair_idx[,
                                                                         2]) | (varpair[1] == pair_idx[, 2] & varpair[2] == pair_idx[,
                                                                                                                                     1]))
      tmplist[[i]] <- pcs_list[[idx]]
    }
    copulalist[[tree]] <- tmplist
  }

  # rename sub.matrix
  recoded.index <- sort(unique(as.vector(sub.matrix)))[-1]
  recoded.sub <- as.vector(sub.matrix)
  recoded.pair.index <- pair.index
  for (i in (1:length(recoded.index)))
  {
    recoded.sub[which(recoded.sub == recoded.index[i])] <- i
    recoded.pair.index[which(recoded.pair.index == recoded.index[i])] <- i
  }
  sub.matrix <- matrix(recoded.sub, dim(sub.matrix)[1], dim(sub.matrix)[2])

  # create subvine
  subvine <- rvinecopulib::vinecop_dist(pair_copulas = copulalist, structure = sub.matrix)
  root.function <- function(theta, Nsim)
  {
    if(subvine$pair_copulas[[num_trees]][[1]]$family == "indep")
      return(- sigma.target[pair.index[1],
                            pair.index[2]])

    subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(theta)
    sim.sample <- rvinecopulib::rvinecop(Nsim, subvine, cores = cores)[,
                                                                       recoded.pair.index]
    marg <- margins[[pair.index[1]]]
    par <- marg[names(marg) != "distr"]
    par[[length(par) + 1]] <- sim.sample[, 1]
    var1 <-do.call(paste0("q", marg$distr), par)

    marg <- margins[[pair.index[2]]]
    par <- marg[names(marg) != "distr"]
    par[[length(par) + 1]] <- sim.sample[, 2]
    var2 <-do.call(paste0("q", marg$distr), par)

    stats::cov(cbind(var1, var2))[1, 2] - sigma.target[pair.index[1],
                                                        pair.index[2]]
  }

  curr_family <- subvine$pair_copulas[[num_trees]][[1]]$family
  curr_rotation <- subvine$pair_copulas[[num_trees]][[1]]$rotation
  log.message <- "The specified bicopula is", curr_family, "with rotation", curr_rotation, ".\n")

  bb <- get_lowerupper(subvine$pair_copulas[[length(subvine$pair_copulas)]][[1]]$family)
  root <- tryCatch(rootSearch(root.function, thetaLower = bb[1], thetaUpper = bb[2],
                              Nmax=Nmax, numrootpoints, conflevel, numpoints),
                   error = function(err)
                   {
                     print(err)
                     return(NA)
                   })
  if (!is.na(root[1]))
  {
    log.message <- c(log.message, root[2])
    log.message <- c(log.message , "        Calibrated copula parameter is", round(root,3), "\n"))
    subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(root)
    varpair <- pair.index
    idx <- which((pair.index[1] == pair_idx[, 1] & pair.index[2] ==
                    pair_idx[, 2]) | (pair.index[1] == pair_idx[, 2] & pair.index[2] ==
                                        pair_idx[, 1]))

    return(list(idx, subvine$pair_copulas[[num_trees]], log.message))
  }

  log.message <- c(log.message, "           Specified bicopula not feasible \n")
  continue <- TRUE
  if (isFALSE(curr_family %in% c("gauss", "frank", "indep")))
  {
    log.message <- c(log.message,"           Trying to rotate ")
    for (addrot in c(90, 180, 270))
    {
      subvine$pair_copulas[[num_trees]][[1]]$rotation <- (curr_rotation +
                                                            addrot)%%360

      log.message <- c(log.message, paste("-", (curr_rotation + addrot)%%360))
      root.function <- function(theta, Nsim)
      {
        if(subvine$pair_copulas[[num_trees]][[1]]$family == "indep")
          return(- sigma.target[pair.index[1],
                                pair.index[2]])
        subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(theta)
        sim.sample <- rvinecopulib::rvinecop(Nsim, subvine, cores = cores)[,
                                                                                             recoded.pair.index]
        marg <- margins[[pair.index[1]]]
        par <- marg[names(marg) != "distr"]
        par[[length(par) + 1]] <- sim.sample[, 1]
        var1 <-do.call(paste0("q", marg$distr), par)

        marg <- margins[[pair.index[2]]]
        par <- marg[names(marg) != "distr"]
        par[[length(par) + 1]] <- sim.sample[, 2]
        var2 <-do.call(paste0("q", marg$distr), par)

        stats::cov(cbind(var1, var2))[1, 2] - sigma.target[pair.index[1],
                                                    pair.index[2]]
      }

      root <- tryCatch(rootSearch(root.function, thetaLower = bb[1],
                                  thetaUpper = bb[2], Nmax=Nmax,
                                  numrootpoints, conflevel, numpoints), error = function(err)
                                  {
                                    print(err)
                                    return(NA)
                                  })
      if (!is.na(root[1]))
      {
        log.message <- c(log.message, root[2])
        log.message <- c(log.message,paste("            Calibrated copula parameter", round(root,3), "\n"))
        continue <- FALSE
        break  # out of addrot
      }
    }
  }
  if (continue)
  {
    # look for other families
    for (family in family_set[curr_family != family_set])
    {
      log.message <- c(log.message, paste("\n          Switching to family", family, "\n"))
      bb <- get_lowerupper(family)
      upper <- bb[2]
      lower <- bb[1]
      pc <- rvinecopulib::bicop_dist(family)
      pc$parameters <- matrix(mean(c(upper, lower)))
      copulalist <- subvine$pair_copulas
      d <- subvine$structure$d
      sub.matrix <- rvinecopulib::as_rvine_matrix(rvinecopulib::get_structure(subvine))[1:d,
                                                            ]
      copulalist[[num_trees]][[1]] <- pc
      subvine <- rvinecopulib::vinecop_dist(pair_copulas = copulalist, structure = sub.matrix)

      rotations <- c(0, 90, 180, 270)
      # exception for gauss:
      if (family %in% c("gauss", "frank"))
        rotations <- c(0)

      for (addrot in rotations)
      {

        log.message <- c(log.message,paste("            with rotation ", addrot))
        subvine$pair_copulas[[num_trees]][[1]]$rotation <- addrot
        root.function <- function(theta, Nsim)
        {
          if(subvine$pair_copulas[[num_trees]][[1]]$family == "indep")
            return(- sigma.target[pair.index[1],
                                  pair.index[2]])
          subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(theta)
          sim.sample <- rvinecopulib::rvinecop(Nsim, subvine, cores = cores)[,
                                                                                               recoded.pair.index]
          marg <- margins[[pair.index[1]]]
          par <- marg[names(marg) != "distr"]
          par[[length(par) + 1]] <- sim.sample[, 1]
          var1 <-do.call(paste0("q", marg$distr), par)

          marg <- margins[[pair.index[2]]]
          par <- marg[names(marg) != "distr"]
          par[[length(par) + 1]] <- sim.sample[, 2]
          var2 <-do.call(paste0("q", marg$distr), par)

          stats::cov(cbind(var1, var2))[1, 2] - sigma.target[pair.index[1],
                                                      pair.index[2]]
        }

        root <- tryCatch(rootSearch(root.function, thetaLower = bb[1],
                                    thetaUpper = bb[2], Nmax=Nmax,
                                    numrootpoints, conflevel, numpoints), error = function(err)
                                    {
                                      print(err)
                                      return(NA)
                                    })
        if (!is.na(root[1]))
        {
          log.message <- c(log.message, root[2])
          log.message <- c(log.message,paste("           Calibrated copula parameter", round(root,3), "\n"))
          continue <- FALSE
          break  # out of addrot
        }
      }  #addrot end
      if (!continue)
        break
    }  #end families
  }

  if (continue)
    return(NA)


  subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(root)
  varpair <- pair.index
  idx <- which((pair.index[1] == pair_idx[, 1] & pair.index[2] == pair_idx[,
                                                                           2]) | (pair.index[1] == pair_idx[, 2] & pair.index[2] == pair_idx[,                                                                                                                                    1]))

  return(list(idx, subvine$pair_copulas[[num_trees]], log.message))
}


rootSearch <- function(root.function, thetaLower, thetaUpper,Nmax, numrootpoints , conflevel, numpoints){
  Nsmall <- 1.5*10^3
  get_root <- function(n, lower, upper)
  {
    res <- tryCatch(stats::uniroot(root.function, lower = lower, upper = upper,
                            Nsim = n), error = function(err)
                            {
                              return(NA)
                            })
    if (is.na(res[1]))
      return(NA)
    else
      return(res$root)
  }

  # Estimate standard error of root
  roots <- replicate(numrootpoints, get_root(Nsmall,thetaLower, thetaUpper))
  solution_exists <- sum(is.na(roots)) < 0.2*numrootpoints
  if(solution_exists)# remove the few NAs
    roots <- roots[stats::complete.cases(roots)]

  while (!solution_exists)
  {
      if(sum(is.na(roots)) > 0.7*length(roots)){#give up
        return(NA)
      } else { # increase precision
        roots <- replicate(numrootpoints, get_root(20*Nsmall,thetaLower, thetaUpper))
        if(sum(is.na(roots)) > 0.3*length(roots)){
          return(NA)
        } else {
          solution_exists <- TRUE
          roots <- roots[!is.na(roots)]
        }
      }
  }

  #remove outliers
  exclude <- which(roots < stats::median(roots) -2*stats::IQR(roots) | roots > stats::median(roots) +2*stats::IQR(roots) )

  if(length(exclude)>0)
    roots <- roots[-exclude]
  confint <- stats::t.test(roots, conf.level = conflevel)$conf.int

  lowerStart <- max(thetaLower, confint[1])
  upperStart <- min(thetaUpper, confint[2])
  search.message <-paste("        Searching in interval (", round(lowerStart, 2),
        ",", round(upperStart, 2), ")\n")

  x <- seq(lowerStart, upperStart, length.out = numpoints)
  y <- sapply(x, root.function, Nsim = Nmax)
  form <- stats::lm(formula = y ~ poly(x, 2, raw = TRUE))
  coefs <- stats::coef(form)
  poly_fun <- function(x)
      {
        coefs[1] + coefs[2] * x + coefs[3] * x^2
      }

  rootFinal <- tryCatch(stats::uniroot(poly_fun, lower = lowerStart, upper = upperStart),
                        error = function(err)
                        {
                          return(NA)
                        })
  if(is.na(rootFinal[1]))
    return(NA)
  else
    return(list(rootFinal$root, search.message)
}
