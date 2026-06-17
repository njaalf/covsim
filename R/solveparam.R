solve_param <- function(sigma.target, pair.index, cond.index, Matrix, margins, pair_idx,
                        pcs_list, family_set, Nmax, numrootpoints, conflevel, numpoints, cores)
{
  # construct subvine
  I <- c(pair.index, cond.index)
  sub.matrix <- create.submatrix(I, Matrix)$sub.matrix
  sub.matrix <- sub.matrix[rev(1:ncol(sub.matrix)), ]
  sub.pcs <- get_pair_idx(sub.matrix)
  num_trees <- max(sub.pcs[, "tree"])
  copulalist <- vector(mode = "list", length = num_trees)
  for (tree in 1:num_trees)
  {
    tmp <- matrix(sub.pcs[sub.pcs[, 3] == tree, 1:2], ncol = 2)
    tmplist <- vector(mode = "list", length = nrow(tmp))
    for (i in 1:nrow(tmp))
    {
      varpair <- tmp[i, 1:2]
      idx <- which((varpair[1] == pair_idx[, 1] & varpair[2] == pair_idx[, 2]) |
                   (varpair[1] == pair_idx[, 2] & varpair[2] == pair_idx[, 1]))
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

  # ---- Pre-draw uniforms once for the whole pair-copula calibration ----
  # The lower-tree pair-copulas are fixed during this calibration; only the
  # top copula's family/rotation/parameter varies. By reusing the same
  # uniform inputs across all candidate evaluations, the root-finding
  # objective becomes a deterministic function of theta, removing Monte
  # Carlo noise and allowing a single uniroot() call.
  d_sub <- ncol(sub.matrix)
  Nsim <- Nmax
  U_cached <- matrix(stats::runif(Nsim * d_sub), Nsim, d_sub)

  # Precompute marginal inverse-CDF closures for the two variables of interest
  marg1 <- margins[[pair.index[1]]]
  par1 <- marg1[names(marg1) != "distr"]
  qfun1_name <- paste0("q", marg1$distr)
  marg2 <- margins[[pair.index[2]]]
  par2 <- marg2[names(marg2) != "distr"]
  qfun2_name <- paste0("q", marg2$distr)
  target_cov <- sigma.target[pair.index[1], pair.index[2]]

  apply_marg <- function(sim.sample) {
    par1[[length(par1) + 1]] <- sim.sample[, 1]
    par2[[length(par2) + 1]] <- sim.sample[, 2]
    cbind(do.call(qfun1_name, par1), do.call(qfun2_name, par2))
  }

  make_root_fn <- function() {
    function(theta) {
      if (subvine$pair_copulas[[num_trees]][[1]]$family == "indep")
        return(-target_cov)
      subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(theta)
      sim.sample <- rvinecopulib::inverse_rosenblatt(U_cached, subvine,
                                                     cores = cores)[, recoded.pair.index]
      vars <- apply_marg(sim.sample)
      stats::cov(vars)[1, 2] - target_cov
    }
  }

  # Deterministic uniroot wrapper. Returns numeric root or NA.
  try_root <- function() {
    fam <- subvine$pair_copulas[[num_trees]][[1]]$family
    bb <- get_lowerupper(fam)
    if (fam == "indep") {
      return(if (isTRUE(all.equal(target_cov, 0))) 0 else NA_real_)
    }
    root.function <- make_root_fn()
    # Evaluate endpoints; require sign change for uniroot
    f_lo <- tryCatch(root.function(bb[1]), error = function(e) NA_real_)
    f_hi <- tryCatch(root.function(bb[2]), error = function(e) NA_real_)
    if (is.na(f_lo) || is.na(f_hi) || f_lo * f_hi > 0)
      return(NA_real_)
    res <- tryCatch(stats::uniroot(root.function, lower = bb[1], upper = bb[2],
                                   f.lower = f_lo, f.upper = f_hi,
                                   tol = .Machine$double.eps^0.25),
                    error = function(e) NA)
    if (is.na(res[1])) NA_real_ else res$root
  }

  curr_family <- subvine$pair_copulas[[num_trees]][[1]]$family
  curr_rotation <- subvine$pair_copulas[[num_trees]][[1]]$rotation
  log.message <- paste("The specified bicopula is", curr_family,
                       "with rotation", curr_rotation, ".\n")

  root <- try_root()

  if (is.na(root)) {
    log.message <- c(log.message, "           Specified bicopula not feasible \n")
  }

  # Try alternative rotations of current family
  if (is.na(root) && !curr_family %in% c("gauss", "frank", "indep")) {
    log.message <- c(log.message, "           Trying to rotate ")
    for (addrot in c(90, 180, 270)) {
      subvine$pair_copulas[[num_trees]][[1]]$rotation <- (curr_rotation + addrot) %% 360
      log.message <- c(log.message, paste("-", (curr_rotation + addrot) %% 360))
      root <- try_root()
      if (!is.na(root)) break
    }
  }

  # Try alternative families
  if (is.na(root)) {
    for (family in family_set[family_set != curr_family]) {
      log.message <- c(log.message, paste("\n          Switching to family", family, "\n"))
      pc <- rvinecopulib::bicop_dist(family)
      bb_init <- get_lowerupper(family)
      pc$parameters <- matrix(mean(bb_init))
      copulalist <- subvine$pair_copulas
      d <- subvine$structure$d
      sub.matrix <- rvinecopulib::as_rvine_matrix(
        rvinecopulib::get_structure(subvine))[1:d, ]
      copulalist[[num_trees]][[1]] <- pc
      subvine <- rvinecopulib::vinecop_dist(pair_copulas = copulalist,
                                            structure = sub.matrix)

      rotations <- if (family %in% c("gauss", "frank")) 0 else c(0, 90, 180, 270)
      for (addrot in rotations) {
        log.message <- c(log.message, paste("            with rotation ", addrot))
        subvine$pair_copulas[[num_trees]][[1]]$rotation <- addrot
        root <- try_root()
        if (!is.na(root)) break
      }
      if (!is.na(root)) break
    }
  }

  if (is.na(root)) {
    warning(log.message)
    return(NA)
  }

  subvine$pair_copulas[[num_trees]][[1]]$parameters <- matrix(root)
  idx <- which((pair.index[1] == pair_idx[, 1] & pair.index[2] == pair_idx[, 2]) |
               (pair.index[1] == pair_idx[, 2] & pair.index[2] == pair_idx[, 1]))

  return(list(idx, subvine$pair_copulas[[num_trees]]))
}


# Legacy multi-stage stochastic root finder. Retained for backward compatibility
# but no longer called by solve_param, which now uses a deterministic single
# uniroot() over cached uniforms via inverse_rosenblatt().
rootSearch <- function(root.function, thetaLower, thetaUpper, Nmax, numrootpoints,
                       conflevel, numpoints, cores) {
  Nsmall <- 1.5 * 10^3
  get_root <- function(n, lower, upper) {
    res <- tryCatch(stats::uniroot(root.function, lower = lower, upper = upper,
                                   Nsim = n), error = function(err) NA)
    if (is.na(res[1])) NA else res$root
  }
  roots <- replicate(numrootpoints, get_root(Nsmall, thetaLower, thetaUpper))
  solution_exists <- sum(is.na(roots)) < 0.2 * numrootpoints
  if (solution_exists)
    roots <- roots[stats::complete.cases(roots)]
  while (!solution_exists) {
    if (sum(is.na(roots)) > 0.7 * length(roots)) return(NA)
    roots <- replicate(numrootpoints, get_root(20 * Nsmall, thetaLower, thetaUpper))
    if (sum(is.na(roots)) > 0.3 * length(roots)) return(NA)
    solution_exists <- TRUE
    roots <- roots[!is.na(roots)]
  }
  exclude <- which(roots < stats::median(roots) - 2 * stats::IQR(roots) |
                   roots > stats::median(roots) + 2 * stats::IQR(roots))
  if (length(exclude) > 0) roots <- roots[-exclude]
  confint <- stats::t.test(roots, conf.level = conflevel)$conf.int
  lowerStart <- max(thetaLower, confint[1])
  upperStart <- min(thetaUpper, confint[2])
  x <- seq(lowerStart, upperStart, length.out = numpoints)
  y <- sapply(x, root.function, Nsim = Nmax)
  form <- stats::lm(formula = y ~ poly(x, 2, raw = TRUE))
  coefs <- stats::coef(form)
  poly_fun <- function(x) coefs[1] + coefs[2] * x + coefs[3] * x^2
  rootFinal <- tryCatch(stats::uniroot(poly_fun, lower = lowerStart, upper = upperStart),
                        error = function(err) NA)
  if (is.na(rootFinal[1])) NA else rootFinal$root
}
