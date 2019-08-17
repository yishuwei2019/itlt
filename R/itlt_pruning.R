NodeCount <- function(nd) {
  ifelse(is.null(nd) || is.terminal(nd),
         0,
         1 + NodeCount(nd$kids[[1]]) + NodeCount(nd$kids[[2]]))
}


GetGt <- function(nd, nknot.transform) {
  if (is.null(nd) || is.null(nd$split)) {
    return(0)
  } else if (nknot.transform != -1) {
    score <- Chi2Normalize(nd$split$info$score, nknot.transform)
  } else {
    score <- nd$split$info$score
  }

  score + GetGt(nd$kids[[1]], nknot.transform) + GetGt(nd$kids[[2]], nknot.transform)
}


GetGtNew <- function(nd, data, fm, tf, model, nknot, p.value = 1) {
    if (is.null(nd) || is.null(nd$split)) {
      return(0)
    }

    grp <- as.numeric(data[, nd$split$varid] > nd$split$breaks)
    if (model != "avg" && substr(model, 1, 4) != "last") {
      D.traj <-
        data[, sapply(tf, function(x) {
          return(paste0(all.vars(fm)[1], x))
        })]
      if (model %in% c("multest.acat", "multest.bonferroni", "mmreg", "mdistance") && nknot < 0) {
        D.traj <- spline.preprocess(D.traj, tf, nknot = nknot, alpha = 2)
        tf <- 1:nknot
      }
      score <- do.call(model,
                list(
                  D.traj = D.traj,
                  D.cov = cbind(data[, c("id", "treat")], grp),
                  tf = tf,
                  nknot = nknot,
                  p.value = p.value
                ))
      if (is.na(score)) {
        return(0)
      }
      if (model %in% c("bsgee", "bslme", "mmreg")) {
        score <- Chi2Normalize(score, nknot)
      }
    } else {
      score <- it.scalar.score(data, fm, grp, p.value)
    }

    left <- GetGtNew(
      nd = nd$kids[[1]],
      data = data[grp == 0, ],
      fm = fm,
      tf = tf,
      nknot = nknot,
      model = model,
      p.value = p.value
    )
    right <- GetGtNew(
      nd = nd$kids[[2]],
      data = data[grp == 1, ],
      fm = fm,
      tf = tf,
      nknot = nknot,
      model = model,
      p.value = p.value
    )
    score + ifelse(is.na(left), 0, left) + ifelse(is.na(right), 0, right)
  }


GetGtNewWrap <-
  function(nd,
           data,
           fm,
           tf,
           model,
           nknot,
           p.value = 1,
           preprocess = TRUE) {
    # get gt from new dataset
    if (is.null(nd) || is.null(nd$split)) {
      return(0)
    }

    if (preprocess) {
      pprocess <- preprocess(data, fm, tf, nknot, model)
      data = pprocess$data
      fm = pprocess$fm
      tf = pprocess$tf
    }

    GetGtNew(
      nd = nd,
      data = data,
      fm = fm,
      tf = tf,
      model = model,
      nknot = nknot,
      p.value = p.value
    )
  }


DeleteKids <- function(nd, id) {
  if (nd$id == id || is.terminal(nd)) {
    nd$kids <- nd$surrogates <- nd$split <- nd$info <- NULL
    return(nd)
  } else {
    nd$kids[[1]] <- DeleteKids(nd$kids[[1]], id)
    nd$kids[[2]] <- DeleteKids(nd$kids[[2]], id)
    return(nd)
  }
}


AlphaCalculate <- function(nd, nknot.transform) {
  if (is.null(nd) || is.null(nd$split)) {
    return(list())
  }
  left <- AlphaCalculate(nd$kids[[1]], nknot.transform)
  right <- AlphaCalculate(nd$kids[[2]], nknot.transform)
  record <- list(list(node = nd, alpha = GetGt(nd, nknot.transform) / NodeCount(nd)))
  if (length(left) > 0) {
    record <- c(record, left)
  }
  if (length(right) > 0) {
    record <- c(record, right)
  }
  return(record)
}


SplinePrune <- function(tr, nknot.transform) {
  # create a sequence of pruned nodes
  nd <- tr$node
  record <- list(list(node = nd, alpha = GetGt(nd, nknot.transform) / NodeCount(nd)))

  while (depth(nd) > 1) {
    seq <- AlphaCalculate(nd, nknot.transform)
    weak <-
      seq[[which.min(lapply(seq, function(x) {
        return(x$alpha)
      }))]]
    nd <- DeleteKids(nd, weak$node$id)
    if (is.null(nd$split$breaks)) {
      break
    }
    record <- c(record, list(list(
      node = nd, alpha = weak$alpha
    )))
  }
  record
}


OptimalSelf <- function(prune.seq, lambda, nknot.transform) {
  targets <- lapply(prune.seq, function(x) {
    return(GetGt(x$node, nknot.transform) - lambda * NodeCount((x$node)))
  })
  prune.seq[which.max(targets)]
}


# entry point
BtPrune <-
  function(data,
           tr,
           fm,
           tf,
           split.covs,
           nknot,
           maxdepth,
           model,
           p.value,
           nCutpoints = 20,
           alpha = 2,  # spline penalization
           lambdas = c(0, 2, 3, 4, log(nrow(data)), log(nrow(data) * length(tf))),
           B = 20) {
    pp <- preprocess(data, fm, tf, nknot, model)
    data.p <- pp$data
    fm.p <- pp$fm
    tf.p <- tf

    # knot selection for spline methods
    if (model %in% c("multest.acat", "multest.bonferroni", "mmreg", "mdistance")) {
      D.traj <-
        data[, sapply(tf, function(x) {
          return(paste0(all.vars(fm)[1], x))
        })]
      if (nknot <= 0) {
        nknot <- knot.select(D.traj, tf, alpha = alpha)
      }
    } else if (model == "avg" || substr(model, 1, 4) == "last") {
      nknot <- -1
    }

    nknot.transform <- ifelse(model %in% c("bsgee", "bslme", "mmreg"), nknot, -1)
    prune.seq <- SplinePrune(tr, nknot.transform)

    if (length(prune.seq) == 1) {
      return(tr)
    }
    alpha <- unlist(lapply(prune.seq, function(x) {
      return(x$alpha)
    }))
    alpha2 <-
      c(sqrt(alpha[1:(length(alpha) - 1)] * alpha[2:length(alpha)]), alpha[length(alpha)])
    bias <- matrix(0, B, length(prune.seq))
    for (ii in 1:B) {
      data.b <- data[sample(nrow(data), round(2 / 3 * nrow(data))), ]
      data.b.test <- data[!data$id %in% data.b$id, ]
      tr.b <- ItltTree(
        data = data.b,
        fm = fm,
        tf = tf,
        split.covs = split.covs,
        nknot = nknot,
        maxdepth = maxdepth,
        nCutpoints = nCutpoints,
        model = model,
        details = FALSE,
        p.value = p.value
      )
      prune.seq.b <- SplinePrune(tr.b, nknot.transform)
      for (jj in 1:length(alpha2)) {
        try({
          nd.b.jj <- OptimalSelf(prune.seq.b, alpha2[jj], nknot.transform)[[1]]$node
          # bias[ii, jj] <- GetGtNewWrap(nd.b.jj, data, model) - GetGtNewWrap(nd.b.jj, data.b, model)
          bias[ii, jj] <-
            GetGtNewWrap(
              nd = nd.b.jj,
              data = preprocess(data.b.test, fm, tf, nknot, model)$data,
              fm = fm.p,
              tf = tf.p,
              nknot = nknot,
              model = model,
              p.value = p.value,
              preprocess = FALSE
            ) - GetGt(nd.b.jj, nknot.transform)
        })
      }
    }

    gt <-
      unlist(lapply(prune.seq, function(x) {
        return(
          GetGtNewWrap(
            nd = x$node,
            data = data.p,
            fm = fm.p,
            tf = tf.p,
            nknot = nknot,
            model = model,
            p.value = p.value,
            preprocess = FALSE
          )
        )
      })) + colMeans(bias, na.rm = TRUE)
    penalty <-
      unlist(lapply(prune.seq, function(x) {
        return(NodeCount(x$node))
      }))

    ret <- list(gt = gt, penalty = penalty)
    for (lambda in lambdas) {
      node <- prune.seq[[which.max(gt - penalty * lambda)]]$node
      fitted <- fitted_node(node, data = data)
      ret[[toString(lambda)]] <- list(party(
        node,
        data = data,
        fitted = data.frame("(fitted)" = fitted,
                            check.names = FALSE),
        terms = NULL
      ))
    }

    ret
  }
