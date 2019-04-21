ItltPartition <- function(
           data,
           fm,
           tf,
           split.covs,
           nknot,
           p.value,
           delta,
           nCutpoints,
           size.limit,
           details,
           alpha,
           model) {
    # categorical variable must be 0/1

    if (model != "avg" && substr(model, 1, 4) != "last") {
      D.traj <- data[, paste0(all.vars(fm)[1], tf)]
      if (model %in% c("multest.acat", "multest.bonferroni", "mmreg", "mdistance")) {
        D.traj <- spline.preprocess(D.traj, tf, nknot = nknot, alpha = alpha)
        tf <- 1:nknot
      }
    }

    best.var <- best.cut <- NA
    max.score <- -Inf
    if (details) {
      print(paste("the sample size:", nrow(data),
        "and number of treated:", sum(data[, "treat"])))
    }

    # variables having less than five values are considered categorical
    cat.list <-
      apply(data[, split.covs], 2, function(x) {
        return(length(unique(x)))
      })
    cat.list <- names(cat.list[cat.list <= 5])

    for (ii in split.covs) {
      z <- data[, ii]
      if (length(unique(z)) == 1) {
        next
      }
      if (is.element(ii, cat.list)) {
        zcuts <- sort(unique(z))
        zcuts <- zcuts[1:(length(zcuts) - 1)]
      } else {
        zcuts <-
          quantile(z,
                   probs = seq(delta, 1 - delta, length = nCutpoints),
                   na.rm = TRUE)
        zcuts <- unique(zcuts)
      }

      for (jj in zcuts) {
        if (details) {
          print(paste("cov:", ii, "cut:", jj))
        }
        score <- NA
        if (is.element(ii, cat.list)) {
          # grp <- sign(z %in% jj)
          grp <- sign(z <= jj)
          cut1 <- paste(jj, collapse = " ")
        } else {
          grp <-
            sign(z <= jj)  # 2 - kidids_split(split, data) = grp   because TRUE == 1
          cut1 <- as.character(jj)
        }
        grp <- matrix(grp, dimnames = list(NULL, "grp"))

        if (sum(grp) >= size.limit &&
            sum(grp) <= length(grp) - size.limit) {
          if (model != "avg" && substr(model, 1, 4) != "last") {
            score <-
              do.call(
                model,
                list(
                  D.traj = D.traj,
                  D.cov = cbind(data[, c("id", "treat")], grp),
                  tf = tf,
                  nknot = nknot,
                  p.value = p.value
                )
              )
          } else {
            score <- it.scalar.score(data, fm, grp, p.value)
          }


          if (!is.na(score) && score > max.score) {
            max.score <- score
            best.var <- ii
            best.cut <- jj
          }

          if (details) {
            print(cbind(
              best.var = best.var,
              best.cut = best.cut,
              score = score,
              max.score = max.score
            ))
          }
        }
      }
    }

    if (is.infinite(max.score)) {
      return(NULL)
    }

    # this is for pruning process
    if (model == "bsgee" || model == "bslme") {
      max.score <- Chi2Nomalize(max.score, nknot)
    }

    partysplit(
      varid = which(names(data) == best.var),
      breaks = best.cut,
      index = 1:2,
      info = list(
        size = nrow(data),
        best.var = best.var,
        score = max.score
      )
    )
  }
