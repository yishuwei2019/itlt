#' Calculate test statistic using mixed effects model (LRT based model)
#'
#' Mixed effects based model tests trajectories by testing model fitting.
#' Trajectories are assumed to follow mixed effects models. Under the null
#' hypothesis there are only one treatment effect curve for both nodes. Under
#' the alternative, separate treatment effect curves are assumed for both nodes
#'
#' @param D.traj dataset of purely longitudinal observations in wide format
#' @param D.cov covariate matrix of intercept, treat, grp, treat * grp
#' @param tf time frame
#' @param nknot knot number
#' @param p.value significance value required for a node to be splitted, default
#' is 1
#' @param degree degree of spline coefficients, default is 3
#'
#' @return (LRT) chi^2 test statistic of df(nknot) under the mixed effects model
bslme <- function(
  D.traj,
  D.cov,
  tf,
  nknot = 5,
  p.value = .05,
  degree = 3) {
  # mixed effects based method (LRT based method)
  assertthat::are_equal(ncol(D.traj), length(tf))
  assertthat::are_equal(nrow(D.traj), nrow(D.cov))
  covs <- colnames(D.cov)[3:ncol(D.cov)]
  if (!"id" %in% covs) {
    covs <- c("id", covs)
  }

  N.t <- length(tf)
  # there are inner.knot - 1 intervals inside time frame
  # gap <- (tf[N.t] - tf[1]) / (nknot - 1)
  # knots <-
  #   c(seq(tf[1] - degree * gap, tf[N.t] + degree * gap, by = gap))
  # B <- splineDesign(knots, tf, degree + 1)
  # B <- B[, 2:(ncol(B) - 1)]
  # B0 <- matrix(0, N.t, ncol(B))
  # new basis system
  gap <- (tf[N.t] - tf[1]) / (nknot - 1)
  knots <- seq(tf[1], tf[N.t], by = gap)
  B <- bs(tf, knots = knots)
  B <-  B[, 2:(ncol(B) - 1)]
  B0 <- matrix(0, N.t, ncol(B))


  # organize data by group indicator
  g <- rep(0, nrow(D.cov)) # 0 is left control(grp == 1, treat == 0)
  g[D.cov$grp == 1 & D.cov$treat == 1] <- 1 # left treat
  g[D.cov$grp == 0 & D.cov$treat == 0] <- 2 # right control
  g[D.cov$grp == 0 & D.cov$treat == 1] <- 3 # right treat
  data <- cbind(D.cov, D.traj, g)
  data <- data[order(data$g, decreasing = F),]

  # organize data to long format
  y <- c(t(data[, colnames(D.traj)]))
  D.cov <-
    apply(data[, covs], 2, function(x) {
      return(rep(x, each = N.t))
    })
  # design matrix under null
  V.lc <-
    do.call(rbind, replicate(sum(g == 0), cbind(B, B0, B0), simplify = F))
  V.lt <-
    do.call(rbind, replicate(sum(g == 1), cbind(B, B, B0), simplify = F))
  V.rc <-
    do.call(rbind, replicate(sum(g == 2), cbind(B0, B0, B), simplify = F))
  V.rt <-
    do.call(rbind, replicate(sum(g == 3), cbind(B0, B, B), simplify = F))
  V <- rbind(V.lc, V.lt, V.rc, V.rt)
  colnames(V) <- paste0("v", 1:ncol(V))

  # design matrix under alternative
  U.lc <-
    do.call(rbind, replicate(sum(g == 0), cbind(B, B0, B0, B0), simplify = F))
  U.lt <-
    do.call(rbind, replicate(sum(g == 1), cbind(B, B, B0, B0), simplify = F))
  U.rc <-
    do.call(rbind, replicate(sum(g == 2), cbind(B0, B0, B, B0), simplify = F))
  U.rt <-
    do.call(rbind, replicate(sum(g == 3), cbind(B0, B, B, B), simplify = F))
  U <- rbind(U.lc, U.lt, U.rc, U.rt)
  colnames(U) <- paste0("u", 1:ncol(U))

  D.null <- cbind(D.cov, y, V)
  D.alt <- cbind(D.cov, y, U)
  D.null <- D.null[rowSums(is.na(D.null)) == 0,]
  D.alt <- D.alt[rowSums(is.na(D.alt)) == 0,]
  fm.null <-
    as.formula(paste0("y~", paste(colnames(V), collapse = "+")))
  fm.alt <-
    as.formula(paste0("y~", paste(colnames(U), collapse = "+")))

  ## mixed effect model
  fm.random = as.formula("random =~ 1|id")
  lrt <- NA
  try({
    null.lme <-
      lme(fm.null,
          data = data.frame(D.null),
          fm.random,
          na.action = na.omit)
    alt.lme <-
      lme(fm.alt,
          data = data.frame(D.alt),
          fm.random,
          na.action = na.omit)
    lrt <- -2 * (null.lme$logLik - alt.lme$logLik)
  })

  lrt <- max(lrt, 0)
  pv = 1 - pchisq(lrt, ncol(B))
  if (is.na(pv) || pv > p.value) {
    return(NA)
  } else {
    return(lrt)
  }
}


bsgee <- function(D.traj,
                  D.cov,
                  tf,
                  nknot,
                  p.value,
                  degree = 3) {
  assertthat::are_equal(ncol(D.traj), length(tf))
  assertthat::are_equal(nrow(D.traj), nrow(D.cov))
  N.t <- length(tf)
  # there are inner.knot - 1 intervals inside time frame
  gap <- (tf[N.t] - tf[1]) / (nknot - 1)
  knots <-
    c(seq(tf[1] - degree * gap, tf[N.t] + degree * gap, by = gap))
  B <- splineDesign(knots, tf, degree + 1)
  B <- B[, 2:(ncol(B) - 1)]
  B0 <- matrix(0, N.t, ncol(B))

  # organize data by group indicator
  g <- rep(0, nrow(D.cov)) # 0 is left control(grp == 1, treat == 0)
  g[D.cov$grp == 1 & D.cov$treat == 1] <- 1 # left treat
  g[D.cov$grp == 0 & D.cov$treat == 0] <- 2 # right control
  g[D.cov$grp == 0 & D.cov$treat == 1] <- 3 # right treat
  data <- cbind(D.cov, D.traj, g)
  data <- data[order(data$g, decreasing = F),]

  # organize data to long format
  y <- c(t(data[, colnames(D.traj)]))
  D.cov <-
    apply(data[, !(colnames(data) %in% colnames(D.traj)), drop = F],
          2, function(x) {
            return(rep(x, each = N.t))
          })
  U.lc <-
    do.call(rbind, replicate(sum(g == 0), cbind(B, B0, B0, B0), simplify = F))
  U.lt <-
    do.call(rbind, replicate(sum(g == 1), cbind(B, B, B0, B0), simplify = F))
  U.rc <-
    do.call(rbind, replicate(sum(g == 2), cbind(B0, B0, B, B0), simplify = F))
  U.rt <-
    do.call(rbind, replicate(sum(g == 3), cbind(B0, B, B, B), simplify = F))
  U <- rbind(U.lc, U.lt, U.rc, U.rt)
  colnames(U) <- paste0("u", 1:ncol(U))

  gee.fm <-
    as.formula(paste0("y~", paste(colnames(U), collapse = "+")))
  gee.data <- data.frame(cbind(D.cov, y, U))
  gee.data <- gee.data[rowSums(is.na(gee.data)) == 0,]
  W2 <- NA
  try({
    gee.fit <- gee(gee.fm, id = gee.data$id, data = gee.data, silent = TRUE)
    beta <- gee.fit$coefficients
    beta <- beta[(2 + 3 * ncol(B)):length(beta)]
    omega <- gee.fit$robust.variance
    omega <-
      omega[(2 + 3 * ncol(B)):nrow(omega), (2 + 3 * ncol(B)):nrow(omega)]
    W2 <- t(matrix(beta)) %*% solve(omega) %*% matrix(beta)
  })

  pv = 1 - pchisq(W2, ncol(B))
  if (is.na(pv) || pv > p.value) {
    return(NA)
  } else {
    return(W2)
  }
}


multest.acat <- function(D.traj, D.cov, tf, nknot, p.value) {
  # multiple test based method with ACAT

  assertthat::are_equal(ncol(D.traj), length(tf))
  assertthat::are_equal(nrow(D.traj), nrow(D.cov))
  ND <- cbind(D.traj, D.cov)
  # ND <- ND[rowSums(is.na(ND)) == 0,]

  t.acat <- p.acat <- NA
  try({
    ss <- sapply(1:nknot, function(x) {
      ff <- paste0(colnames(D.traj)[x],
                   "~ 1 + treat + grp + treat * grp")

      rr <- summary(lm(as.formula(ff), data = ND))$coefficients

      return(rr["treat:grp", "Pr(>|t|)"])
    })
    t.acat <- mean(tan((0.5 - ss) * pi))
    p.acat <- .5 - atan(t.acat) / pi
  })


  if (is.na(p.acat) || p.acat > p.value) {
    return(NA)
  } else {
    return(t.acat)
  }
}


multest.bonferroni <- function(D.traj, D.cov, tf, nknot, p.value) {
  # multiple test based method with Bonferroni correction
  # (default multiple test based method)

  assertthat::are_equal(ncol(D.traj), length(tf))
  assertthat::are_equal(nrow(D.traj), nrow(D.cov))
  ND <- cbind(D.traj, D.cov)
  # ND <- ND[rowSums(is.na(ND)) == 0,]

  tstat <- pv <- NA
  try({
    ss <- sapply(1:nknot, function(x) {
      ff <- paste0(colnames(D.traj)[x],
                   "~ 1 + treat + grp + treat * grp")

      rr <- summary(lm(as.formula(ff), data = ND))$coefficients

      return(rr["treat:grp", c("t value", "Pr(>|t|)")])
    })
    tstat <- max(ss["t value", ] ^ 2)
    pv <- min(ss["Pr(>|t|)", ])
  })


  if (is.na(pv) || pv > p.value) {
    return(NA)
  } else {
    # Boforroni correction of test statistic
    # pv <- min(pv * nknot, 1)
    # return(abs(qnorm(pv)))
    return(tstat)
  }
}

mmreg <- function(D.traj, D.cov, tf, nknot, p.value) {
  # multivariate multiple regression

  assertthat::are_equal(ncol(D.traj), length(tf))
  assertthat::are_equal(nrow(D.traj), nrow(D.cov))
  ND <- cbind(D.traj, D.cov)
  # ND <- ND[rowSums(is.na(ND)) == 0,]

  ss <- NA
  try({
    full.fm <- as.formula(paste0(
      "cbind(",
      paste(colnames(D.traj), collapse = ","),
      ")~ treat + grp + treat * grp"
    ))

    fit <- mvlm(full.fm,
                data = ND,
                contr.factor = c("treat:grp"))
    ss <- summary(fit)[[1]]["treat:grp",]
    ss <- ifelse(summary(fit)[[4]]["treat:grp",] <= p.value, ss ^ 2, NA)
  })

  ss
}


it.scalar.score <- function(data, fm, grp, p.value) {
  y.name <- all.vars(fm)[1]

  ND <- cbind(data[, c(y.name, "treat")], grp)
  colnames(ND) <- c("y", "treat", "grp")

  t.stat <- NA
  p_value <- NA
  reg <- lm(y ~ 1 + treat + grp + treat * grp, data = ND)
  try({
    t.stat <- summary(reg)$coefficients[4, 3]
    p_value <- summary(reg)$coefficients[4, 4]
  })

  ifelse(p_value <= p.value, t.stat ^ 2, NA)
}
