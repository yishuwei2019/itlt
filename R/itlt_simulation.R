#' Simulation for longitudinal datasets
#' @param N sample size
#' @param tf time frame
#' @param e error standard deviation
#' @param p missing value probability
#' @param asigma sd of subject specific random effect
#' @param seed seed for random number
#' @param type simulation setting (1: non-linearity, 2: periodicity,
#' 3: non-linearity, 4: no heteorogeneity in treatment effect)
#' @examples
#' ItltSimulation1(250, 1:12)
#' @return a simulated longitudinal dataset

ItltSimulation1 <-
  function(N,
           tf,
           e = .05,
           p = 0,
           asigma = .1,
           seed = "1234",
           type = 1) {
    set.seed(seed)
    if (type == 1) {
      mu10 <- .55 + .45 * exp(-tf / 2) # treat0
      mu00 <- exp(-tf / 16) # control 0
      mu11 <- exp(-tf / 8) # treat 1
      mu01 <- .3 + .7 * exp(-tf / 2) # control 1
    } else if (type == 2) {
      mu10 <-
        .7 + .3 * cos(tf * 1.2) - 0.03 * tf #treat 0
      mu10 <- sapply(mu10, function(x) {
        return(max(x, 0))
      })
      mu00 <- .7 + .3 * cos(tf * .9) # control 0
      mu11 <- exp(-tf / 8) # treat 1
      mu01 <- .25 + .75 * exp(-tf / 2) # control 1
    } else if (type == 3) {
      mu10 <- .5 + .5 * exp(-tf) # treatment 0
      mu00 <- exp(-tf / 8) # control 0
      mu11 <- exp(-tf / 4) # treatment 1
      mu01 <- .25 + .75 * exp(-tf) # control 1
    } else if (type == 4) {
      # two groups with same trajectory (group is meaning less)
      mu10 <- .5 + .5 * exp(-tf) # treatment 0
      mu00 <- exp(-tf / 8) # control 0
      mu11 <- .5 + .5 * exp(-tf) # treatment 1
      mu01 <- exp(-tf / 8) # control 1
    }

    id <- 1:N
    treatment <-
      c(rep(0, N / 2), rep(1, N / 2)) # half control half treatment
    X1 <- rbinom(N, 1, 0.5)
    X2 <- rbinom(N, 1, 0.5)
    X3 <- rbinom(N, 1, 0.5)
    X4 <- runif(N, 0, 1)
    X5 <- runif(N, 0, 1)
    X6 <- runif(N, 0, 1)
    X7 <- rbinom(N, 1, 0.5)
    X8 <- rbinom(N, 1, 0.5)
    X9 <- rbinom(N, 1, 0.5)
    X10 <- rbinom(N, 1, 0.5)

    group <- X4 > .7

    covs <- cbind(treatment, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, group)

    ymatrix <- matrix(0, N, length(tf))
    ymatrix[treatment == 1 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 0), mu10, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 0), mu00, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 1), mu11, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 1), mu01, simplify = FALSE))
    ymatrix <-
      ymatrix + matrix(rnorm(N * length(tf), 0, e), nrow = N) +
      matrix(rep(rnorm(N, 0, asigma), each = length(tf)),
             nrow = N,
             byrow = TRUE)

    ymatrix <-
      t(apply(ymatrix, 1, function(x) {
        ind <-
          c(0, rbinom(length(x) - 1, 1, p))
        x[ind == 1] <- NA
        return(x)
      }))
    colnames(ymatrix) <- paste0("y", tf)

    data <- data.frame(cbind(covs, ymatrix))
    # data$avg.y <- rowMeans(ymatrix, na.rm = TRUE)
    data$avg.y <-
      rowMeans(ymatrix[, (ncol(ymatrix) - 3):ncol(ymatrix)], na.rm = TRUE) - ymatrix[, 1]
    data <- data[sample(nrow(data)), ]
    data <- cbind(id, data, group)
    data
  }


# sim.depth1.data2 <-
#   function(N,
#            b0,
#            b1,
#            tf,
#            e = .05,
#            seed = "1234",
#            type = 2) {
#     # interaction with slope
#     # b0: intercept, treat1, control1, treat0, control0
#     # b1: slope
#     # e: error sigma
#
#     set.seed(seed)
#     # N: sample
#     design <-
#       do.call(rbind, replicate(N, rep(1, length(tf)), simplify = FALSE))
#     design1 <- do.call(rbind, replicate(N, tf, simplify = FALSE))
#     design2 <-
#       do.call(rbind, replicate(N, tf ^ 2, simplify = FALSE))
#
#     id <- 1:N
#     treatment <-
#       c(rep(0, N / 2), rep(1, N / 2)) # half control half treatment
#     X1 <- rbinom(N, 1, 0.5)
#     X2 <- rbinom(N, 1, 0.5)
#     X3 <- rbinom(N, 1, 0.5)
#     X4 <- runif(N, 0, 1)
#     X5 <- runif(N, 0, 1)
#     X6 <- runif(N, 0, 1)
#     X7 <- runif(N, 0, 1)
#     X8 <- runif(N, 0, 1)
#     X9 <- rbinom(N, 1, 0.5)
#     X10 <- rbinom(N, 1, 0.5)
#     group <- X1 == 1 # trajectory pattern differ for each group
#     covs <- cbind(treatment, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, group)
#
#     ymatrix <- matrix(0, N, length(tf))
#     if (type == 2) {
#       ind1 <- treatment == 0 & X4 <= .7
#       ind2 <- treatment == 1 & X4 <= .7
#       ind3 <- treatment == 0 & X4 > .7
#       ind4 <- treatment == 1 & X4 > .7
#       ymatrix[ind1, ] <- 1 - 0.0 * design1[ind1, ] - 0.002 * design2[ind1, ] * covs[ind1, "X4"]
#       ymatrix[ind2, ] <- 1 - 0.0 * design1[ind2, ] - 0.008 * design2[ind2, ] * covs[ind2, "X4"]
#       ymatrix[ind3, ] <- .8 - 0.0 * design1[ind3, ] - 0.002 * design2[ind3, ] * covs[ind3, "X4"]
#       ymatrix[ind4, ] <- .8 - 0.0 * design1[ind4, ] - 0.004 * design2[ind4, ] * covs[ind4, "X4"]
#     }
#     ymatrix <-
#       ymatrix + matrix(rnorm(N * length(tf), 0, e), nrow = N)
#     colnames(ymatrix) <-
#       sapply(tf, function(x) {
#         return(paste0("y", x))
#       })
#
#     data <- cbind(covs, ymatrix)
#     data <- data.frame(data)
#     # data$avg.y <- rowMeans(ymatrix)
#     data$avg.y <-
#       rowMeans(ymatrix[, (ncol(ymatrix) - 3):ncol(ymatrix)], na.rm = TRUE) - ymatrix[, 1]
#     data <- data[sample(nrow(data)), ]
#     data <- cbind(id, data, group)
#
#     data
#   }


sim.depth2.data1 <-
  function(N,
           tf,
           e = .05,
           p = 0,
           asigma = .1,
           seed = "1234",
           type = 1) {
    set.seed(seed)
    mu10 <- exp(-tf / 4) # treatment 0
    mu00 <- .25 + .75 * exp(-tf) # control 0
    mu11 <- .5 + .5 * exp(-tf) # treatment 1
    mu01 <- exp(-tf / 8) # control 1
    mu12 <- .5 + .5 * exp(-tf / 2) # treatment 2
    mu02 <- exp(-tf / 8) # control 2
    mu13 <- exp(-tf / 2)  # treatment 3
    mu03 <- .25 + .75 * exp(-tf)  # control 3

    id <- 1:N
    treatment <-
      c(rep(0, N / 2), rep(1, N / 2)) # half control half treatment
    X1 <- rbinom(N, 1, 0.5)
    X2 <- rbinom(N, 1, 0.5)
    X3 <- rbinom(N, 1, 0.5)
    X4 <- runif(N, 0, 1)
    X5 <- runif(N, 0, 1)
    X6 <- runif(N, 0, 1)
    X7 <- rbinom(N, 1, 0.5)
    X8 <- rbinom(N, 1, 0.5)
    X9 <- rbinom(N, 1, 0.5)
    X10 <- rbinom(N, 1, 0.5)
    group <- rep(0, N)  # X1 == 1 & X4 <= .7
    group[X1 == 1 & X4 > .7] <- 1
    group[X1 == 0 & X4 > .6] <- 2
    group[X1 == 0 & X4 <= .6] <- 3
    covs <- cbind(treatment, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, group)

    ymatrix <- matrix(0, N, length(tf))
    ymatrix[treatment == 1 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 0), mu10, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 0), mu00, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 1), mu11, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 1), mu01, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 2, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 2), mu12, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 2, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 2), mu02, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 3, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 3), mu13, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 3, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 3), mu03, simplify = FALSE))
    ymatrix <-
      ymatrix + matrix(rnorm(N * length(tf), 0, e), nrow = N) +
      matrix(rep(rnorm(N, 0, asigma), each = length(tf)),
             nrow = N,
             byrow = TRUE)
    ymatrix <-
      t(apply(ymatrix, 1, function(x) {
        ind <-
          c(0, rbinom(length(x) - 1, 1, p))

        x[ind == 1] <- NA
        return(x)
      }))
    colnames(ymatrix) <- paste0("y", tf)

    data <- data.frame(cbind(covs, ymatrix))
    data <- data[sample(nrow(data)), ]
    data <- cbind(id, data, group)
    data
  }


sim.depth2.data2 <-
  function(N,
           tf,
           e = .05,
           p = 0,
           asigma = .1,
           seed = "1234",
           type = 1) {
    # N: sample
    # p: missing value
    # include setting 1, 3, 4, 4_re(5)
    set.seed(seed)
    mu10 <-
      .7 + .3 * cos(tf * 1.2) - 0.03 * tf #treat 0
    mu10 <- sapply(mu10, function(x) {
      return(max(x, 0))
    })
    mu00 <- .7 + .3 * cos(tf * .9) # control 0
    mu11 <- exp(-tf / 8) # treat 1
    mu01 <- .25 + .75 * exp(-tf / 2) # control 1
    mu12 <- .5 + .5 * exp(-tf / 2) # treatment 2
    mu02 <- exp(-tf / 8) # control 2
    mu13 <- exp(-tf / 2)  # treatment 3
    mu03 <- .25 + .75 * exp(-tf)  # control 3

    id <- 1:N
    treatment <-
      c(rep(0, N / 2), rep(1, N / 2)) # half control half treatment
    X1 <- rbinom(N, 1, 0.5)
    X2 <- rbinom(N, 1, 0.5)
    X3 <- rbinom(N, 1, 0.5)
    X4 <- runif(N, 0, 1)
    X5 <- runif(N, 0, 1)
    X6 <- runif(N, 0, 1)
    X7 <- rbinom(N, 1, 0.5)
    X8 <- rbinom(N, 1, 0.5)
    X9 <- rbinom(N, 1, 0.5)
    X10 <- rbinom(N, 1, 0.5)
    group <- rep(0, N)  # X4 > .7 & X1 == 1
    group[X4 > .7 & X1 == 0] <- 1
    group[X4 <= .7 & X2 == 1] <- 2
    group[X4 <= .7 & X2 == 0] <- 3
    covs <- cbind(treatment, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, group)

    ymatrix <- matrix(0, N, length(tf))
    ymatrix[treatment == 1 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 0), mu10, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 0, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 0), mu00, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 1), mu11, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 1, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 1), mu01, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 2, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 2), mu12, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 2, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 2), mu02, simplify = FALSE))
    ymatrix[treatment == 1 &
              group == 3, ] <-
      do.call(rbind, replicate(sum(treatment == 1 &
                                     group == 3), mu13, simplify = FALSE))
    ymatrix[treatment == 0 &
              group == 3, ] <-
      do.call(rbind, replicate(sum(treatment == 0 &
                                     group == 3), mu03, simplify = FALSE))
    ymatrix <-
      ymatrix + matrix(rnorm(N * length(tf), 0, e), nrow = N) +
      matrix(rep(rnorm(N, 0, asigma), each = length(tf)),
             nrow = N,
             byrow = TRUE)
    ymatrix <-
      t(apply(ymatrix, 1, function(x) {
        ind <-
          c(0, rbinom(length(x) - 1, 1, p))

        x[ind == 1] <- NA
        return(x)
      }))
    colnames(ymatrix) <- paste0("y", tf)

    data <- data.frame(cbind(covs, ymatrix))
    data <- data[sample(nrow(data)), ]
    data <- cbind(id, data, group)
    data
  }
