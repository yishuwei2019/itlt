#' Selection of number of knots in multiple test base and mmp based methods
#'
#' For intensive longitudinal dataset (multiple test based method and mmp based
#' method), a preprocessing step using smoothing splines is used.Generalized
#' cross validation (GCV) is used for knot number selection.
#'
#' @param D.traj dataset consists purely of longitudinal observations
#' @param tf time frame
#' @param alpha penalization term in GCV
#'
#' @return mean and standard deviation of GCV of each possible knot number
knot.select <- function(D.traj, tf, alpha = 2) {
  gcvs <- matrix(Inf, 2, ncol(D.traj) / 2)
  knot_min <- 3
  for (kk in knot_min:(ncol(D.traj) / 2)) {
    crits <- apply(D.traj, 1, function(row) {
      return(smooth.spline(
        x = tf[!is.na(row)],
        y = row[!is.na(row)],
        nknots = kk,
        cv = FALSE,
        lambda = alpha
      )$cv.crit)
    })
    gcvs[1, kk] <- mean(crits)
    gcvs[2, kk] <- sd(crits)
  }
  # nknot <- which.min(round(gcvs, 5))
  return(which.min(gcvs[1, ]))
}
