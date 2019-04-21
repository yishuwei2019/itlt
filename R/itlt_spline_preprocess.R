#' Use smoothing spline to preprocess original longitudinal data
#'
#' For multiple test based and MMP based methods, preprocess using smoothing
#' splines are the first step
#'
#' @param D.traj dataset consists purely of original longitudinal observations
#' @param tf time frame
#' @param nknot knot number
#' @param alpha penalization term in GCV
#'
#' @return new longitudinal dataset with fitted spline coefficients
spline.preprocess <- function(D.traj, tf, nknot, alpha) {
  D.traj <- t(apply(D.traj, 1, function(row) {
    tt <- smooth.spline(
        x = tf[!is.na(row)],
        y = row[!is.na(row)],
        nknots = nknot,
        lambda = alpha,
        cv = FALSE
      )
    return(tt$fit$coef)
  }))

  D.traj <- D.traj[, 3:ncol(D.traj)]
  colnames(D.traj) <- paste0("y", 1:ncol(D.traj))
  return(D.traj)
}
