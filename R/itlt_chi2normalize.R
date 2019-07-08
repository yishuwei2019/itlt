#' Nomalizing a \chi^2 square statistic from df(n) to df(1)
#'
#' The pruning procedure is proposed under the assumption that test statistics at
#' each node are \chi^2(1) distributed.
#'
#' @param t original test statistic
#' @param n original degree of freedom
#' @examples
#' Chi2Normalize(3, 5)
#' @return \chi^2 test statistic of df(1) with the same significance level
Chi2Normalize <- function(t, n) {
  if (n == 1) {
    return(t)
  }

  w1 = (sqrt(2 * t) - sqrt(2 * n - 1) + 1) ^ 2 / 2
  w2 = max(0, (7 / 9 + sqrt(n) * ((t / n) ^ (1 / 3) - 1 + 2 / 9 / n)) ^ 3)

  if (t < n + 10 * sqrt(2 * n)) {
    return(w2)
  } else if ((t >= n + 10 * sqrt(2 * n)) && w2 < t) {
    return((w1 + w2) / 2)
  } else {
    return(w1)
  }
}
