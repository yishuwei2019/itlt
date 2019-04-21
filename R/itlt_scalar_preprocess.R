preprocess <- function(data, fm, tf, nknot, model) {
  if (!"id" %in% colnames(data)) {
    colnames(data)[1] <- "id"
  }
  colnames(data) <- gsub(all.vars(fm)[2], "treat", colnames(data))
  y.names <-
    sapply(tf, function(x) {
      return(paste0(all.vars(fm)[1], x))
    })
  D.traj <- data[, y.names]

  if (model == "avg") {
    scaler <- rowMeans(D.traj, na.rm = TRUE)
    data <- cbind(data, scaler)
    fm <- update(fm, scaler ~ .)
  } else if (substr(model, 1, 4) == "last") {
    last <- strtoi(substr(model, 5, nchar(model)))
    if (last == 1) {
      scaler <- D.traj[, ncol(D.traj)]
    } else {
      scaler <- rowMeans(D.traj[, (ncol(D.traj) - last + 1):ncol(D.traj)], na.rm = TRUE)
    }
    data <- cbind(data, scaler)
    fm <- update(fm, scaler ~ .)
  }

  list(data = data, fm = fm)
}
