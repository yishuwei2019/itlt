#' Implement tree structured subgroup identification on longitudinal datasets
#'
#' Methods implemented include: mixed effects based method (LRT based method),
#' GEE based method, multiple test based method, multivariate multiple regression
#' (mmp) based method. mixed effects and GEE are for regular longitudinal datasets and
#' multiple test based and mmp based method are for intensive longitudinal datasets.
#' This function is based on partynode package. Also interaction trees by scalar
#' based methods (avg, last-k) are also implemented.But the scalar based method
#' needs to have the same data format (see data, fm, tf for details), otherwise
#' regular code for interaction tree could be used
#'
#'
#' @param data dataset in wide format, please see fm and time frame for column
#' name restrictions. Also data must have a column named "id", otherwise the first
#' column will be specified to be "id"
#' @param fm formula of the form "response name ~ treatment name". data must
#' contain column names of the response with all its time frames. For instance,
#' if fm = "y~treatment", tf = 1:12, then columns names in the dataset must
#' include y1, y2,..., y12. And these longitudinal observations will be used
#' as reponse for subgroup identfication
#' @param tf time frame
#' @param split.covs splitting covariates, vector of strings
#' @param nknot knot numbers, for multiple test based and mmp based methods,
#' nknot <= 0 means automatic selection of knot numbers using GCV
#' @param p.value significance level for a node to be considered, default is 1
#' @param delta lower bound of the ratio of the size of the child node and the
#' parent node, default is 0.1
#' @param nCutpoints number of cutoff points for a continuous splitting covariate,
#' default is 20
#' @param minsplit minimum number of samples for a node to be considered for
#' splitting, default is 30
#' @param size.limit minimum size for a node, default is delta * minsplit
#' @param maxdepth maximum depth of the tree, default is 4
#' @param details whether details in the treee building process will be printed,
#' default is FALSE
#' @param alpha tuning parametering in the selection of knot numbers using
#' smoothing splines, specific for multiple test based method and MMP based method,
#' default is 2
#' @param model method to be implemented.
#' @return identified tree of partynode type
ItltTree <-
  function(data,  # wide format and has the format of "y.name1, y.name2,....
           fm,  # fm is of the format y.name ~ treatment.name + base covariates
           tf,
           split.covs,  # could be both numerical indicator or variable names
           nknot = 5,
           p.value = 1,
           delta = 0.1,
           nCutpoints = 20,
           minsplit = 20,
           size.limit = 10,
           maxdepth = 3,
           details = "TRUE",
           alpha = 2.5,
           model = "bslme") {
    # data must have a column named "id"
    if (is.numeric(split.covs[1])) {
      split.covs <-
        sapply(split.covs, function(x) {
          return(colnames(data)[x])
        })
    }

    if (model == "bsgee") {
      size.limit = 50
    }

    # deal with scalar methods
    pprocess <- preprocess(data, fm, tf, nknot, model)
    data <- pprocess$data
    fm <- pprocess$fm

    # knot selection for spline methods
    if (model %in% c("multest.acat", "multest.bonferroni", "mmreg", "mdistance")) {
      D.traj <-
        data[, sapply(tf, function(x) {
          return(paste0(all.vars(fm)[1], x))
        })]
      if (nknot <= 0) {
        nknot <- knot.select(D.traj, tf, alpha = alpha)
      }
    }

    node <-
      ItltGrow(
        data = data,
        fm = fm,
        tf = tf,
        split.covs = split.covs,
        nknot = nknot,
        p.value = p.value,
        name = "1",
        depth = 1,
        delta = delta,
        nCutpoints = nCutpoints,
        minsplit = minsplit,
        size.limit = size.limit,
        maxdepth = maxdepth,
        details = details,
        alpha = alpha,
        model = model
      )

    fitted <- fitted_node(node, data = data)
    ret <- party(
      node,
      data = data,
      fitted = data.frame("(fitted)" = fitted,
                          check.names = FALSE),
      terms = NULL
    )

    ret
    # as.constparty(ret)
  }
