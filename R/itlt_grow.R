ItltGrow <-
  function(data,
           fm,
           tf,
           split.covs,
           nknot,
           p.value,
           name = "1",
           depth = 1,
           delta = 0.2,
           nCutpoints = 20,
           minsplit = 50,
           size.limit = 30,
           maxdepth = 3,
           details = FALSE,
           alpha = 2,
           model = "bslme") {
    if (depth > maxdepth || length(unique(data[, "id"])) < minsplit) {
      return(partynode(id = strtoi(gsub("2", "0", name), base = 2)))
    }

    if (details) {
      print(paste0("node: ", name))
    }
    split <-
      ItltPartition(
        data,
        fm,
        tf = tf,
        split.covs = split.covs,
        nknot = nknot,
        p.value = p.value,
        delta = delta,
        nCutpoints = nCutpoints,
        size.limit = size.limit,
        details = details,
        alpha = alpha,
        model = model
      )

    if (is.null(split)) {
      return(partynode(id = strtoi(gsub("2", "0", name), base = 2)))
    }

    kidids <- kidids_split(split, data)
    kids <- vector(mode = "list", 2)

    kids[[1]] <-
      ItltGrow(
        data = data[kidids == 1,],
        fm = fm,
        tf = tf,
        split.covs = split.covs,
        nknot = nknot,
        p.value = p.value,
        name = paste0(name, "1"),
        depth = depth + 1,
        delta = delta,
        nCutpoints = nCutpoints,
        minsplit = minsplit,
        size.limit = size.limit,
        maxdepth = maxdepth,
        details = details,
        alpha = alpha,
        model = model
      )

    kids[[2]] <-
      ItltGrow(
        data = data[kidids == 2,],
        fm = fm,
        tf = tf,
        split.covs = split.covs,
        nknot = nknot,
        p.value = p.value,
        name = paste0(name, "2"),
        depth = depth + 1,
        delta = delta,
        nCutpoints = nCutpoints,
        minsplit = minsplit,
        size.limit = size.limit,
        maxdepth = maxdepth,
        details = details,
        alpha = alpha,
        model = model
      )

    partynode(id = strtoi(gsub("2", "0", name), base = 2),
              split = split,
              kids = kids)
  }
