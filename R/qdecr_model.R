qdecr_model <- function(model, prep_fun, md, id, vertex, margs, dir_tmp2) {
  prepvw <- list(model = model, data = md, id = id, 
                 vertex = vertex, margs = margs, path = dir_tmp2)
  class(prepvw) <- "prepvw"
  do.call2(prep_fun, list(prepvw = prepvw))
}

qdecr_backing_path <- function(prepvw, to_keep, to_remove) {
  prepvw$so <- c(to_keep, to_remove)
  prepvw$backing <- paste0(prepvw$path, "_", prepvw$so[1:4], "_backend")
  prepvw$backing_to_remove <- paste0(prepvw$path, "_", prepvw$so[5], "_backend")
  return(prepvw)
}

prep_fastlm <- function(prepvw) {
  to_keep <- c("coef", "se", "t", "p")
  to_remove <- "resid"
  prepvw <- qdecr_backing_path(prepvw, to_keep, to_remove)
  mf <- prepvw$margs
  if (is.null(mf$formula)) 
    stop("No `formula` set in margs.")
  iii <- c("formula", "data", "method")
  mf[iii] <- mf[!sapply(mf, is.symbol)][iii]
  mfz <- mf[match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)]
  if (!mf$method %in% 0:5) stop("The specified `method` for fastLm is not defined as a number between 0 and 5.")
  nr <- nrow(prepvw$data[[1]])
  mx <- lapply(prepvw$data, function(x) {
    mfz$data <- x
    mfz$data[, prepvw$vertex] <- 999
    do.call2("stats::model.frame", mfz)
  })
  nn <- length(prepvw$data)
  mt <- attr(mx[[nn]], "terms")
  w <- as.vector(model.weights(mx[[nn]]))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  mfz2 <- mfz
  mfz2$data <- prepvw$data[[nn]]
  mfz2$data[, prepvw$vertex] <- 888
  mx_test <- do.call2("stats::model.frame", mfz2)
  mx_test[, prepvw$vertex] <- 999
  if (!identical(mx[[nn]], mx_test)) stop ("Somewhere in your formula you specified a special term related to your vertex measure", 
                                           " (interaction, polynomial, AsIs, etc); `qdecr_fastlm` currently does not support this.")
  y <- stats::model.response(mx[[nn]], "numeric")
  if (nrow(mx[[nn]]) != nr) stop("The data that you are putting into the regression has missings! \n",
                                 "QDECR can't handle that yet; we will fix this soon!")
  ys <- if(identical(unname(y), rep(999, nrow(mx[[nn]])))) "LHS" else "RHS"
  if (ys == "LHS") { 
    mx_test2 <- mx_test
    mx_test2b <- stats::model.matrix(mx_test2, object = mt)
    if (Matrix::rankMatrix(mx_test2b) < ncol(mx_test2b)) stop ("The design matrix is NOT full rank. Please check if you have collinear columns in your data.")
  }
  if (stats::is.empty.model(mt)) stop("The provided model (to fastLm) is empty. Check your data + formula.")
  mm <- NULL
  prepvw$ff <- "vw_fastlm_slow"
  if (prepvw$vertex %in% colnames(attr(mt, "factors")) || ys == "LHS"){
    mm <- lapply(mx, stats::model.matrix, object = mt)
    ff <- "vw_fastlm"
    vw <- list(mm = mm, mf = mx[[1]], ff = ff, formula = mf$formula, vertex = prepvw$vertex, y = y, ys = ys, w = w, method = mf$method, backing = prepvw$backing, backing_to_remove = prepvw$backing_to_remove, so = prepvw$so)
  } else {
    warning("Your formula for `fastLm` contains complicated terms. \n",
            "We will rely on the slower implementation of our fastLm.")
    vw <- prepvw
  }
  class(vw) <- c(vw$ff, "vw")
  vw
  
}
