
qdecr_model <- function(model, md, id, vertex, margs, dir_tmp2){
  
  prepvw <- list(model = model, data = md, id = id, 
             vertex = vertex, margs = margs)
  class(prepvw) <- "prepvw"
  
  if (prepvw$model == "QDECR::megha") {
    model_qdecr_megha(prepvw)
  } else if (prepvw$model == "RcppEigen::fastLm") {
      prepvw$so <- c("coef", "se", "t", "p", "resid")
      prepvw$backing <- paste0(dir_tmp2, "_", prepvw$so[1:4], "_backend")
      prepvw$backing_to_remove <- paste0(dir_tmp2, "_", prepvw$so[5], "_backend")
      model_RcppEigen_fastLm(prepvw) 
  } else if (prepvw$model == "stats::glm") {
      model_stats_glm(prepvw)
  } else if (prepvw$model == "stats::lm") {
      model_stats_lm(prepvw)
  } else if (prepvw$model == "survival::coxph") {
      model_survival_coxph(prepvw)
  } else if (prepvw$model == "default") {
      prepvw$so <- c("coef", "se", "t", "p", "resid")
      prepvw$backing <- paste0(dir_tmp2, "_", prepvw$so[1:4], "_backend")
      prepvw$backing_to_remove <- paste0(dir_tmp2, "_", prepvw$so[5], "_backend")
      model_default(prepvw)
  }
}

model_default <- function(prepvw) {
  prepvw$ff <- "vw_default"
  prepvw$formula <- prepvw$margs$formula
  vw <- prepvw
  class(vw) <- c(vw$ff, "vw")
  vw
}

model_qdecr_megha <- function(prepvw){
  stop("`QDECR::megha` has not been implemented yet.")
}

model_RcppEigen_fastLm <- function(prepvw){
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
  
  mfz2 <- mfz
  mfz2$data <- prepvw$data[[nn]]
  mfz2$data[, prepvw$vertex] <- 888
  mx_test <- do.call2("stats::model.frame", mfz2)
  mx_test[, prepvw$vertex] <- 999
  
  if (!identical(mx[[nn]], mx_test)) stop ("Somewhere in your formula you specified a special term related to your vertex measure", 
    " (interaction, polynomial, AsIs, etc); `qdecr_fastlm` currently does not support this.")

  y <- model.response(mx[[nn]], "numeric")
  
  if (nrow(mx[[nn]]) != nr) stop("The data that you are putting into the regression has missings! \n",
                                 "QDECR can't handle that yet; we will fix this soon!")
  
  ys <- if(identical(unname(y), rep(999, nrow(mx[[nn]])))) "LHS" else "RHS"
   
  if (ys == "LHS") { 
    mx_test2 <- mx_test
    mx_test2b <- model.matrix(mx_test2, object = mt, contrasts)
    if (Matrix::rankMatrix(mx_test2b) < ncol(mx_test2b)) stop ("The design matrix is NOT full rank. Please check if you have collinear columns in your data.")
  }
  
  if (is.empty.model(mt)) stop("The provided model (to fastLm) is empty. Check your data + formula.")
  mm <- NULL
  prepvw$ff <- "vw_fastlm_slow"
  if (prepvw$vertex %in% colnames(attr(mt, "factors")) || ys == "LHS"){
    mm <- lapply(mx, model.matrix, object = mt, contrasts)
    ff <- "vw_fastlm"
    vw <- list(mm = mm, mf = mx[[1]], ff = ff, formula = mf$formula, vertex = prepvw$vertex, y = y, ys = ys, method = mf$method, backing = prepvw$backing, backing_to_remove = prepvw$backing_to_remove, so = prepvw$so)
  } else {
    warning("Your formula for `fastLm` contains complicated terms. \n",
            "We will rely on the slower implementation of our fastLm.")
    vw <- prepvw
  }
  class(vw) <- c(vw$ff, "vw")
  vw
}

model_stats_glm <- function(prepvw){
  stop("`stats::glm` has not been implemented yet.")
}

model_stats_lm <- function(prepvw){
  mf <- prepvw$margs
  ret.x <- mf$x
  ret.y <- mf$y
  if (is.null(mf$drop.unused.levels)) 
    mf$drop.unused.levels <- TRUE
  if (is.null(mf$formula)) 
    stop("No `formula` set in margs.")
  iii <- c("formula", "data", "subset", "weights", "na.action", "offset")
  mf[iii] <- mf[!sapply(mf, is.symbol)][iii]
  mfz <- mf[match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)]
  nr <- nrow(prepvw$data[[1]])
  mx <- lapply(prepvw$data, function(x) {
    mfz$data <- x
    mfz$data[, prepvw$vertex] <- 999
    do.call2("stats::model.frame", mfz)
  })
  nn <- length(prepvw$data)
  if (mf$method == "model.frame")
    stop("`margs$method` is equal to `model.frame` (for the lm call). Try `qr`.") else 
      if (mf$method != "qr")
    stop(gettextf("margs$method = '%s' is not supported in lm. Try 'qr'",
                     method), domain = NA)
  mt <- attr(mx[[nn]], "terms")
  y <- model.response(mx[[nn]], "numeric")
  ys <- if(identical(unname(y), rep(999, nr))) "LHS" else "RHS"
  w <- as.vector(model.weights(mx[[nn]]))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'margs$weights' must be a numeric vector (for the lm call)")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations) (in lm call)", 
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) stop("The provided model (to lm) is empty. Check your data + formula.")
  test <- model.matrix(mt, mx[[1]], contrasts) 
  test2 <- model.matrix(mt, mx[[nn]], contrasts)
  mm <- NULL
  prepvw$ff <- "vw_lm_fit_slow"
  if (prepvw$vertex %in% colnames(attr(mt, "factors")) || ys == "LHS"){
    # mfz2 <- mfz
    # mfz2$data <- prepvw$data[[1]]
    # mfz2$data[, prepvw$vertex] <- mfz$data[, prepvw$vertex]
    # mm2 <- do.call2("stats::model.frame", mfz2)
    if(identical(test, test2)){
      mm <- lapply(mx, model.matrix, object = mt, contrasts)
      ff <- if (is.null(w)) "vw_lm_fit" else "vw.lm_wfit"
      vw <- list(mm = mm, ff = ff, vertex = prepvw$vertex, ys = ys,
                 y = y, w = w, offset = offset,
                 singular.ok = prepvw$margs$singular.ok)
    } else {
      warning("Your formula for `lm` contains computed terms. \n", 
              "We will rely on the slower implementation of our lm.")
      vw <- prepvw
    }
  } else {
    warning("Your formula for `lm` contains complicated terms. \n",
            "We will rely on the slower implementation of our lm.")
    vw <- prepvw
  }
  class(vw) <- c(vw$ff, "vw")
  vw
}

model_survival_coxph <- function(prepvw){
  stop("`survival::coxph` has not been implemented yet.")
}

