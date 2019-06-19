vertexwise <- function(vw, analysis_fun, chunk = 1000) do.call2(analysis_fun, list(vw = vw, chunk = chunk))

analysis_chunkedlm <- function(vw, chunk) {
  
  # determine the number of vertices/loops and the number of variables
  n <- nrow(vw$mgh)
  nn <- ncol(vw$model$mm[[1]])
  m <- length(vw$model$mm)
  df <- nrow(vw$model$mm[[1]]) - nn
  
  # all your shared ram are belong to us
  c_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = vw$model$backing[1])
  s_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = vw$model$backing[2])
  t_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = vw$model$backing[3])
  p_vw <- bigstatsr::FBM(nn, n, init = 1, backingfile = vw$model$backing[4])
  r_vw <- bigstatsr::FBM(nrow(vw$model$mm[[1]]), n, init = 0, backingfile = vw$model$backing_to_remove[1])
  
  # mask all stuff
  mgh_is_0 <- fbm_row_is_0(vw$mgh, ncores = 1)
  masked <- !as.logical(vw$mask) | mgh_is_0
  
  # make chunk sequence
  iv <- which(!masked)
  cstart <- seq(1, length(iv), chunk)
  cend <- cstart + chunk - 1
  cend[length(cend)] <- length(iv)
  
  # run most of the regression
  if (vw$model$ys != "LHS") stop("The vertex measure has to be on the left hand side for analysis_fastlm.")
  XTX <- lapply(vw$model$mm, function(z) chol2inv(chol(crossprod(z))))
  XTXX <- lapply(1:m, function(z) tcrossprod(XTX[[z]], vw$model$mm[[z]]))
  
  # reduce load per core
  X <- vw$model$mm
  Ya <- vw$mgh
  
  # parallel loop
  cl <- parallel::makeForkCluster(vw$input$n_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  capture.output(pb <- txtProgressBar(0, length(cstart), style = 3), file = "/dev/null")
  
  foreach::foreach (i = seq_along(cstart), .combine = "c") %dopar% {
    setTxtProgressBar(pb, i)
    
    id <- iv[cstart[i]:cend[i]]
    
    Y <- t(Ya[id, ])
    
    # get coefs
    bhat <- lapply(XTXX, function(z) z %*% Y)
    
    # get residuals
    res <- lapply(1:m, function(z) Y - X[[z]] %*% bhat[[z]])
    
    # get se
    s2 <- lapply(res, function(z) colSums(z^2 / df))
    se <- lapply(1:m, function(z) do.call("cbind", lapply(s2[[z]], function(q) sqrt(diag(q * XTX[[z]])))))
    
    # pool and get t
    out <- quick_pool2(bhat, se = se)
    
    # write out
    c_vw[, id] <- out$results
    s_vw[, id] <- out$se
    t_vw[, id] <- out$t
    p_vw[, id] <- -1 * log10(out$p)
    r_vw[, id] <- Reduce("+", res) / length(res)
    
    NULL
  }
  
  parallel::stopCluster(cl)
  on.exit(invisible(NULL))
  
  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")
  return(out)
}