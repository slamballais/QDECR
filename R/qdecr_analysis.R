vertexwise <- function(vw, split, ...) UseMethod("vertexwise", vw)

vertexwise.vw_default <- function(model, input, mask, mgh, n_cores, dir_tmp, project) {
    # determine the number of vertices/loops and the number of variables
  n <- nrow(mgh)
  nn <- 3
  
  # all your shared ram are belong to us
  c_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[1])
  s_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[2])
  t_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[3])
  p_vw <- bigstatsr::FBM(nn, n, init = 1, backingfile = model$backing[4])
  
  catprompt(model$backing_to_remove[1])
  
  r_vw <- bigstatsr::FBM(nrow(model$data[[1]]), n, init = 0, backingfile = model$backing_to_remove[1])
  
  # vw
  cl <- parallel::makeForkCluster(n_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  capture.output(pb <- txtProgressBar(0, n, style = 3), file = "/dev/null")
  
  foreach::foreach (i = seq_len(n), .combine = "c") %dopar% { 
    setTxtProgressBar(pb, i)
    if (mask[i] != 0 && !(0 %in% mgh[i,])) {
      ma <- model$margs
      dd <- ma$data    
      fun <- deparse(ma[[1]])
      ma[[1]] <- NULL
      ma$data <- NULL  
      
      # integrate vertex data
      x <- mgh[i, ]
      for (z in seq_along(dd)) dd[[z]][, input$vertex] <- x
      
      # reset formula
      ma$formula <- as.formula(substitute(ma$formula))
      
      # run model
      fit <- lapply(dd, function(z) do.call(fun, c(list(z), ma)))
      
      # extract coefficients, residuals and standard errors
      rl <- lapply(fit, function(z) z$residuals)
      lsum <- lapply(fit, function(z) summary(z)$coefficients)
      cl <- lapply(lsum, `[`,, 1)
      sl <- lapply(lsum, `[`,, 2)
      
      # apply pooling
      out <- QDECR::summary.pool_mi(miceadds::pool_mi(cl, se = sl))
      
      # fill tables
      c_vw[, i] <- out$results
      s_vw[, i] <- out$se
      t_vw[, i] <- out$t
      p_vw[, i] <- -1 * log10(out$p)
      r_vw[, i] <- Reduce("+", rl) / length(rl)
      
    }
    NULL
  }
  parallel::stopCluster(cl)
  
  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")
  return(out)
}

vertexwise.vw_fastlm <- function(model, input, mask, mgh, n_cores, dir_tmp, project){
  
  # determine the number of vertices/loops and the number of variables
  n <- nrow(mgh)
  nn <- ncol(model$mm[[1]])
  m <- length(model$mm)
  
  # all your shared ram are belong to us
  c_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[1])
  s_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[2])
  t_vw <- bigstatsr::FBM(nn, n, init = 0, backingfile = model$backing[3])
  p_vw <- bigstatsr::FBM(nn, n, init = 1, backingfile = model$backing[4])
  r_vw <- bigstatsr::FBM(nrow(model$mm[[1]]), n, init = 0, backingfile = model$backing_to_remove[1])
  
  # vw
  cl <- parallel::makeForkCluster(n_cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  capture.output(pb <- txtProgressBar(0, n, style = 3), file = "/dev/null")
  
  foreach::foreach (i = seq_len(n), .combine = "c") %dopar% { 
    setTxtProgressBar(pb, i)
    mghi <- mgh[i, ]
    
    if (mask[i] != 0 && !(0 %in% mghi)) {
      
      # integrate vertex data
      if (model$ys == "LHS") {
        model$y <- mghi
      } else {
        for (j in seq_along(model$mm)){
          model$mm[[j]][, input$vertex] <- mghi
        }
      }
      
      # run fastLmPure
      fit <- lapply(model$mm, RcppEigen::fastLmPure, model$y, method = model$method)
      
      # write out the residuals
      r_vw[, i] <- Reduce("+", lapply(fit, function(x) x$residuals)) / m
      
      # extract coefficients and standard errors
      cl <- lapply(fit, coef)
      sl <- lapply(fit, function(x) x$se)
      
      # apply pooling
      out <- quick_pool(cl, se = sl)
      
      # fill tables
      c_vw[, i] <- out$results
      s_vw[, i] <- out$se
      t_vw[, i] <- out$t
      p_vw[, i] <- -1 * log10(out$p)
       
    }
    NULL
  }
  parallel::stopCluster(cl)
  
  out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
  names(out) <- c("coefficients", "standard_errors", "t_values", "p_values", "residuals")
  return(out)
}
# 
# vertexwise.vw_lm_fit <- function(vw, n_cores, dir_tmp, project){
#   
#   # determine the number of vertices/loops and the number of variables
#   n <- ncol(vw$mgh)
#   nn <- ncol(vw$mm[[1]])
#   
#   # all your shared ram are belong to us
#   lock <- tempfile()
#   c_vw <- FBM(n, nn, init = 0, backingfile = paste0(dir_tmp, "/", project, "_coef_backend"))
#   s_vw <- FBM(n, nn, init = 0, backingfile = paste0(dir_tmp, "/", project, "_se_backend"))
#   t_vw <- FBM(n, nn, init = 0, backingfile = paste0(dir_tmp, "/", project, "_t_backend"))
#   p_vw <- FBM(n, nn, init = 0, backingfile = paste0(dir_tmp, "/", project, "_p_backend"))
#   r_vw <- FBM(n, nn, init = 0, backingfile = paste0(dir_tmp, "/", project, "_resid_backend"))
#   
#   # vw
#   cl <- parallel::makeForkCluster(n_cores, outfile = "")
#   doParallel::registerDoParallel(cl)
#   
#   pb <- txtProgressBar(0, n, style = 2)
#   
#   foreach (i = icount(n)) %dopar% { 
#     setTxtProgressBar(pb, i)
#     
#     # integrate vertex data
#     x <- vw$mgh[, i]
#     if (vw$ys == "LHS") {
#       vw$y <- x
#     } else {
#       for (j in seq_along(vw$mm)){
#         vw$mm[[j]][, vw$vertex] <- x
#       }
#     }
#     
#     # run lm.fit
#     fit <- lapply(vw$mm, lm.fit, vw$y, vw$offset, singular.ok = vw$singular.ok)
#     
#     # extract coefficients, residuals and standard errors
#     cl <- lapply(fit, coef)
#     rl <- lapply(fit, function(x) x$residuals)
#     sl <- lapply(seq_along(vw$mm), function(x) {
#       qq <- sum((vw$y - fitted(fit[[x]]))^2)/(nrow(vw$mm[[x]])-ncol(vw$mm[[x]]))
#       sqrt(diag(solve(crossprod(vw$mm[[x]]))))*qq
#     })
#     
#     # apply pooling
#     out <- QDECR::summary.pool_mi(miceadds::pool_mi(cl, se = sl))
#     
#     # fill tables
#     locked <- flock::lock(lock)
#     c_vw[i, ] <- out$results
#     s_vw[i, ] <- out$se
#     t_vw[i, ] <- out$t
#     p_vw[i, ] <- out$p
#     r_vw[i, ] <- Reduce("+", rl) / length(rl)
#     flock::unlock(locked)
#     NULL
#   }
#   parallel::stopCluster(cl)
#   
#   out <- list(c_vw, s_vw, t_vw, p_vw, r_vw)
# }
# 
# vertexwise.vw_lm_wfit <- function(vw){
#   if (vw$ys == "LHS") {
#     y <- vdata
#   } else {
#     for (i in seq_along(vw$mm)){
#       vw$mm[[i]][, vw$vertex] <- vdata
#     }
#   }
#   fit <- lapply(vw$mm, lm.wfit, y, vw$w, vw$offset, singular.ok = vw$singular.ok)
#   ll <- lapply(fit, summary)
#   cl <- lapply(ll, function(x) x$coef[,1])
#   sl <- lapply(ll, function(x) x$coef[,2])
#   out <- miceadds::pool_mi(cl, se = sl)
#   out
# }