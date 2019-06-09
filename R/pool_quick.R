pool_quick <- function(cm, sm){
  m <- nrow(cm)
  k <- ncol(cm)
  cp <- rep(1/m, m) %*% cm
  smn <- rep(1/m, m) %*% (sm^2)
  e <- cm - matrix(cp, nrow = m, ncol = k, byrow = TRUE)
  b <- (t(e) %*% e)/(m - 1)
  t <- smn + (1+1/m) * diag(b)
  sp <- sqrt(t)
  return(list(cp, sp))
}

# modified miceadds::summary.pool_mi that removed printing
summary.pool_mi <- function(object, alpha = 0.05, ...) {
  out <- data.frame(results = object$qbar, se = sqrt(diag(object$t)))
  crit <- stats::qt(alpha/2, object$df, lower.tail = FALSE)
  out$t <- object$tval
  out$p <- object$pval
  out$"(lower" <- out$results - crit * out$se
  out$"upper)" <- out$results + crit * out$se
  out$missInfo <- paste0(round(100 * object$fmi, 1), " %")
  out
}

quick_pool <- function (qhat, se) {
### REWRITE OF mice::pool AND mice::mice_df
  eps <- 1e-100
  m <- length(qhat)
  k <- length(qhat[[1]])
  im <- (1 + 1/m)
  um <- colMeans(do.call("rbind", se)^2)
  
  qhat2 <- do.call("rbind", qhat)
  qbar <- colMeans(qhat2)

  e <- sweep(qhat2, 2, qbar, `-`) 
  bm <- colSums(e^2)/(m - 1 + eps)
  t <- um + im * bm
  se2 <- sqrt(t)
  tval <- qbar/se2
  
  lambda <- im * (bm/t)
  
  eps2 <- 1e-04
  dfcom <- 1e+07
  lambda[lambda < eps2] <- eps2
  dfold <- (m - 1 + eps2)/lambda^2
  dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
  df <- dfold * dfobs/(dfold + dfobs)

  pval <- 2 * stats::pt(-abs(tval), df = df)
  
  list(results = qbar, 
       se = se2, 
       t = tval, 
       p = pval)
}