qdecr_decon <- function(x, y = environment()) {
  y <- as.list(substitute(x, env = y))
  y[-1] <- lapply(y[-1], eval)
  y
}

qdecr_extract <- function(margs, model){
  if (model == "QDECR::megha"){
    stop("MEGHA is not implemented yet.") 
  } else if (model %in% c("RcppEigen::fastLm", "stats::glm", 
                          "stats::lm", "survival::coxph", "default")) {
    if (!is.null(margs) && !is.symbol(margs$data)) margs$data else "fail"
    }
}

qdecr_setnames <- function(margs, model){
  m <- names(margs)
  if(is.null(m)) m <- rep("", length(margs))
  f <- methods::formalArgs(get2(model))
  f2 <- formals(get2(model))
  b <- which(m[-1] == "")
  m[b+1] <- f[!f %in% m[-1]][length(b)]
  names(margs) <- m
  margs <- c(margs, f2[!names(f2) %in% m])
  if(is.symbol(margs$`...`)) margs$`...` <- NULL
  margs
}

do.call2 <- function(what, args, ...) { # Copied from: https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  what <- get2(what)
  do.call(what, as.list(args), ...)
}

get2 <- function(x, ...) {
  if(is.character(x)){
    fn <- strsplit(x, "::")[[1]]
    x <- if(length(fn) == 1) {
      get(fn[[1]], envir = parent.frame(), mode = "function", ...)
    } else {
      get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function", ...)
    }
  }
  x
}