#' \code{message} function with verbose option
#'
#' This is a wrapper for \code{message}. It calls \code{message} if \code{verbose}
#' is set to true.
#'
#' @param verbose Logical. If TRUE, \code{message} is called.
#' @param ... Any input for \code{message}

catprompt <- function(string, n = 60, base = "-", center = TRUE, addSpace = TRUE, start = 5, verbose = TRUE){
  if(addSpace) string <- paste0(" ", string, " ")
  ns <- nchar(string)
  if(ns > n) n <- ns
  phase <- n - ns
  if (center) {
    start <- end <- phase / 2
  } else {
    if (phase < start) n <- n + (start - phase)
    phase <- ns - n
    end <- n - start - ns
  }
  br1 <- paste0("\n", paste(rep(base, n), collapse = ""))
  sr <- paste(rep(base, start), collapse = "")
  er <- paste(rep(base, end), collapse = "")
  main <- paste0(sr, string, er)
  if (nchar(main) < n) main <- paste0(main, rep(base, n - nchar(main)))
  br2 <- paste(rep(base, n), collapse = "")
  cat(br1, "\n")
  cat(main, "\n")
  cat(br2, "\n\n")
}

message2 <- function(..., verbose = TRUE) {
  if (verbose) message(...)
  invisible()
}

collapse <- function(..., collapse = NULL) paste(..., collapse = collapse)

