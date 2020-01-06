##' @name qdecr_read
##' @rdname qdecr_read
##' 
##' @title qdecr_read function family
##' 
##' @param vw Output object from a qdecr analysis
##' @param stack stack number of name of the stack of interest
##' 
##' @return an mgh object of the corresponding map
NULL

##' @rdname qdecr_read
##' @export
qdecr_read_ocn <- function(vw, stack) load.mgh(vw$stack$ocn.mgh[[stack]])

##' @rdname qdecr_read
##' @export
qdecr_read_ocn_mask <- function(vw, stack) qdecr_read_ocn(vw, stack)$x > 0

##' @rdname qdecr_read
##' @export
qdecr_read_coef <- function(vw, stack) load.mgh(vw$stack$coef[[stack]])

##' @rdname qdecr_read
##' @export
qdecr_read_p <- function(vw, stack) load.mgh(vw$stack$p[[stack]])

##' @rdname qdecr_read
##' @export
qdecr_read_t <- function(vw, stack) load.mgh(vw$stack$t[[stack]])

##' @rdname qdecr_read
##' @export
qdecr_read_se <- function(vw, stack) load.mgh(vw$stack$se[[stack]])
