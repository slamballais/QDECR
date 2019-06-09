formula.vw_fastlm <- function(vw) vw$model$formula

nobs.vw_fastlm <- function(vw) vw$results$residuals$nrow

print.vw <- function(vw) qdecr_print_describe(vw$describe, verbose = TRUE)

qdecr_fwhm <- function(vw) return(vw$post$fwhm_est)

stacks <- function(vw) vw$stack$names




