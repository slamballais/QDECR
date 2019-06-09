#' @importFrom bigstatsr rows_along
#' @importFrom bigstatsr cols_along

#' @export
fbm_row_mean <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) rowMeans(X[row.ind[ind], col.ind]), a.combine = "c", ncores = ncores, ind = seq_along(row.ind))
}

#' @export
fbm_col_mean <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) colMeans(X[row.ind, col.ind[ind]]), a.combine = "c", ncores = ncores, ind = seq_along(col.ind))
}

#' @export
fbm_row_sum <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) rowSums(X[row.ind[ind], col.ind]), a.combine = "c", ncores = ncores, ind = seq_along(row.ind))
}

#' @export
fbm_col_sum <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) colSums(X[row.ind, col.ind[ind]]), a.combine = "c", ncores = ncores, ind = seq_along(col.ind))
}

#' @export
fbm_row_sd <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind[ind], col.ind], 1, sd), a.combine = "c", ncores = ncores, ind = seq_along(row.ind))
}

#' @export
fbm_col_sd <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind, col.ind[ind]], 2, sd), a.combine = "c", ncores = ncores, ind = seq_along(col.ind))
}