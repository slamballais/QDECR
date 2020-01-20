#' @name fbm_functions_from_qdecr
#' @rdname fbm_functions_from_qdecr
#' 
#' @title fbm functions from qdecr
#' 
#' @param X the FBM (file-backed matrix) object
#' @param ncores the number of cores to be used
#' @param row.ind the rows to go over (by default, it goes over all rows)
#' @param col.ind the columns to go over (by default, it goes over all columns)
#' @param row.mask a separate mask for the rows (added because this was useful internally)
#' @param col.mask a separate mask for the columns (added because this was useful internally)
#' @importFrom bigstatsr rows_along
#' @importFrom bigstatsr cols_along
#' @return an mgh object of the corresponding map
NULL

#' @rdname fbm_functions_from_qdecr
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

#' @rdname fbm_functions_from_qdecr
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

#' @rdname fbm_functions_from_qdecr
#' @export
fbm_row_is_0 <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind[ind], col.ind], 1, function(q) any(q == 0)), a.combine = "c", ncores = ncores, ind = seq_along(row.ind))
}

#' @rdname fbm_functions_from_qdecr
#' @export
fbm_col_is_0 <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind, col.ind[ind]], 2, function(q) any(q == 0)), a.combine = "c", ncores = ncores, ind = seq_along(col.ind))
}

#' @rdname fbm_functions_from_qdecr
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

#' @rdname fbm_functions_from_qdecr
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

#' @rdname fbm_functions_from_qdecr
#' @export
fbm_row_sd <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind[ind], col.ind], 1, stats::sd), a.combine = "c", ncores = ncores, ind = seq_along(row.ind))
}

#' @rdname fbm_functions_from_qdecr
#' @export
fbm_col_sd <- function(X, ncores = 1, 
                         row.ind = bigstatsr::rows_along(X), col.ind = bigstatsr::cols_along(X), 
                         row.mask = NULL, col.mask = NULL) {
  if (is.numeric(row.mask)) row.mask <- as.logical(row.mask)
  if (is.numeric(col.mask)) col.mask <- as.logical(col.mask)
  if(!is.null(row.mask)) row.ind <- row.ind[row.mask] 
  if(!is.null(col.mask)) col.ind <- col.ind[col.mask] 
  bigstatsr::big_apply(X, a.FUN = function(X, ind) apply(X[row.ind, col.ind[ind]], 2, stats::sd), a.combine = "c", ncores = ncores, ind = seq_along(col.ind))
}