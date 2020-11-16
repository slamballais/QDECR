#' Updates the internally stored paths in a vw object
#'
#' @param vw A `vw` object (QDECR output), e.g. from `qdecr_load`.
#' @param dir_fshome character; path to FREESURFER_HOME to be used.
#' @param dir_subj character; path to SUBJECTS_DIR to be used.
#' @param dir_project character; path to the project directory (defaults to the parent dir of the .rds file).
#' @param mask_path character; if specified and QDECR_mask == FALSE, this will become the mask path.
#' @param qdecr_mask logical (default = TRUE); should the QDECR mask be used as the mask?
#' @param overwrite logical (default = TRUE); should the project .rds file be overwritten?
#' @return The vw object, but with updated paths.
#' @export

qdecr_update_path <- function(vw, dir_fshome, dir_subj, dir_project = dirname(vw$paths$rds), mask_path = NULL, qdecr_mask = TRUE, overwrite = TRUE) {
  to_sub <- vw$paths[["final_path"]]
  path_names <- names(vw$paths)
  path <- dir_project
  vw$paths <- as.list(sub(to_sub, path, vw$paths))
  names(vw$paths) <- path_names
  vw$paths$dir_tmp <- path
  vw$paths$dir_tmp2 <- path
  vw$stack[-1] <- rapply(vw$stack[-1], function(x) sub(to_sub, path, x), how = "replace")
  vw$paths$dir_fshome <- dir_fshome
  vw$paths$dir_subj <- dir_subj
  if(qdecr_mask) mask_path <- system.file("extdata", paste0(vw$input$hemi, ".fsaverage.cortex.mask.mgh"), package = "QDECR")
  if(!is.null(mask_path)) vw$paths$mask_path <- mask_path
  if (overwrite) {
    message("Overwriting .rds file")
    saveRDS(vw, vw$paths$rds)
  }
  return(vw)
}


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


