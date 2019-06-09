#' Saves "vw" objects from the QDECR package
#'
#' @param vw Object of class "vw". Potentially other classes related to specific methods, e.g. "vw_fastlm"
#' @return NULL
#' @export

qdecr_save <- function(vw, ...) UseMethod("qdecr_save", vw)

qdecr_save.vw_fastlm <- function(vw, save_data = TRUE, ...) {
  if (!save_data) {
    vw$input$md <- NULL
    vw$input$margs$data <- NULL
    vw$model$mm <- NULL
  }
  NextMethod()
}

qdecr_save.vw <- function(vw, file = vw$input$project2, project_dir = TRUE, export_vertex = FALSE, save_data = TRUE) {
  if (project_dir) file <- file.path(vw$paths$dir_out, file)
  saveRDS(vw, paste0(file, ".rds"))
  invisible(NULL)
}
