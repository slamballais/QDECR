qdecr_start_describe <- function() {
  emp <- as.matrix(data.frame(type = character(),
                              name = character(),
                              value = character(),
                              stringsAsFactors = FALSE))
  description <- list(call = emp,
                      input = emp,
                      paths = emp,
                      data = emp,
                      model = emp,
                      mask = emp,
                      stack = emp,
                      post = emp)                      
  description
}

qdecr_pre_describe <- function(vw, mask, verbose = TRUE) {
  info <- qdecr_start_describe()
  info$input <- rbind(info$input, 
                      c("input", "Hemisphere", vw$input$hemi),
                      c("input", "Project name", vw$input$project),
                      c("input", "Project final name", vw$input$project2), 
                      c("input", "Number of cores", vw$input$n_cores),
                      c("input", "Target", vw$input$target), 
                      c("input", "dir_out_tree", vw$input$dir_out_tree),
                      c("input", "clobber", vw$input$clobber),
                      c("input", "fwhm", vw$input$fwhm))
  info$paths <- rbind(info$paths,
                      c("paths", "Subjects dir", vw$paths$dir_subj),
                      c("paths", "Freesurfer home dir", vw$paths$dir_fshome),
                      c("paths", "Temp dir", vw$paths$dir_tmp),
                      c("paths", "Output dir", vw$paths$dir_out),
                      c("paths", "Path to default mask", if(file.exists(vw$paths$mask_path)) vw$paths$mask_path else "Not found!"))
  info$data <- rbind(info$data, 
                     c("data", "N subjects", length(unique(vw$input$md[[1]][[vw$input$id]]))),
                     c("data", "N data points", nrow(vw$input$md[[1]])), 
                     c("data", "N datasets", length(vw$input$md)),
                     c("data", "Vertices loaded", nrow(vw$mgh)))
  im <- if (vw$input$model != "default") vw$input$model else vw$model$margs[[1]]
  info$model <- rbind(info$model,
                      c("model", "Model", im), 
                      c("model", "Vertex data", vw$input$measure),
                      c("model", "Formula", paste(deparse(vw$model$formula), collapse = "")))                      
  info$mask <- rbind(info$mask,
                     c("mask", "Mask origin", if (!is.null(mask)) "User defined." else vw$paths$mask_path),  
                     c("mask", "Masked vertices", sum(vw$mask)))       
  if (verbose) qdecr_print_describe(info, verbose = verbose)
  info
}
                
qdecr_print_describe <- function(info, verbose = TRUE) {
  if(!verbose) return(NULL)
  info[] <- lapply(info, as.data.frame, stringsAsFactors = FALSE)
  info <- do.call("rbind", info)
  xchar1 <- nchar(info$type)
  xchar2 <- nchar(info$name)
  print_length1 <- max(xchar1) - xchar1
  print_length2 <- max(xchar2) - xchar2 + 10
  print_out <- paste0("[", info$type,  strrep(" ", print_length1), "] ", info$name, ":", strrep(" ", print_length2), "[", info$value, "]")
  for (n in print_out) message(n)
  invisible(NULL)
}

qdecr_post_describe <- function(vw, verbose = TRUE) {
  info <- vw$describe
  ls <- length(vw$stack$names)
  info$stack <- rbind(info$stack, cbind(rep("stack", ls), paste("Stack", 1:ls), vw$stack$names))
  info$post <- rbind(info$post, 
                     c("post", "Final mask path", vw$paths$final_mask_path),
                     c("post", "Final N vertices", sum(vw$post$final_mask)),
                     c("post", "Estimated fwhm", vw$post$fwhm_est),
                     c("post", paste("Mean", vw$input$measure, "per vertex"), mean(vw$post$mgh_description$vertex_mean)),
                     c("post", paste("SD", vw$input$measure, "per vertex"), mean(vw$post$mgh_description$vertex_sd)),
                     c("post", paste("Mean", vw$input$measure, "per subject"), mean(vw$post$mgh_description$subject_mean)),
                     c("post", paste("SD", vw$input$measure, "per subject"), mean(vw$post$mgh_description$subject_sd)))
  if (verbose) qdecr_print_describe(info, verbose = verbose)
  info
}