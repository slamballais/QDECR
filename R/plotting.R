#'Opens Freeview
#'
#'Opens Freeview to view the surface of a specific stack/contrast
#'
#'This function directly calls Freeview so that you can view the surface of a specific
#'variable. As it opens Freeview directly, you will have access to all the functions
#'of Freeview by default. To continue using R, simply close Freeview.
#'
#'@param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#'@param stack Either a numeric or a string for your variable (use `stacks(vw)`)
#'@return NULL
#'@export

freeview <- function(vw, stack = NULL) {
  if (is.null(stack)) stop("`stack` not defined. Please choose: ", paste(stacks(vw), collapse = ", "))
  if (length(stack) > 1) stop("More than 1 stack specified.")
  if(is.character(stack)) stack <- which(stacks(vw) == stack)
  if(is.null(stack) || stack > length(stacks(vw))) stop("specified `stack` is not present in this dataset. Please choose: ", paste(stacks(vw), collapse = ", "))
  
  temp_mgh <- qdecr_read_coef(vw, stack)
  temp_mgh$x <- temp_mgh$x %MASK% qdecr_read_ocn_mask(vw, stack)
  if (all(temp_mgh$x == 0)) stop("Stack does not contain information (e.g. because of no significant findings), aborting plot.")
  temp_mgh_file <- paste0(vw$paths$dir_tmp2, ".stack", stack, ".temp_mgh.mgh")
  save.mgh(temp_mgh, temp_mgh_file)
  on.exit(file.remove(temp_mgh_file), add = TRUE)
  
  cmdStr <- paste0("freeview --surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated:overlay=", temp_mgh_file)
  system(cmdStr)
  
  invisible(NULL)
}

#'Histograms of the vertex measures
#'
#'Plots a histogram of the mean vertex measure, either by vertex or by subject
#'
#'This function works on top of the base `hist` function. It simply makes a histogram
#'of the mean vertex values either per vertex or per subject (if you specify `qtype = "subject"`)
#'
#'@param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#'@param qtype A string; either "vertex" for vertex-wise means, or "subject" for subject-wise means
#'@param xlab See `?hist`
#'@param main see `?hist`
#'@param ... Further arguments for `hist`
#'@return see `?hist`
#'@export
#'
hist.vw <- function(vw, qtype = c("vertex", "subject"), xlab = NULL, main = NULL, ...) {
  qtype <- match.arg(qtype)
  x <- if(qtype == "vertex") vw$post$mgh_description$vertex_mean else vw$post$mgh_description$subject_mean
  if (is.null(xlab)) xlab <- vw$input$measure
  if (is.null(main)) main <- qtype
  hist(x[], xlab = xlab, main = main, ...)
}

#'Opens Freeview to take snapshots
#'
#'Opens Freeview to take snapshots and compiles them into a plot
#'
#'This function directly calls Freeview and generates files (by default .tiff)
#'and additionally composes the images and plots them. Magick++ is required
#'for the composing step (see www.qdecr.com).
#'
#'@param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#'@param stack Either a numeric or a string for your variable (use `stacks(vw)`)
#'@param ext Extension of the image files that will be stored on disk
#'@param zoom Float that determines how far the brains are zoomed in
#'@param compose Logical; if TRUE, a single compiled image will be made (requires Magick++)
#'@param plot_brain Logical; if TRUE, returns a graphical device with the composed images
#'@param save_plot Logical; if TRUE, saves the composed image
#'@return NULL
#'@export

qdecr_snap <- function(vw, stack = NULL, ext = ".tiff", zoom = 1, compose = TRUE, plot_brain = TRUE, save_plot = TRUE) {
  if (is.null(stack)) stop("`stack` not defined. Please choose: ", paste(stacks(vw), collapse = ", "))
  if (length(stack) > 1) stop("More than 1 stack specified.")
  if(is.character(stack)) stack <- which(stacks(vw) == stack)
  if(is.null(stack) || stack > length(stacks(vw))) stop("specified `stack` is not present in this dataset. Please choose: ", paste(stacks(vw), collapse = ", "))
  if (substring(ext, 1, 1) != ".") ext <- paste0(".", ext)
  if (is.null(zoom) || !is.numeric(zoom)) stop("zoom is not a number. Please specify a zoom ratio (ideally below 1.4).")
  
  hemi <- vw$input$hemi
  name <- paste0(vw$paths$final_path, ".stack", stack)
  snap_order <- c("lateral", "medial", "superior", "inferior")
  if (hemi == "rh") snap_order[1:2] <- snap_order[2:1]
  snap_names <- paste0(name, ".", snap_order, ext)
  
  snap_cmd <- c(qsnap_zoom(zoom),
                qsnap(snap_names[1]),
                qsnap_a(180),
                qsnap(snap_names[2]),
                if (hemi == "lh") qsnap_a(180),
                qsnap_e(90),
                qsnap(snap_names[3]),
                qsnap_e(180),
                qsnap(snap_names[4]),
                "--quit")
  tfile <-  file.path(vw$paths$dir_tmp, "tmp_snapshot_qdecr.txt")
  write.table(snap_cmd, tfile, quote = F, row.names = F, col.names = F)
  on.exit(file.remove(tfile))
  
  temp_mgh <- qdecr_read_coef(vw, stack)
  temp_mgh$x <- temp_mgh$x %MASK% qdecr_read_ocn_mask(vw, stack)
  if (all(temp_mgh$x == 0)) stop("Stack does not contain information (e.g. because of no significant findings), aborting plot.")
  temp_mgh_file <- paste0(vw$paths$dir_tmp2, ".stack", stack, ".temp_mgh.mgh")
  save.mgh(temp_mgh, temp_mgh_file)
  on.exit(file.remove(temp_mgh_file), add = TRUE)
  
  cmdStr <- paste0("freeview --surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated:overlay=", temp_mgh_file, " -cmd ", tfile)
  system(cmdStr)
  
  if (compose) {
    if (!requireNamespace("magick", quietly = TRUE)) {
      stop("Package \"magick\" needed for the compiled image. Please install it. See the introductions on www.qdecr.com",
           call. = FALSE)
    }
    
    for (i in seq_along(snap_names)) {
      tt <- magick::image_read(snap_names[i])
      w <- magick::image_info(tt)$width
      h <- magick::image_info(tt)$height
      r <- magick::image_raster(tt)
      r$col <- r$col == "#000000ff"
      black <- with(r, tapply(col, x, all))
      ll <- rle(as.numeric(black))
      if (length(ll$lengths) != 3) stop("the zoom is too large to handle, aborting plot.")
      hw <- w/2
      cr <- hw - ll$lengths[2]
      crl <- ll$lengths[1] - cr/2
      crop_string <- paste0(hw, "x", h, "+", crl, "+", 0)
      if (i %in% 3:4) {
        black2 <- with(r, tapply(col, y, all))
        ll2 <- rle(as.numeric(black2))
        if (length(ll2$lengths) != 3) stop("the zoom is too large to handle, aborting plot.")
        hh <- h/2
        cr2 <- hh - ll2$lengths[2]
        ccl <- ceiling(cr2/2)
        crop_string <- paste0(hw, "x", hh, "+", crl, "+", ll2$lengths[1] - ccl)
      }
      nn <- magick::image_crop(tt, crop_string)
      nn2 <- if (i == 1) nn else c(nn2, nn)
    }
    ia1 <- magick::image_append(nn2[1:2])
    ia2 <- magick::image_append(nn2[3:4])
    p <- magick::image_append(c(ia1, ia2), stack = TRUE)
    if (plot_brain) print(plot(p))
    if (save_plot) magick::image_write(p, paste0(name, ".", "plot", ext))
    return(p)
  } else {
    return(invisible(NULL))
  }
}

qsnap <- function(x) paste0("--ss ", x)
qsnap_a <- function(n) paste("--camera Azimuth", n)
qsnap_e <- function(n) paste("--camera Elevation", n)
qsnap_zoom <- function(n) paste("--zoom", n)