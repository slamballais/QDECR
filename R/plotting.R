#'Opens Freeview
#'
#'Opens Freeview to view the surface of a specific stack/contrast
#'
#'This function directly calls Freeview so that you can view the surface of a specific
#'variable. As it opens Freeview directly, you will have access to all the functions
#'of Freeview by default. To continue using R, simply close Freeview.
#'
#' @param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#' @param stack Either a numeric or a string for your variable (use `stacks(vw)`)
#' @param type The type of map to plot (currently for OLS). Default is `coef`, others are `se`, `t` and `p`
#' @param sig Logical; if TRUE (default), only the significant clusters are shown
#' @param ... any arguments that freeview (on command line) normally takes for manipulating surface files
#' @return NULL
#' @export

freeview <- function(vw, stack = NULL, type = c("coef", "se", "t", "p"), sig = TRUE, ...) {
  if (is.null(stack)) stop("`stack` not defined. Please choose: ", paste(stacks(vw), collapse = ", "))
  if (length(stack) > 1) stop("More than 1 stack specified.")
  if(is.character(stack)) stack <- which(stacks(vw) == stack)
  if(is.null(stack) || stack > length(stacks(vw))) stop("specified `stack` is not present in this dataset. Please choose: ", paste(stacks(vw), collapse = ", "))
  type <- match.arg(type)
  empty_val <- if(type == "p") 1 else 0

  read_fun <- switch(type, 
                     coef = qdecr_read_coef, 
                     se = qdecr_read_se, 
                     t = qdecr_read_t,
                     p = qdecr_read_p)
  
  temp_mgh <- read_fun(vw, stack)
  if (all(temp_mgh$x == empty_val)) stop("Stack does not contain information (e.g. because of no significant findings), aborting plot.")
  if (sig) {
    p_mgh <- qdecr_read_ocn_mask(vw, stack)
    temp_mgh$x[!p_mgh] <- empty_val
    if (all(temp_mgh$x == empty_val)) stop("No information in the stack passed the threshold (i.e. `p_thr` was set too strict), aborting plot.")
  }

  temp_mgh_file <- paste0(vw$paths$dir_tmp2, ".stack", stack, ".temp_mgh.mgh")
  save.mgh(temp_mgh, temp_mgh_file)
  on.exit(file.remove(temp_mgh_file), add = TRUE)
  
  surface <- paste0("--surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated")

  freeview_args <- list(...)
  freeview_args$surface <- surface
  freeview_args$overlay <- temp_mgh_file
  
  if (is.null(freeview_args$overlay_method)) {
    freeview_args$overlay_method <- "linearopaque"
    temp_vals <- temp_mgh$x[temp_mgh$x != empty_val]
    freeview_args$overlay_threshold <- c(min(abs(temp_vals)), max(abs(temp_vals)))
  }
  
  surf_cmd <- do.call(freeview_surf_cmd, freeview_args)
  cmd_str <- paste("freeview", surf_cmd)
  cat("Opening Freeview on the command line. Pop-up incoming.\nCall:", cmd_str, "\n")
  system(cmd_str)
  
  invisible(NULL)
}

freeview_surf_cmd <- function(surface,
                              curvature = NULL, 
                              curvature_method = NULL,
                              overlay = NULL,
                              overlay_reg = NULL,
                              overlay_method = NULL,
                              overlay_color = NULL,
                              overlay_threshold = NULL,
                              overlay_frame = NULL,
                              correlation = NULL,
                              color = NULL,
                              edgecolor = NULL,
                              edgethickness = NULL,
                              annot = NULL,
                              annot_outline = NULL,
                              name = NULL,
                              offset = NULL,
                              visible = NULL,
                              vector = NULL,
                              target_surf = NULL,
                              label = NULL,
                              label_outline = NULL,
                              label_color = NULL,
                              label_centroid = NULL,
                              label_visible = NULL,
                              spline = NULL,
                              vertex = NULL,
                              vertexcolor = NULL,
                              goto = NULL,
                              hide_in_3d = NULL,
                              all = NULL) {
  string <- surface
  if(!is.null(curvature)) string <- paste0(string, ":curvature=", curvature)
  if(!is.null(curvature_method)) string <- paste0(string, ":curvature_method=", curvature_method)
  if(!is.null(overlay)) string <- paste0(string, ":overlay=", overlay)
  if(!is.null(overlay_reg)) string <- paste0(string, ":overlay_reg=", overlay_reg)
  if(!is.null(overlay_method)) string <- paste0(string, ":overlay_method=", overlay_method)
  if(!is.null(overlay_color)) string <- paste0(string, ":overlay_color=", collapse(overlay_color, collapse = ","))
  if(!is.null(overlay_threshold)) string <- paste0(string, ":overlay_threshold=", collapse(overlay_threshold, collapse = ","))
  if(!is.null(overlay_frame)) string <- paste0(string, ":overlay_frame=", overlay_frame)
  if(!is.null(correlation)) string <- paste0(string, ":correlation=", correlation)
  if(!is.null(color)) string <- paste0(string, ":color=", collapse(color, collapse = ","))
  if(!is.null(edgecolor)) string <- paste0(string, ":edgecolor=", collapse(edgecolor, collapse = ","))
  if(!is.null(edgethickness)) string <- paste0(string, ":edgethickness=", edgethickness)
  if(!is.null(annot)) string <- paste0(string, ":annot=", annot)
  if(!is.null(annot_outline)) string <- paste0(string, ":annot_outline=", annot_outline)
  if(!is.null(name)) string <- paste0(string, ":name=", name)
  if(!is.null(offset)) string <- paste0(string, ":offset=", collapse(offset, collapse = ","))
  if(!is.null(visible)) string <- paste0(string, ":visible=", visible)
  if(!is.null(vector)) string <- paste0(string, ":vector=", vector)
  if(!is.null(target_surf)) string <- paste0(string, ":target_surf=", target_surf)
  if(!is.null(label)) string <- paste0(string, ":label=", label)
  if(!is.null(label_outline)) string <- paste0(string, ":label_outline=", label_outline)
  if(!is.null(label_color)) string <- paste0(string, ":label_color=", label_color)
  if(!is.null(label_centroid)) string <- paste0(string, ":label_centroid=", label_centroid)
  if(!is.null(label_visible)) string <- paste0(string, ":label_visible=", label_visible)
  if(!is.null(spline)) string <- paste0(string, ":spline=", spline)
  if(!is.null(vertex)) string <- paste0(string, ":vertex=", vertex)
  if(!is.null(vertexcolor)) string <- paste0(string, ":vertexcolor=", collapse(vertexcolor, collapse = ","))
  if(!is.null(goto)) string <- paste0(string, ":goto=", goto)
  if(!is.null(hide_in_3d)) string <- paste0(string, ":hide_in_3d=", hide_in_3d)
  if(!is.null(all)) string <- paste0(string, ":all=", all)
  string
}

#'Histograms of the vertex measures
#'
#'Plots a histogram of the mean vertex measure, either by vertex or by subject
#'
#'This function works on top of the base `hist` function. It simply makes a histogram
#'of the mean vertex values either per vertex or per subject (if you specify `qtype = "subject"`)
#'
#' @param x The output object of a qdecr call (e.g. qdecr_fastlm)
#' @param qtype A string; either "vertex" for vertex-wise means, or "subject" for subject-wise means
#' @param xlab See `?hist`
#' @param main see `?hist`
#' @param ... Further arguments for `hist`
#' @return see `?hist`
#' @method hist vw
#' @export
#'
hist.vw <- function(x, qtype = c("vertex", "subject"), xlab = NULL, main = NULL, ...) {
  qtype <- match.arg(qtype)
  y <- if(qtype == "vertex") x$post$mgh_description$vertex_mean else x$post$mgh_description$subject_mean
  if (is.null(xlab)) xlab <- x$input$measure
  if (is.null(main)) main <- qtype
  graphics::hist(y[], xlab = xlab, main = main, ...)
}

#'Opens Freeview to take snapshots
#'
#'Opens Freeview to take snapshots and compiles them into a plot
#'
#'This function directly calls Freeview and generates files (by default .tiff)
#'and additionally composes the images and plots them. Magick++ is required
#'for the composing step (see www.qdecr.com).
#'
#' @param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#' @param stack Either a numeric or a string for your variable (use `stacks(vw)`)
#' @param type The type of map to plot (currently for OLS). Default is `coef`, others are `se`, `t` and `p`
#' @param ext Extension of the image files that will be stored on disk
#' @param zoom Float that determines how far the brains are zoomed in
#' @param compose Logical; if TRUE, a single compiled image will be made (requires Magick++)
#' @param plot_brain Logical; if TRUE, returns a graphical device with the composed images
#' @param save_plot Logical; if TRUE, saves the composed image
#' @param sig Logical; if TRUE, only the significant clusters are shown
#' @param ... any arguments that freeview (on command line) normally takes for manipulating surface files
#' @return NULL
#' @export

qdecr_snap <- function(vw, stack = NULL, type = c("coef", "se", "t", "p"), ext = ".tiff", zoom = 1, compose = TRUE, plot_brain = TRUE, save_plot = TRUE, sig = TRUE, ...) {
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
  
  snap_cmd <- c("--viewport 3d",
                qsnap_zoom(zoom),
                qsnap(snap_names[1]),
                qsnap_a(180),
                qsnap(snap_names[2]),
                if (hemi == "lh") qsnap_a(180),
                qsnap_e(90),
                qsnap(snap_names[3]),
                qsnap_e(180),
                qsnap(snap_names[4]),
                "--quit")
  tfile <- file.path(vw$paths$dir_tmp, "tmp_snapshot_qdecr.txt")
  utils::write.table(snap_cmd, tfile, quote = F, row.names = F, col.names = F)
  on.exit(file.remove(tfile))
  
  type <- match.arg(type)
  empty_val <- if(type == "p") 1 else 0
  
  read_fun <- switch(type, 
                     coef = qdecr_read_coef, 
                     se = qdecr_read_se, 
                     t = qdecr_read_t,
                     p = qdecr_read_p)
  
  temp_mgh <- read_fun(vw, stack)
  if (all(temp_mgh$x == empty_val)) stop("Stack does not contain information (e.g. because of no significant findings), aborting plot.")
  if (sig) {
    p_mgh <- qdecr_read_ocn_mask(vw, stack)
    temp_mgh$x[!p_mgh] <- empty_val
    if (all(temp_mgh$x == empty_val)) stop("No information in the stack passed the threshold (i.e. `p_thr` was set too strict), aborting plot.")
  }

  temp_mgh_file <- paste0(vw$paths$dir_tmp2, ".stack", stack, ".temp_mgh.mgh")
  save.mgh(temp_mgh, temp_mgh_file)
  on.exit(file.remove(temp_mgh_file), add = TRUE)
  
  surface <- paste0("--surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated")
  
  freeview_args <- list(...)
  freeview_args$surface <- surface
  freeview_args$overlay <- temp_mgh_file
  
  if (is.null(freeview_args$overlay_method)) {
    freeview_args$overlay_method <- "linearopaque"
    temp_vals <- temp_mgh$x[temp_mgh$x != empty_val]
    freeview_args$overlay_threshold <- c(min(abs(temp_vals)), max(abs(temp_vals)))
  }
  
  surf_cmd <- do.call(freeview_surf_cmd, freeview_args)
  cmd_str <- paste("freeview", surf_cmd, "-cmd", tfile)
  cat("Opening Freeview on the command line. Pop-up incoming.\nCall:", cmd_str, "\n")
  cat(tfile, "contains:", snap_cmd, "\n")
  system(cmd_str)
  
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
    if (plot_brain) print(graphics::plot(p))
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