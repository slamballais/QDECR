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


qdecr_read_ocn <- function(vw, stack) load.mgh(vw$stack$ocn.mgh[[stack]])
qdecr_read_ocn_mask <- function(vw, stack) qdecr_read_ocn(vw, stack)$x > 0

qdecr_read_coef <- function(vw, stack) load.mgh(vw$stack$coef[[stack]])
qdecr_read_p <- function(vw, stack) load.mgh(vw$stack$p[[stack]])
qdecr_read_t <- function(vw, stack) load.mgh(vw$stack$t[[stack]])
qdecr_read_se <- function(vw, stack) load.mgh(vw$stack$se[[stack]])

`%AND%` <- function(map1, map2) map1 & map2
`%OR%` <- function(map1, map2) map1 | map2
`%XOR%` <- function(map1, map2) (map1 & !map2) | (!map1 & map2)
`%MASK%` <- function(map1, map2) map1 * map2

AND <- function(map1, map2) map1 %AND% map2
OR <- function(map1, map2) map1 %OR% map2
XOR <- function(map1, map2) map1 %XOR% map2
MASK <- function(map1, map2) map1 %MASK% map2

hist.vw <- function(vw, qtype = c("vertex", "subject"), xlab = NULL, main = NULL, ...) {
  qtype <- match.arg(qtype)
  x <- if(qtype == "vertex") vw$post$mgh_description$vertex_mean else vw$post$mgh_description$subject_mean
  if (is.null(xlab)) xlab <- vw$input$measure
  if (is.null(main)) main <- qtype
  hist(x[], xlab = xlab, main = main, ...)
}

qdecr_clusters <- function(vw, name = "aparc.annot") {
  file <- paste0(vw$input$hemi, ".", name)
  path <- file.path(vw$paths$dir_subj, vw$input$target, "label", file)

  if (!file.exists(path)) stop("Provided annotation file does not exist.")

  annot <- load.annot(path)

  ocn_annot <- lapply(seq_along(vw$stack$names), function(x) {
    ocn <- qdecr_read_ocn(vw, x)
    if (any(ocn$x > 0)) annots <- lapply(seq_len(max(ocn$x)), function(i) annot$vd_label[ocn$x == i])
  })
  ocn_annot2 <- do.call(c, ocn_annot[!sapply(ocn_annot, is.null)])
  ocn_annot3 <- lapply(ocn_annot2, table)

  ta <- table(annot$vd_label)
  nta <- names(ta)

  ocn_annot4 <- lapply(ocn_annot3, function(x) {
    nx <- names(x)
    new <- nta[!nta %in% nx]
    lnew <- length(new)
    if (lnew > 0) {
      x <- setNames(c(x, rep(0, lnew)), c(nx, new))
      x <- x[match(nta, names(x))]
    }
    x
  })

  ocn_annot5 <- lapply(ocn_annot4, function(x) {
    names(x) <- annot$LUT$LUT_labelname[match(names(x), annot$LUT$LUT_value)]
    na_nx <- is.na(names(x))
    if (sum(na_nx) == 1) names(x)[na_nx] <- "missing label"
    x
  })

  ocn_p1 <- lapply(ocn_annot5, function(x) round(prop.table(x) * 100, 2))
  ocn_p2 <- lapply(ocn_annot5, function(x) round(x / ta * 100, 2))

  ocn_p <- lapply(seq_along(ocn_p1), function(i) {
    x <- data.frame(area = names(ocn_p1[[i]]),
                    size = ocn_annot5[[i]],
                    to_cluster = ocn_p1[[i]],
                    to_area = as.vector(ocn_p2[[i]])
                    )
    x <- x[x$to_cluster != 0, ]
    x[order(x$to_cluster, decreasing = TRUE), ]
  })

  # return
  ocn_p

}

demo_plot <- function(vw, stack = NULL, software = c("tksurfer", "freeview"), verbose = TRUE){
  if (is.null(stack)) stop("`stack` not defined. Please choose: ", paste(stacks(vw), collapse = ", "))
  software <- match.arg(software)
  if (software == "freeview"){
		#freeview --surface $SUBJECTS_DIR/fsaverage/surf/${hemi}.inflated:overlay=./${OBN2}.cache.th13.abs.sig.cluster.mgh:overlay_threshold=1.3,4 ${cmds}
    cmdStr <- paste0("freeview --surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated:overlay=", vw$stack$cluster.mgh[[stack]])
    system(cmdStr, ignore.stdout = !verbose)
  } else if (software == "tksurfer") {
    #cmdStr <- paste("tksurfer", vw$input$target, vw$input$hemi, "inflated", "-overlay", vw$stack$coef[[stack]], "-overlay", vw$stack$cluster.mgh[[stack]])
    cmdStr <- paste("tksurfer", vw$input$target, vw$input$hemi, "inflated", "-overlay", vw$stack$cluster.mgh[[stack]])
    message2(verbose = verbose, cmdStr)
    system(cmdStr, ignore.stdout = !verbose)
  }
}

qsnap <- function(x) paste0("--ss ", x)
qsnap_a <- function(n) paste("--camera Azimuth", n)
qsnap_e <- function(n) paste("--camera Elevation", n)
qsnap_zoom <- function(n) paste("--zoom", n)

qdecr_snap <- function(vw, stack = NULL, ext = ".tiff", zoom = 1.4, plot_brain = TRUE, save_plot = TRUE) {
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
  p
}

summary.vw_fastlm <- function(vw, verbose = FALSE, annot = FALSE, file = "aparc.annot", regions = 3){

  if(!is.numeric(regions)) stop("Provided `annot_length` is not a number.")
  if(!(regions > 0)) stop("Provided `annot_length` is not a positive number.")

  # Retrieve the stacks to be loaded
  s <- vw$stack$names

  # Load the cluster/ocn data
  cl <- list()
  for (i in seq_along(s)) cl[[i]] <- qdecr_read_ocn(vw, i)

  # For each stack, deconstruct which indices each cluster is composed of
  for (i in seq_along(s)) cl[[i]] <- lapply(seq_len(max(cl[[i]]$x)), function(x) which(cl[[i]]$x == x))

  # Make an empty data frame
  nr <- sum(sapply(cl, length))
  cm <- paste0("mean_", vw$input$measure)
  nc <- c("variable", "cluster", "n_vertices", cm, "mean_coefficient", "mean_se")
  cs <- lapply(cl, function(x) {
    y <- as.data.frame(matrix(NA, nrow = length(x), ncol = length(nc)))
    names(y) <- nc
    y})

  # Get info per stack
  cl2 <- list()
  ct_m <- rep(FALSE, length(vw$post$final_mask))
  for (i in seq_along(s)){

    message2(paste0("Summarizing data for `", s[i], "`"), verbose = verbose)

    # Load coef and se
    coo <- qdecr_read_coef(vw, i)$x
    cse <- qdecr_read_se(vw, i)$x

    for (j in seq_along(cl[[i]])){
      ct <- cl[[i]][[j]]

      ct_mask <- ct_m
      ct_mask[ct] <- TRUE

      cs[[i]][j, "variable"] <- s[i]
      cs[[i]][j, "cluster"] <- j
      cs[[i]][j, "n_vertices"] <- length(ct)
      cs[[i]][j, cm] <- mean(vw$post$mgh_description$vertex_mean[ct_mask], na.rm = TRUE)
      cs[[i]][j, "mean_coefficient"] <- mean(coo[ct])
      cs[[i]][j, "mean_se"] <- mean(cse[ct])
    }

  }

  cs2 <- do.call("rbind", cs)

  # Get cluster annotation information
  if (annot) {
    ca <- qdecr_clusters(vw, name = file)

    ca2 <- lapply(ca, function(x) {
      nr2 <- nrow(x)
      if (nr2 > regions) nr2 <- regions
      x[seq_len(nr2), ]
    })

    ca3 <- lapply(ca2, function(y) apply(y, 1, function(x) paste0(x[1], " (", trimws(x[3]), "%, ", trimws(x[4]), "%)")))
    ca4 <- lapply(ca3, function(x) {
      lt <- length(x)
      if (lt < regions) x[(lt+1):regions] <- NA
      x
    })
    ca5 <- do.call("rbind", ca4)
    colnames(ca5) <- paste0("top_region", seq_len(regions))

    cs2 <- cbind(cs2, ca5)
  }

  cs2
}


