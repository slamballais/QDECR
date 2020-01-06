#' Combine all MGH files into one object
#'
#' @param input.list path to parent directory
#' @param files.list file names within `input.list`
#'
#' @return matrix with vertex data (columns) for each file (rows)
#'
#'
qdecr_prep_mgh <- function(input_path, 
                           files_list = list.files(input_path),
                           mask, hemi, measure, fwhmc, target, 
                           n_cores, dir_tmp, project, backing, verbose) {
  measure2 <- measure
  if(measure2 == "w_g.pct") measure2 <- "w-g.pct"
  new_files <- file.path(input_path, files_list, "surf", paste(hemi, measure, fwhmc, target, "mgh", sep = "."))
  n <- length(new_files)

  temp <- load.mgh(new_files[1])

  m <- bigstatsr::FBM(temp$ndim1, n, backingfile = backing)
  cl <- if(verbose) parallel::makeForkCluster(n_cores, outfile = "") else parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  capture.output(pb <- txtProgressBar(0, n, style = 3), file = "/dev/null")
  
  foreach(i = seq_len(n)) %dopar% {
    setTxtProgressBar(pb, i)
    m[,i] <- load.mgh(new_files[i])$x
    NULL
  }
  parallel::stopCluster(cl)
  
  if (length(mask) != nrow(m)) stop("Length of mask does not equal number of vertices from data.")
  m
}



#' Load an MGH file into memory
#'
#' @param input.file full path to the mgh file
#'
#' @return mgh object with data and various header elements
#'
#' @export
#' 
load.mgh <- function(input.file) {
  # written by Heath Pardoe, heath.pardoe at nyumc.org, 09/12/2013
  to.read <- file(input.file, "rb")
  v <- readBin(to.read, integer(), endian = "big")
  ndim1 <- readBin(to.read, integer(), endian = "big")
  ndim2 <- readBin(to.read, integer(), endian = "big")
  ndim3 <- readBin(to.read, integer(), endian = "big")
  nframes <- readBin(to.read, integer(), endian = "big")
  type <- readBin(to.read, integer(), endian = "big")
  dof <- readBin(to.read, integer(), endian = "big")
  close(to.read)

  to.read <- file(input.file,"rb")
  dump <- readBin(to.read, double(), size = 4, n = 71, endian = "big")
  #THIS SEEMED TO ONLY READ IN THE FIRST CHUNK OF DATA, and ignored if additional subjects were merged in
  #added ndim1*nframes
  x <- readBin(to.read ,double(), size = 4, n = ndim1*nframes, endian = "big")
  close(to.read)
  list(x = x, v = v, ndim1 = ndim1, ndim2 = ndim2, ndim3 = ndim3, nframes =
         nframes, type = type, dof = dof)
}

#' Load an .annot file into memory
#'
#' @param input.file full path to the .annot file
#'
#' @return annot object with data and various header elements
#'
#' @export
#' 
load.annot <- function(input.file) {
  # rewrite of load.mgh
  to.read <- file(input.file, "rb")
  on.exit (close(to.read))
  a_vtxct <- readBin(to.read, size = 4, integer(), endian = "big")
  a_vd <- readBin(to.read, size = 4, integer(), n = a_vtxct * 2, endian = "big")
  a_tag <- readBin(to.read, size = 4, integer(), endian = "big")
  a_ctabversion <- readBin(to.read, size = 4, integer(), endian = "big")
  a_maxstruc <- readBin(to.read, size = 4, integer(), endian = "big")
  a_len <- readBin(to.read, size = 4, integer(), endian = "big")
  a_fname <- readChar(to.read, a_len)
  a_num_entries <- readBin(to.read, size = 4, integer(), endian = "big")
  
  LUT_Label <- LUT_len <- LUT_labelname <- LUT_red <- LUT_green <- LUT_blue <- LUT_transp <- rep(NA, a_num_entries)
  for (i in seq_len(a_num_entries)) {
    LUT_Label[i] <- readBin(to.read, size = 4, integer(), endian = "big")
    LUT_len[i] <- readBin(to.read, size = 4, integer(), endian = "big")
    LUT_labelname[i] <- readChar(to.read, LUT_len[i])
    LUT_red[i] <- readBin(to.read, size = 4, integer(), endian = "big")
    LUT_green[i] <- readBin(to.read, size = 4, integer(), endian = "big")
    LUT_blue[i] <- readBin(to.read, size = 4, integer(), endian = "big")
    LUT_transp[i] <- readBin(to.read, size = 4, integer(), endian = "big")
  }
  close(to.read)
  on.exit()
  
  a_vtxct2 <- seq(1, a_vtxct * 2, 2)
  a_vd_vno <- a_vd[a_vtxct2]
  a_vd_label <- a_vd[a_vtxct2 + 1]

  list(vtxct = a_vtxct, 
       vd_vno = a_vd_vno, 
       vd_label = a_vd_label,
       tag = a_tag,
       ctabversion = a_ctabversion,
       maxstruc = a_maxstruc,
       len = a_len,
       fname = a_fname,
       num_entries = a_num_entries,
       LUT = data.frame(LUT_Label = LUT_Label,
                        LUT_value = (LUT_blue * 256 ^ 2) + (LUT_green * 256) + LUT_red,
                        LUT_len = LUT_len,
                        LUT_labelname = LUT_labelname,
                        LUT_red = LUT_red,
                        LUT_green = LUT_green,
                        LUT_blue = LUT_blue,
                        LUT_transp = LUT_transp,
                        stringsAsFactors = FALSE)
       )
}

#' Save out an MGH file from a big matrix
#'
#' @param vol MGH object (as from load.mgh)
#' @param fname file name to be used to save out the data
#'
#' @export
#' 
bsfbm2mgh <-function(fbm, fname, filter = NULL) {

  # R translation of save_mgh.m
  # written by Heath Pardoe, heath.pardoe at nyumc.org, 09/12/2013
  # modified by Ryan Muetzel to handle nframes/stacked data
  # modified by Sander Lamballais to convert .bk files (from bigstatsr) to .mgh

  MRI.UCHAR <-  0
  MRI.INT <-    1
  MRI.LONG <-   2
  MRI.FLOAT <-  3
  MRI.SHORT <-  4
  MRI.BITMAP <- 5
  MRI.TENSOR <- 6
  slices <- c(1:256)

  fid <- file(fname, open = "wb", blocking = TRUE)

  if (!is.null(filter)){
    if (max(filter) > fbm$nrow) stop("In bsfbm2mgh, the defined filter exceeds the bounds.")
    if (min(filter) < 1) stop("In bsfbm2mgh, the defined filter contains non-positive numbers.")
  } else {
    filter <- seq_len(fbm$nrow)
  }

  width <- fbm$ncol
  height <- 1
  depth <- 1
  nframes <- length(filter)

  writeBin(as.integer(1), fid, size = 4, endian = "big")
  writeBin(as.integer(width), fid, size = 4, endian = "big")
  writeBin(as.integer(height), fid, size = 4, endian = "big")
  writeBin(as.integer(depth), fid, size = 4, endian = "big")
  writeBin(as.integer(nframes), fid, size = 4, endian = "big")
  writeBin(as.integer(MRI.FLOAT), fid, size = 4, endian = "big")
  writeBin(as.integer(1), fid, size = 4, endian = "big")
  
  UNUSED.SPACE.SIZE <- 256
  USED.SPACE.SIZE <- (3 * 4 + 4 * 3 * 4)
  unused.space.size <- UNUSED.SPACE.SIZE - 2

  writeBin(as.integer(0), fid, size = 2, endian = "big")
  writeBin(as.integer(rep.int(0, unused.space.size)), fid, size = 1)
  bpv <- 4
  nelts <- width * height
  for (i in filter) writeBin(fbm[i, ], fid, size = 4, endian = "big")
  close(fid)
  NULL
}

#' Save out an MGH file from memory
#'
#' @param vol MGH object (as from load.mgh)
#' @param fname file name to be used to save out the data
#'
#' @export
#'
save.mgh <-function(vol,fname) {

  # R translation of save_mgh.m
  # written by Heath Pardoe, heath.pardoe at nyumc.org, 09/12/2013

  MRI.UCHAR <-  0
  MRI.INT <-    1
  MRI.LONG <-   2
  MRI.FLOAT <-  3
  MRI.SHORT <-  4
  MRI.BITMAP <- 5
  MRI.TENSOR <- 6
  slices <- c(1:256)

  fid <- file(fname, open = "wb", blocking = TRUE)
  on.exit(close(fid))

  width <- vol$ndim1
  height <- vol$ndim2
  depth <- vol$ndim3
  #RLM added
  nframes <- vol$nframes

  writeBin(as.integer(1), fid, size = 4, endian = "big")
  writeBin(as.integer(width), fid, size = 4, endian = "big")
  writeBin(as.integer(height), fid, size = 4, endian = "big")
  writeBin(as.integer(depth), fid, size = 4, endian = "big")
  #here we replace the default of 1 frame
  writeBin(as.integer(nframes), fid, size = 4, endian = "big")
  #writeBin(as.integer(1),fid,size = 4, endian = "big")

  # HP note: I've ignored all the diffusion tensor stuff
  writeBin(as.integer(MRI.FLOAT), fid, size = 4, endian = "big")

  writeBin(as.integer(1), fid, size = 4, endian = "big")
  # dof = fread(fid, 1, 'int');
  ## HP note: ignored ^this line^ from save_mgh.m

  UNUSED.SPACE.SIZE <- 256
  USED.SPACE.SIZE <- (3 * 4 + 4 * 3 * 4)  # space for ras transform

  unused.space.size <- UNUSED.SPACE.SIZE - 2

  # ignored all the stuff about "M" - could probably do it if necessary so let me know
  #if (nargin > 2)
  ## fwrite(fid, 1, 'short')        # ras.good.flag <- 0
  # writeBin(1,fid,size = 2, endian = "big")
  # unused.space.size <- unused.space.size - USED.SPACE.SIZE
  ## fwrite(fid, sizes(1), 'float32')  # xsize
  ## fwrite(fid, sizes(2), 'float32')  # ysize
  ## fwrite(fid, sizes(3), 'float32')  # zsize
  #
  # fwrite(fid, M(1,1), 'float32')   # x.r
  # fwrite(fid, M(2,1), 'float32')   # x.a
  # fwrite(fid, M(3,1), 'float32')   # x.s

  # fwrite(fid, M(1,2), 'float32')   # y.r
  # fwrite(fid, M(2,2), 'float32')   # y.a
  # fwrite(fid, M(3,2), 'float32')   # y.s

  # fwrite(fid, M(1,3), 'float32')   # z.r
  # fwrite(fid, M(2,3), 'float32')   # z.a
  # fwrite(fid, M(3,3), 'float32')   # z.s

  # fwrite(fid, M(1,4), 'float32')   # c.r
  # fwrite(fid, M(2,4), 'float32')   # c.a
  # fwrite(fid, M(3,4), 'float32')   # c.s
  #else
  # fwrite(fid, 0, 'short')        # ras.good.flag <- 0
  writeBin(as.integer(0), fid, size = 2, endian = "big")

  #   } #

  writeBin(as.integer(rep.int(0, unused.space.size)), fid, size = 1)
  bpv <- 4    # bytes/voxel
  nelts <- width * height   # bytes per slice
  #writeBin(vol$x, fid, size = 4, endian = "big")
  writeBin(vol$x, fid, size = 4, endian = "big")
}


#' Create an object that is structured like an mgh object
#'
#' @param x the vertex-wise values
#' @param v (to be added)
#' @param ndim1 (to be added)
#' @param ndim2 (to be added)
#' @param ndim3 (to be added)
#' @param nframes (to be added)
#' @param type (to be added)
#' @param dof (to be added)
#'
#' @export
#'

as_mgh <- function(x, v = NULL, ndim1 = NULL, ndim2 = NULL, ndim3 = NULL, nframes, type, dof) {
  
  if (is.vector(x)) {
    v <- 1L
    ndim1 <- as.integer(length(x))
    ndim2 <- 1L
    ndim3 <- 1L
    type <- 3L
    nframes <- 1L
    dof <- 0L
  } else {
    stop("as_mgh only support objects of class `vector` right now.")
  }
  
  out <- list(x = x, 
              v = v, 
              ndim1 = ndim1, 
              ndim2 = ndim2, 
              ndim3 = ndim3, 
              nframes = nframes, 
              type = type, 
              dof = dof)
  
  return(out)
  
}