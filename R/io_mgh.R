#' Combine all MGH files into one object
#'
#' @param input_path path to parent directory
#' @param files_list file names within `input.list`
#' @param mask mask to be used in later steps (i.e. analysis)
#' @param hemi "lh" or "rh"
#' @param measure vertex measure (e.g. "thickness")
#' @param fwhmc the fwhm value in full notation (i.e. 20, 23, 30, etc)
#' @param target the target template (usually "fsaverage")
#' @param n_cores the number of cores to be used
#' @param dir_tmp directory to store the temporary big matrices; useful for shared memory; defaults to `dir_out`
#' @param project the base name you want to assign to the output files
#' @param backing the path to the backing file
#' @param verbose if TRUE, writes out standard log; if FALSE, no output is generated
#' @return matrix with vertex data (columns) for each file (rows)
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
  utils::capture.output(pb <- utils::txtProgressBar(0, n, style = 3), file = "/dev/null")
  i <- NULL
  foreach(i = seq_len(n)) %dopar% {
    utils::setTxtProgressBar(pb, i)
    m[,i] <- load.mgh(new_files[i])$x
    NULL
  }
  parallel::stopCluster(cl)
  if (length(mask) != nrow(m)) stop("Length of mask does not equal number of vertices from data.")
  m
}

#' Load an MGH file into memory
#' 
#' Originally written by Heath Perdoe (09/12/2013)
#'
#' @param input.file full path to the mgh file
#'
#' @return mgh object with data and various header elements
#'
#' @export
#' 
load.mgh <- function(input.file) {
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
#' Modified version of save_mgh.m (by Heath Pardoe, 09/12/2013)
#'
#' @param fbm a file backed matrix
#' @param fname file name to be used to save out the data
#' @param filter a vector of indices of rows to include of fbm
#'
#' @export
#' 
bsfbm2mgh <-function(fbm, fname, filter = NULL) {
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
#' #' R translation of save_mgh.m (by Heath Pardoe, 09/12/2013)
#'
#' @param vol MGH object (as from load.mgh)
#' @param fname file name to be used to save out the data
#'
#' @export
#'
save.mgh <-function(vol,fname) {
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
  nframes <- vol$nframes
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
  writeBin(vol$x, fid, size = 4, endian = "big")
}

#' Create an object that is structured like an mgh object
#'
#' @param x the vertex-wise values
#' @param v Version (default = 1)
#' @param ndim1 Width / 1st dimension
#' @param ndim2 Height / 2nd dimension
#' @param ndim3 Depth / 3rd dimension
#' @param nframes Number of scalar components
#' @param type Data type, can be UCHAR (0), SHORT (4), INT (1) or FLOAT (3)
#' @param dof Degrees of freedom
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