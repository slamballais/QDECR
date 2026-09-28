#!/usr/bin/env Rscript
# Puts the example's results into the site:
#
#   Rscript export.R <repository root>
#
# export.sh calls it. From results/ (run.sh) and data/ (download.sh) under
# QDECR_EXAMPLE_ROOT it writes, relative to website/:
#
# - src/data/example/run.json: the run as src/lib/example.ts reads it: the sample, the
#   software, the model, and per hemisphere the stacks and the significant clusters with
#   their size, cluster-wise p-value, peak and regions. The site's tables and figures
#   captions come from here.
# - src/data/example/subjects.csv: the data frame the analysis read.
# - src/data/example/output/: what R printed, as text, for the site's output blocks: the
#   run's log, print(out), stacks(out), summary(out, annot = TRUE), and FreeSurfer's own
#   cluster summary for the age stack.
# - src/assets/example/: the histograms and the qdecr_snap() images as PNG.
# - public/viewer/: for the interactive viewer, the inflated surface
#   and curvature of fsaverage6, and the age stack's t-statistic and cluster maps on it.
#   fsaverage6's 40,962 vertices are the first 40,962 of fsaverage (the icosahedra nest),
#   so a map is downsampled by taking its first 40,962 values; the script checks the
#   nesting on the spheres before relying on it.
#
# Every derived file carries ABIDE's CC BY-NC-SA licence; the credit is
# in run.json for the pages to print.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: export.R <repository root>")
repo <- args[1]
root <- Sys.getenv("QDECR_EXAMPLE_ROOT")
if (!nzchar(root)) stop("QDECR_EXAMPLE_ROOT is not set: run this through export.sh, or source env.sh first.")
fshome <- Sys.getenv("FREESURFER_HOME")

suppressPackageStartupMessages(library(QDECR))

results <- file.path(root, "results")
staged <- file.path(results, "site")
website <- file.path(repo, "website")
data_dir <- file.path(website, "src", "data", "example")
output_dir <- file.path(data_dir, "output")
assets_dir <- file.path(website, "src", "assets", "example")
viewer_dir <- file.path(website, "public", "viewer")
for (d in c(data_dir, output_dir, assets_dir, viewer_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

hemis <- c("lh", "rh")
runs <- lapply(hemis, function(h) jsonlite::read_json(file.path(staged, paste0(h, ".run.json"))))
names(runs) <- hemis

# ---------- the sample ----------

subjects <- read.csv(file.path(root, "data", "phenotypes.csv"))
excluded <- read.csv(file.path(root, "data", "excluded.csv"))
file.copy(file.path(root, "data", "phenotypes.csv"), file.path(data_dir, "subjects.csv"), overwrite = TRUE)

dataset <- list(
  name = "ABIDE I, as preprocessed with FreeSurfer 5.1 by the Preprocessed Connectomes Project",
  site = unique(subjects$site),
  n = nrow(subjects),
  sex = list(female = sum(subjects$sex == "female"), male = sum(subjects$sex == "male")),
  age = list(
    min = round(min(subjects$age), 1),
    max = round(max(subjects$age), 1),
    mean = round(mean(subjects$age), 1),
    median = round(median(subjects$age), 1)
  ),
  excluded = excluded
)

# ---------- the software ----------

os <- {
  release <- readLines("/etc/os-release")
  sub('^PRETTY_NAME="?([^"]*)"?$', "\\1", grep("^PRETTY_NAME=", release, value = TRUE))
}
platform <- if (any(grepl("microsoft", readLines("/proc/version"), ignore.case = TRUE))) {
  "Windows Subsystem for Linux 2"
} else {
  "Linux"
}
software <- list(
  qdecr = runs$lh$qdecr,
  r = runs$lh$r,
  freesurfer = runs$lh$freesurfer,
  os = os,
  platform = platform
)

# ---------- each hemisphere ----------

# FreeSurfer's binary triangle format: a three-byte magic number, a "created by" line
# and an empty one, then the counts and the coordinates, all big-endian. Only the
# vertices are needed here, to check that fsaverage6 nests in fsaverage.
read_vertices <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  stopifnot(identical(as.integer(bytes[1:3]), c(255L, 255L, 254L)))
  newline <- as.raw(10)
  header_end <- which(bytes[-1] == newline & bytes[-length(bytes)] == newline)[1] + 1
  con <- rawConnection(bytes[-seq_len(header_end)])
  on.exit(close(con))
  n_vertices <- readBin(con, "integer", size = 4, endian = "big")
  readBin(con, "integer", size = 4, endian = "big")
  matrix(readBin(con, "numeric", n = n_vertices * 3, size = 4, endian = "big"), ncol = 3, byrow = TRUE)
}

gzip_to <- function(from, to) {
  bytes <- readBin(from, "raw", file.size(from))
  con <- gzfile(to, "wb")
  on.exit(close(con))
  writeBin(bytes, con)
}

# The number of a stack by its name, as qdecr_snap() and the file names count them.
stack_number <- function(out, name) which(stacks(out) == name)

# mri_surfcluster's summary table for a stack: comment lines, then one row per cluster.
# A stack with no cluster has no rows, which read.table reports as an error.
read_cluster_summary <- function(path) {
  columns <- c("cluster", "max", "vtxMax", "sizeMm2", "mniX", "mniY", "mniZ", "cwp", "cwpLow", "cwpHi", "nVtxs", "wghtVtx", "annot")
  rows <- tryCatch(
    read.table(path, comment.char = "#", header = FALSE, stringsAsFactors = FALSE, col.names = columns),
    error = function(e) NULL
  )
  if (is.null(rows)) as.data.frame(setNames(replicate(length(columns), logical(0), simplify = FALSE), columns)) else rows
}

clean_log <- function(lines) {
  # A progress bar redraws its line with carriage returns; keep what it showed last.
  lines <- sub(".*\r", "", lines)
  # Drop the empty run of lines the bar leaves behind.
  lines[!(lines == "" & c(TRUE, lines[-length(lines)] == ""))]
}

export_hemisphere <- function(hemi) {
  run <- runs[[hemi]]
  out <- qdecr_load(file.path(results, run$project))
  n6 <- 40962

  # -- the text the site quotes --
  writeLines(clean_log(readLines(file.path(results, paste0(hemi, ".log.txt")))), file.path(output_dir, paste0(hemi, ".log.txt")))
  for (name in c("print", "stacks", "summary")) {
    file.copy(file.path(staged, paste0(hemi, ".", name, ".txt")), file.path(output_dir, paste0(hemi, ".", name, ".txt")), overwrite = TRUE)
  }
  age <- stack_number(out, "age")
  file.copy(out$stack$cluster.summary[[age]], file.path(output_dir, paste0(hemi, ".age.cluster.summary.txt")), overwrite = TRUE)

  # -- the figures --
  for (png in list.files(staged, pattern = paste0("^", hemi, "\\..*\\.png$"))) {
    file.copy(file.path(staged, png), file.path(assets_dir, png), overwrite = TRUE)
  }

  # -- the clusters --
  # QDECR's summary gives the means over each cluster and the regions it lies in;
  # mri_surfcluster's gives the size, the cluster-wise p-value and the peak. The
  # regions come from the same internal function summary() uses, whole rather than as
  # the strings it prints.
  summary_rows <- summary(out, annot = FALSE)
  regions <- QDECR:::qdecr_clusters(out)
  clusters <- list()
  for (i in seq_len(nrow(summary_rows))) {
    row <- summary_rows[i, ]
    stack <- stack_number(out, row$variable)
    fs <- read_cluster_summary(out$stack$cluster.summary[[stack]])
    fs <- fs[fs$cluster == row$cluster, ]
    stopifnot(nrow(fs) == 1, fs$nVtxs == row$n_vertices)
    top <- head(regions[[i]], 3)
    clusters[[i]] <- list(
      stack = row$variable,
      cluster = row$cluster,
      nVertices = row$n_vertices,
      sizeMm2 = round(fs$sizeMm2, 1),
      cwp = fs$cwp,
      peak = list(value = round(fs$max, 3), vertex = fs$vtxMax, region = fs$annot),
      meanThickness = signif(row$mean_thickness, 4),
      meanCoefficient = signif(row$mean_coefficient, 4),
      meanSe = signif(row$mean_se, 4),
      regions = data.frame(name = top$area, ofCluster = top$to_cluster, ofRegion = top$to_area)
    )
  }

  # -- the viewer's maps, on fsaverage6 --
  sphere7 <- read_vertices(file.path(fshome, "subjects", "fsaverage", "surf", paste0(hemi, ".sphere")))
  sphere6 <- read_vertices(file.path(fshome, "subjects", "fsaverage6", "surf", paste0(hemi, ".sphere")))
  stopifnot(nrow(sphere6) == n6, isTRUE(all.equal(sphere7[seq_len(n6), ], sphere6)))
  for (surface in c("inflated", "curv")) {
    file.copy(
      file.path(fshome, "subjects", "fsaverage6", "surf", paste0(hemi, ".", surface)),
      file.path(viewer_dir, paste0(hemi, ".", surface)),
      overwrite = TRUE
    )
  }
  for (map in c("t", "ocn")) {
    full <- if (map == "t") qdecr_read_t(out, age) else qdecr_read_ocn(out, age)
    temp <- tempfile(fileext = ".mgh")
    save.mgh(as_mgh(full$x[seq_len(n6)]), temp)
    gzip_to(temp, file.path(viewer_dir, paste0(hemi, ".age.", map, ".mgz")))
    unlink(temp)
  }

  describe <- as.data.frame(out$describe$data, stringsAsFactors = FALSE)
  list(
    project = run$project,
    vertices = list(
      loaded = as.integer(describe$value[describe$name == "Vertices loaded"]),
      analysed = sum(out$post$final_mask)
    ),
    fwhmEstimate = as.numeric(out$post$fwhm_est),
    seconds = run$seconds,
    stacks = data.frame(number = seq_along(stacks(out)), name = stacks(out)),
    clusters = clusters,
    # Not in the schema, for the record: the call as QDECR printed it.
    call = out$describe$call[out$describe$call[, "name"] == "qdecr_fastlm call", "value"]
  )
}

hemispheres <- lapply(hemis, export_hemisphere)
names(hemispheres) <- hemis

# ---------- the model ----------

lh <- qdecr_load(file.path(results, runs$lh$project))
mcz <- as.numeric(sub(".*\\.th(\\d+)\\..*", "\\1", lh$stack$cluster.summary[[1]]))
call <- hemispheres$lh$call
cwp <- if (grepl("cwp_thr *= *", call)) as.numeric(sub(".*cwp_thr *= *([0-9.]+).*", "\\1", call)) else 0.025
model <- list(
  formula = paste(deparse(formula(lh)), collapse = ""),
  measure = lh$input$measure,
  fwhm = lh$input$fwhm,
  mczThr = mcz,
  cwpThr = cwp,
  nCores = lh$input$n_cores
)

run <- list(
  date = substr(runs$lh$started, 1, 10),
  dataset = dataset,
  software = software,
  model = model,
  hemispheres = hemispheres,
  credit = list(
    licence = "CC BY-NC-SA 3.0",
    abide = "https://fcon_1000.projects.nitrc.org/indi/abide/",
    pcp = "http://preprocessed-connectomes-project.org/abide/",
    note = "Every figure and number derived from these data carries ABIDE's licence."
  )
)
jsonlite::write_json(run, file.path(data_dir, "run.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

cat(sprintf(
  "Exported: %d subjects, %d clusters (lh), %d clusters (rh); run.json, %d output files, %d figures, %d viewer files.\n",
  dataset$n, length(hemispheres$lh$clusters), length(hemispheres$rh$clusters),
  length(list.files(output_dir)), length(list.files(assets_dir)), length(list.files(viewer_dir))
))
