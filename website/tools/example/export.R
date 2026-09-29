#!/usr/bin/env Rscript
# Puts the example's results into the site:
#
#   Rscript export.R <repository root>
#
# export.sh calls it. From results/ (run.sh) and data/ (download.sh) under
# QDECR_EXAMPLE_ROOT it writes, relative to website/:
#
# - src/data/example/run.json: the run as src/lib/example.ts reads it, and nothing the
#   schema does not name, since it is strict: the sample, the software, the model, the
#   poster's colour scale, the credit, and per hemisphere the stacks and the significant clusters with their size,
#   cluster-wise p-value, peak and regions. The site's tables and figure captions come
#   from here. The names are the glossary's (src/data/glossary.md).
# - src/data/example/subjects.csv: the data frame the analysis read.
# - src/data/example/output/: what R printed, as text, for the site's output blocks: the
#   run's log, head() of the data frame, print(out), stacks(out), summary(out, annot =
#   TRUE), and FreeSurfer's own cluster summary for the age stack.
# - src/assets/example/: the histograms and the qdecr_snap() images as PNG.
# - public/viewer/: for the home page's viewer, per hemisphere the
#   inflated surface of fsaverage6, where its sulci are, and the poster's map on it (the
#   age stack's -log10(p) on its significant clusters), as gzipped MZ3, the compact format
#   NiiVue reads. fsaverage6's 40,962 vertices are the first 40,962 of fsaverage (the
#   icosahedra nest), so a map is downsampled by taking its first 40,962 values; the
#   script checks the nesting on the spheres before relying on it.
#
# Every derived file carries ABIDE's CC BY-NC-SA licence: the credit is
# in run.json for the pages to print, and as LICENCE.txt beside the figures and the
# viewer's files.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: export.R <repository root>")
repo <- args[1]
root <- Sys.getenv("QDECR_EXAMPLE_ROOT")
if (!nzchar(root)) stop("QDECR_EXAMPLE_ROOT is not set: run this through export.sh, or source env.sh first.")
fshome <- Sys.getenv("FREESURFER_HOME")
if (!nzchar(fshome)) stop("FREESURFER_HOME is not set: run this through export.sh, or source env.sh first.")

suppressPackageStartupMessages(library(QDECR))

results <- file.path(root, "results")
# What run.R staged for the site.
staged <- file.path(results, "site")
website <- file.path(repo, "website")
data_dir <- file.path(website, "src", "data", "example")
output_dir <- file.path(data_dir, "output")
assets_dir <- file.path(website, "src", "assets", "example")
viewer_dir <- file.path(website, "public", "viewer")
for (d in c(data_dir, output_dir, assets_dir, viewer_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
# Everything in public/viewer/ is the export's, so it starts empty: a file the export
# no longer writes would otherwise stay in the site.
unlink(list.files(viewer_dir, full.names = TRUE))

hemis <- c("lh", "rh")
runs <- lapply(hemis, function(h) jsonlite::read_json(file.path(staged, paste0(h, ".run.json"))))
names(runs) <- hemis

# ---------- the credit ----------
# ABIDE I's terms: non-commercial research use under CC BY-NC-SA, the dataset named,
# and its funding acknowledged; the PCP asks for its abstract to be cited.

credit <- list(
  licence = "CC BY-NC-SA 3.0",
  abide = list(
    url = "https://fcon_1000.projects.nitrc.org/indi/abide/",
    cite = paste(
      "Di Martino A, Yan C-G, Li Q, et al. (2014). The autism brain imaging data exchange:",
      "towards a large-scale evaluation of the intrinsic brain architecture in autism.",
      "Molecular Psychiatry 19, 659-667. https://doi.org/10.1038/mp.2013.78"
    )
  ),
  pcp = list(
    # http, not https: the site's certificate is issued for another host, so an https link
    # shows the reader a certificate error (checked 29 Sep 2026).
    url = "http://preprocessed-connectomes-project.org/abide/",
    cite = paste(
      "Craddock C, Benhajali Y, Chu C, et al. (2013). The Neuro Bureau Preprocessing Initiative:",
      "open sharing of preprocessed neuroimaging data and derivatives.",
      "Frontiers in Neuroinformatics, Neuroinformatics 2013. https://doi.org/10.3389/conf.fninf.2013.09.00041"
    )
  ),
  funding = paste(
    "Primary support for the work by Adriana Di Martino was provided by the NIMH (K23MH087770)",
    "and the Leon Levy Foundation. Primary support for the work by Michael P. Milham and the INDI",
    "team was provided by gifts from Joseph P. Healy and the Stavros Niarchos Foundation to the",
    "Child Mind Institute, as well as by an NIMH award to MPM (R03MH096321)."
  )
)

credit_text <- c(
  "The files in this directory are derived from ABIDE I, as preprocessed by the",
  "Preprocessed Connectomes Project, and carry that data's licence, CC BY-NC-SA 3.0,",
  "whatever the site's own.",
  "",
  paste0("ABIDE: ", credit$abide$url),
  paste0("  ", credit$abide$cite),
  paste0("PCP: ", credit$pcp$url),
  paste0("  ", credit$pcp$cite),
  "",
  credit$funding
)
for (d in c(assets_dir, viewer_dir)) writeLines(credit_text, file.path(d, "LICENCE.txt"))

# ---------- the sample ----------

subjects <- read.csv(file.path(root, "data", "phenotypes.csv"))
excluded <- read.csv(file.path(root, "data", "excluded.csv"))
invisible(file.copy(file.path(root, "data", "phenotypes.csv"), file.path(data_dir, "subjects.csv"), overwrite = TRUE))

# The data frame's first rows as the quick start prints them, read the way run.R reads it.
pheno <- subjects
pheno$sex <- factor(pheno$sex, levels = c("female", "male"))
writeLines(capture.output(head(pheno)), file.path(output_dir, "pheno.head.txt"))

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
# and an empty one, then the counts, the coordinates and the triangles (by vertex, from
# 0), all big-endian.
read_surface <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  stopifnot(identical(as.integer(bytes[1:3]), c(255L, 255L, 254L)))
  newline <- as.raw(10)
  header_end <- which(bytes[-1] == newline & bytes[-length(bytes)] == newline)[1] + 1
  con <- rawConnection(bytes[-seq_len(header_end)])
  on.exit(close(con))
  n_vertices <- readBin(con, "integer", size = 4, endian = "big")
  n_faces <- readBin(con, "integer", size = 4, endian = "big")
  list(
    vertices = matrix(readBin(con, "numeric", n = n_vertices * 3, size = 4, endian = "big"), ncol = 3, byrow = TRUE),
    faces = matrix(readBin(con, "integer", n = n_faces * 3, size = 4, endian = "big"), ncol = 3, byrow = TRUE)
  )
}

# FreeSurfer's curvature format: a three-byte magic number, the counts of vertices, of
# faces and of values per vertex (1), then a big-endian float per vertex.
read_curvature <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  stopifnot(identical(as.integer(readBin(con, "raw", 3)), c(255L, 255L, 255L)))
  counts <- readBin(con, "integer", n = 3, size = 4, endian = "big")
  stopifnot(counts[3] == 1)
  readBin(con, "numeric", n = counts[1], size = 4, endian = "big")
}

# MZ3, Surf Ice's mesh format, which NiiVue reads: a 16-byte header (the magic number
# "MZ", a bit field of what follows, the counts of triangles and vertices, and a count of
# bytes to skip, 0), then the triangles as 0-based int32 triples, the vertices as float32
# triples, and one float32 value per vertex, each part optional and all little-endian,
# the whole gzipped. A file of values alone is a map for a mesh loaded from another.
write_mz3 <- function(path, n_vertices, faces = NULL, vertices = NULL, values = NULL) {
  stopifnot(is.null(vertices) || nrow(vertices) == n_vertices, is.null(values) || length(values) == n_vertices)
  # The bits for triangles (1), vertices (2) and values (8). Summed from a mask: in R, `!`
  # binds more loosely than `*` and `+`, so 1L * !is.null(faces) + ... would negate the rest.
  attributes <- sum(c(1L, 2L, 8L)[c(!is.null(faces), !is.null(vertices), !is.null(values))])
  con <- gzfile(path, "wb")
  on.exit(close(con))
  writeBin(c(23117L, attributes), con, size = 2, endian = "little")
  writeBin(c(if (is.null(faces)) 0L else nrow(faces), as.integer(n_vertices), 0L), con, size = 4, endian = "little")
  if (!is.null(faces)) writeBin(as.integer(t(faces)), con, size = 4, endian = "little")
  if (!is.null(vertices)) writeBin(as.numeric(t(vertices)), con, size = 4, endian = "little")
  if (!is.null(values)) writeBin(as.numeric(values), con, size = 4, endian = "little")
}

# The number of a stack by its name, as qdecr_snap() and the file names count them.
stack_number <- function(out, name) which(stacks(out) == name)

# mri_surfcluster's summary table for a stack: comment lines, then one row per cluster,
# none for a stack without clusters. Only a malformed row is an error.
cluster_summary_columns <- c(
  "cluster", "max", "vtxMax", "sizeMm2", "mniX", "mniY", "mniZ",
  "cwp", "cwpLow", "cwpHi", "nVtxs", "wghtVtx", "annot"
)
read_cluster_summary <- function(path) {
  lines <- readLines(path)
  rows <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  if (length(rows) == 0) {
    empty <- as.data.frame(setNames(replicate(length(cluster_summary_columns), logical(0), simplify = FALSE), cluster_summary_columns))
    return(empty)
  }
  read.table(text = rows, header = FALSE, stringsAsFactors = FALSE, col.names = cluster_summary_columns)
}

# The home page's poster: the age stack's -log10(p) on its
# significant clusters, in Freeview's heat colours, which match the site's orange; the
# map has no sign, so the page says the effect is thinning. Drawn by Freeview at a size
# of the site's choosing rather than its default window, twice over, and trimmed to the
# brain. The colours run from the cluster-forming threshold, p = 0.001 (3 on this scale),
# to p = 1e-10 (10), where the peaks saturate; Sander picked that range by eye over the
# full one, which leaves the map almost all red. Lateral and medial views, as
# <hemi>.age.p.<view>.png beside the figures. Needs a display, which export.sh provides.
hero_scale <- c(3, 10)
# The map the poster draws and the viewer shows: -log10(p), 0 off the clusters.
age_on_clusters <- function(out, age) {
  p <- qdecr_read_p(out, age)
  p$x[!qdecr_read_ocn_mask(out, age)] <- 0
  p
}
render_hero <- function(hemi, out, age) {
  p <- age_on_clusters(out, age)
  overlay <- tempfile(fileext = ".mgh")
  commands <- tempfile(fileext = ".txt")
  on.exit(unlink(c(overlay, commands)))
  save.mgh(p, overlay)
  surface <- sprintf(
    "%s/fsaverage/surf/%s.inflated:overlay=%s:overlay_method=linearopaque:overlay_threshold=%s,%s",
    Sys.getenv("SUBJECTS_DIR"), hemi, overlay, hero_scale[1], hero_scale[2]
  )
  # Freeview's first view is the lateral side of the left hemisphere and the medial side
  # of the right; a half turn shows the other, as in qdecr_snap().
  views <- if (hemi == "lh") c("lateral", "medial") else c("medial", "lateral")
  shot <- function(view) sprintf("--ss %s 2 1", file.path(assets_dir, sprintf("%s.age.p.%s.png", hemi, view)))
  writeLines(c("--viewport 3d", "--viewsize 1200 900", "--zoom 1", shot(views[1]), "--camera Azimuth 180", shot(views[2]), "--quit"), commands)
  status <- system2("freeview", c("--surface", shQuote(surface), "-cmd", commands), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("Freeview failed to draw the ", hemi, " hero images (exit ", status, ")")
}

clean_log <- function(path) {
  # A progress bar redraws its line with carriage returns; keep what it showed last.
  # Read whole and split on newlines only: readLines would take each carriage return
  # as a line ending and keep every redraw.
  lines <- strsplit(readChar(path, file.size(path), useBytes = TRUE), "\n", fixed = TRUE)[[1]]
  lines <- sub(".*\r", "", lines)
  # Drop the empty run of lines the bar leaves behind.
  lines[!(lines == "" & c(TRUE, lines[-length(lines)] == ""))]
}

export_hemisphere <- function(hemi) {
  run <- runs[[hemi]]
  out <- qdecr_load(file.path(results, run$project))
  n6 <- 40962

  # -- the text the site quotes --
  writeLines(clean_log(file.path(results, paste0(hemi, ".log.txt"))), file.path(output_dir, paste0(hemi, ".log.txt")))
  for (name in c("print", "stacks", "summary")) {
    file.copy(file.path(staged, paste0(hemi, ".", name, ".txt")), file.path(output_dir, paste0(hemi, ".", name, ".txt")), overwrite = TRUE)
  }
  age <- stack_number(out, "age")
  file.copy(out$stack$cluster.summary[[age]], file.path(output_dir, paste0(hemi, ".age.cluster.summary.txt")), overwrite = TRUE)

  # -- the figures --
  for (png in list.files(staged, pattern = paste0("^", hemi, "\\..*\\.png$"))) {
    file.copy(file.path(staged, png), file.path(assets_dir, png), overwrite = TRUE)
  }
  render_hero(hemi, out, age)

  # -- the clusters --
  # QDECR's summary gives the means over each cluster and the regions it lies in;
  # mri_surfcluster's gives the size, the cluster-wise p-value and the peak. The
  # regions come from the same internal function summary() uses, whole rather than as
  # the strings it prints. Both list the clusters in the same order: stack by stack,
  # numbered as the cluster map numbers them.
  summary_rows <- summary(out, annot = FALSE)
  regions <- QDECR:::qdecr_clusters(out)
  clusters <- list()
  for (i in seq_len(nrow(summary_rows))) {
    row <- summary_rows[i, ]
    stack <- stack_number(out, row$variable)
    summary_fs <- read_cluster_summary(out$stack$cluster.summary[[stack]])
    summary_fs <- summary_fs[summary_fs$cluster == row$cluster, ]
    stopifnot(nrow(summary_fs) == 1, summary_fs$nVtxs == row$n_vertices)
    top <- head(regions[[i]], 3)
    clusters[[i]] <- list(
      stack = row$variable,
      cluster = row$cluster,
      nVertices = row$n_vertices,
      sizeMm2 = round(summary_fs$sizeMm2, 1),
      clusterwiseP = summary_fs$cwp,
      # The peak's -log10(p) is infinite where p underflowed to zero, as it does for the
      # intercept (thickness is never zero); JSON has no infinity, so that is null.
      peak = list(
        value = if (is.finite(summary_fs$max)) round(summary_fs$max, 3) else NA_real_,
        vertex = summary_fs$vtxMax,
        region = summary_fs$annot
      ),
      meanThickness = signif(row$mean_thickness, 4),
      meanCoefficient = signif(row$mean_coefficient, 4),
      meanSe = signif(row$mean_se, 4),
      regions = data.frame(name = top$area, ofCluster = top$to_cluster, ofRegion = top$to_area)
    )
  }

  # -- the viewer's files, on fsaverage6 --
  surf <- function(subject, name) file.path(fshome, "subjects", subject, "surf", paste0(hemi, ".", name))
  sphere7 <- read_surface(surf("fsaverage", "sphere"))$vertices
  sphere6 <- read_surface(surf("fsaverage6", "sphere"))$vertices
  stopifnot(nrow(sphere6) == n6, isTRUE(all.equal(sphere7[seq_len(n6), ], sphere6)))
  inflated <- read_surface(surf("fsaverage6", "inflated"))
  write_mz3(file.path(viewer_dir, paste0(hemi, ".inflated.mz3")), n6, faces = inflated$faces, vertices = inflated$vertices)
  # Freeview draws the folds in two greys, by the sign of the curvature (positive in a
  # sulcus), and so does the viewer, so the sign is all it needs: 1 in a sulcus, 0 on a
  # gyrus. A few kilobytes gzipped, where the curvature itself is 160.
  sulci <- as.numeric(read_curvature(surf("fsaverage6", "curv")) > 0)
  write_mz3(file.path(viewer_dir, paste0(hemi, ".sulci.mz3")), n6, values = sulci)
  write_mz3(file.path(viewer_dir, paste0(hemi, ".age.p.mz3")), n6, values = age_on_clusters(out, age)$x[seq_len(n6)])

  describe <- as.data.frame(out$describe$data, stringsAsFactors = FALSE)
  list(
    project = run$project,
    vertices = list(
      loaded = as.integer(describe$value[describe$name == "Vertices loaded"]),
      analysed = sum(out$post$final_mask)
    ),
    smoothness = as.numeric(out$post$fwhm_est),
    seconds = run$seconds,
    stacks = data.frame(number = seq_along(stacks(out)), name = stacks(out)),
    clusters = clusters
  )
}

hemispheres <- lapply(hemis, export_hemisphere)
names(hemispheres) <- hemis

# ---------- the model ----------

lh <- qdecr_load(file.path(results, runs$lh$project))

# FreeSurfer codes the cluster-forming threshold in the file names as th13 to th40, the
# p-values QDECR accepts for mcz_thr (R/qdecr_check.R).
cluster_forming_thresholds <- c("13" = 0.05, "20" = 0.01, "23" = 0.005, "30" = 0.001, "33" = 0.0005, "40" = 0.0001)
threshold_code <- sub(".*\\.th(\\d+)\\..*", "\\1", lh$stack$cluster.summary[[1]])
cluster_forming <- cluster_forming_thresholds[[threshold_code]]

# The result does not store cwp_thr: QDECR keeps it only in the call, which the result
# holds as text. So it is read back from there: a number, or the default when the call
# left it out. A variable in its place is not a number, and that stops the export rather
# than record a guess.
call_text <- lh$describe$call[lh$describe$call[, "name"] == "qdecr_fastlm call", "value"]
clusterwise <- if (grepl("cwp_thr *= *", call_text)) {
  suppressWarnings(as.numeric(sub(".*cwp_thr *= *([^,)]+).*", "\\1", call_text)))
} else {
  0.025
}
if (is.na(clusterwise)) stop("cwp_thr in the call is not a number: ", call_text)

model <- list(
  formula = paste(deparse(formula(lh)), collapse = ""),
  measure = lh$input$measure,
  fwhm = lh$input$fwhm,
  clusterFormingThreshold = cluster_forming,
  clusterwiseThreshold = clusterwise,
  nCores = lh$input$n_cores
)

run <- list(
  date = substr(runs$lh$started, 1, 10),
  dataset = dataset,
  software = software,
  model = model,
  hemispheres = hemispheres,
  poster = list(scale = list(from = hero_scale[1], to = hero_scale[2])),
  credit = credit
)
jsonlite::write_json(run, file.path(data_dir, "run.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")

cat(sprintf(
  "Exported: %d subjects, %d clusters (lh), %d clusters (rh); run.json, %d output files, %d figures, %d viewer files.\n",
  dataset$n, length(hemispheres$lh$clusters), length(hemispheres$rh$clusters),
  length(list.files(output_dir)), length(list.files(assets_dir)), length(list.files(viewer_dir))
))
