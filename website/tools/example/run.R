#!/usr/bin/env Rscript
# One hemisphere of the example analysis, and what the site shows of it:
#
#   Rscript run.R lh
#
# run.sh calls this for both hemispheres, with a display for Freeview, and keeps the log.
# It fits qdecr_thickness ~ age + sex on the subjects download.sh fetched, then does what
# a reader of the tutorials does with the result: prints it, summarises its clusters,
# draws the histograms and takes the snapshots. What R prints is saved as text, and the
# figures as PNG, under results/site/ for export.R to pick up.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% c("lh", "rh")) stop("usage: run.R lh|rh")
hemi <- args[1]

root <- Sys.getenv("QDECR_EXAMPLE_ROOT")
if (!nzchar(root)) stop("QDECR_EXAMPLE_ROOT is not set: run this through run.sh, or source env.sh first.")
results <- file.path(root, "results")
# What the site takes is staged in results/site/, for export.R.
staged <- file.path(results, "site")
dir.create(staged, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(QDECR))
# The cluster summary is a wide table; keep its rows whole in the saved text.
options(width = 200)

# ---------- the data ----------
# The data frame as Get started describes it: one row per subject, the id column holding
# the subject's directory in SUBJECTS_DIR, sex as a factor with female as the reference,
# no missing values.
pheno <- read.csv(file.path(root, "data", "phenotypes.csv"))
pheno$sex <- factor(pheno$sex, levels = c("female", "male"))
stopifnot(!anyNA(pheno[, c("id", "age", "sex")]))

# ---------- the analysis ----------
# The call the quick start shows. n_cores = 4: the machine had 16; the site does not
# need more, and this keeps the example from claiming a large machine. The results go
# to results/<hemi>.age_sex.thickness/; a rerun overwrites them.
started <- Sys.time()
out <- qdecr_fastlm(
  qdecr_thickness ~ age + sex,
  data = pheno,
  id = "id",
  hemi = hemi,
  project = "age_sex",
  dir_out = results,
  dir_tmp = "/dev/shm",
  n_cores = 4,
  clobber = TRUE
)
seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

# ---------- what R prints ----------
# print() writes through message(), so it is captured from the message stream.
save_text <- function(lines, name) writeLines(lines, file.path(staged, paste0(hemi, ".", name, ".txt")))
save_text(capture.output(print(out), type = "message"), "print")
save_text(capture.output(print(stacks(out))), "stacks")
clusters <- summary(out, annot = TRUE)
save_text(capture.output(print(clusters, row.names = FALSE)), "summary")
write.csv(clusters, file.path(staged, paste0(hemi, ".summary.csv")), row.names = FALSE)

# ---------- the histograms ----------
# hist(out) as a reader sees it, at print resolution.
figure <- function(name, draw) {
  png(file.path(staged, paste0(hemi, ".", name, ".png")), width = 1800, height = 1200, res = 220, type = "cairo")
  on.exit(dev.off())
  draw()
}
figure("hist-vertex", function() hist(out))
figure("hist-subject", function() hist(out, qtype = "subject"))

# ---------- the snapshots ----------
# qdecr_snap() opens Freeview on the inflated surface with the stack's map on its
# significant clusters, screenshots four views and composes them into one image, which
# it writes as TIFF next to the output directory. The site takes the composed image as
# PNG. The age stack's snapshots are the site's main figures, so anything wrong with
# them stops the run; a stack with no significant cluster makes qdecr_snap() stop, and
# for the sex stack that is recorded rather than fatal, since its effect may well be
# empty in one hemisphere.
snapshot <- function(stack, type, required) {
  name <- paste0(hemi, ".", stack, ".", type)
  image <- tryCatch(
    qdecr_snap(out, stack = stack, type = type, plot_brain = FALSE),
    error = function(e) {
      if (required) stop("No snapshot for ", name, ": ", conditionMessage(e), call. = FALSE)
      message("No snapshot for ", name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(image)) magick::image_write(image, file.path(staged, paste0(name, ".png")), format = "png")
}
snapshot("age", "coef", required = TRUE)
snapshot("age", "t", required = TRUE)
snapshot("sexmale", "coef", required = FALSE)

# ---------- the record ----------
# What export.R needs beyond the output directory: how long the fit took, on how many
# cores, and the software, for the site's description of the run.
jsonlite::write_json(
  list(
    hemi = hemi,
    project = out$input$project2,
    seconds = round(seconds, 1),
    nCores = out$input$n_cores,
    started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    qdecr = as.character(packageVersion("QDECR")),
    r = R.version.string,
    freesurfer = readLines(file.path(Sys.getenv("FREESURFER_HOME"), "build-stamp.txt"), n = 1),
    blas = sessionInfo()$BLAS
  ),
  file.path(staged, paste0(hemi, ".run.json")),
  auto_unbox = TRUE, pretty = TRUE
)
save_text(capture.output(sessionInfo()), "sessionInfo")
message("\n", hemi, ": done in ", round(seconds), " s; ", nrow(clusters), " significant clusters.")
