#!/usr/bin/env bash
# Installs the R packages the example needs, into the user's own library, the way the
# Get started page tells a reader to: pak first, then QDECR from GitHub through it. Two
# more for the scripts here: magick, which qdecr_snap() needs to compose its snapshots,
# and jsonlite, which export.R uses to write the results for the site.
#
#   bash tools/example/setup-r.sh
#
# Runs after setup-system.sh (R itself, the compilers, ImageMagick's headers). Takes some
# minutes: QDECR's dependencies compile C++ code, and pak from CRAN is built from source.
# Safe to run again; pak skips what is already installed.

set -euo pipefail
source "$(dirname "$0")/env.sh"

# R only creates the user library when asked interactively; Rscript would fall back to
# the system library and fail without root.
mkdir -p "$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"

Rscript - <<'R'
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak(c("slamballais/QDECR", "magick", "jsonlite"))
R

echo
echo "Installed:"
Rscript -e 'for (p in c("QDECR", "bigstatsr", "RcppEigen", "magick", "jsonlite")) cat(sprintf("  %-10s %s\n", p, as.character(packageVersion(p))))'
Rscript -e 'cat("  BLAS:", sessionInfo()$BLAS, "\n")'
