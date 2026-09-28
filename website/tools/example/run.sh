#!/usr/bin/env bash
# Runs the example analysis, both hemispheres, and makes everything the site shows of it:
#
#   bash tools/example/run.sh
#
# After setup-system.sh, setup-freesurfer.sh (with the licence in place), setup-r.sh and
# download.sh. Each hemisphere is one Rscript run of run.R, with everything it prints
# kept in results/<hemi>.log.txt for the site's output blocks. A rerun overwrites the
# results.
#
# Freeview draws the snapshots in a window. By default that window is on a virtual
# screen from Xvfb, so the run needs no desktop and the snapshots come out the same size
# every time; that is how the site's were made. QDECR_EXAMPLE_DISPLAY=desktop uses the
# real display instead: under WSLg, Freeview then opens on the Windows desktop for a few
# seconds per snapshot, which is how a reader sees it.

set -euo pipefail
source "$(dirname "$0")/env.sh"

if [ ! -f "$FREESURFER_HOME/license.txt" ]; then
  echo "No $FREESURFER_HOME/license.txt: FreeSurfer's tools refuse to run without it. See setup-freesurfer.sh." >&2
  exit 1
fi
if [ ! -e "$SUBJECTS_DIR/fsaverage" ]; then
  ln -sfn "$FREESURFER_HOME/subjects/fsaverage" "$SUBJECTS_DIR/fsaverage"
fi

results="$QDECR_EXAMPLE_ROOT/results"
mkdir -p "$results"

runner=()
if [ "${QDECR_EXAMPLE_DISPLAY:-xvfb}" = xvfb ]; then
  if ! command -v xvfb-run > /dev/null; then
    echo "xvfb-run is not installed (setup-system.sh), and QDECR_EXAMPLE_DISPLAY is not 'desktop'." >&2
    exit 1
  fi
  runner=(xvfb-run -a -s "-screen 0 1600x1200x24")
fi

for hemi in lh rh; do
  echo "==== $hemi ===="
  "${runner[@]}" Rscript "$EXAMPLE_DIR/run.R" "$hemi" 2>&1 | tee "$results/$hemi.log.txt"
done

echo
echo "Done. Results in $results; now bash tools/example/export.sh puts them in the site."
