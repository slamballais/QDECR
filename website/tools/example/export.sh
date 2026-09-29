#!/usr/bin/env bash
# Puts the example's results into the site, and reports what it added and how big:
#
#   bash tools/example/export.sh
#
# After run.sh. The work is export.R's; this finds the repository from where the script
# is and shows the size of what will be committed, so the viewer's assets can be checked
# against their budget (under 1 MB a hemisphere) before they are. export.R also has
# Freeview draw the home page's poster, so it gets a display the way run.sh does: Xvfb,
# or the desktop with QDECR_EXAMPLE_DISPLAY=desktop.

set -euo pipefail
source "$(dirname "$0")/env.sh"

runner=()
if [ "${QDECR_EXAMPLE_DISPLAY:-xvfb}" = xvfb ]; then
  runner=(xvfb-run -a -s "-screen 0 1600x1200x24")
fi
"${runner[@]}" Rscript "$EXAMPLE_DIR/export.R" "$REPO_DIR"

website="$REPO_DIR/website"
echo
echo "Sizes:"
du -sh "$website/src/data/example" "$website/src/assets/example" "$website/public/viewer" | sed 's/^/  /'
echo
echo "Viewer files:"
ls -l "$website/public/viewer" | awk 'NR > 1 { printf "  %8d  %s\n", $5, $9 }'
