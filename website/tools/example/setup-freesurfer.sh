#!/usr/bin/env bash
# A slim FreeSurfer 7.4.1 install for the example, in the user's own directory, with no
# root needed:
#
#   bash tools/example/setup-freesurfer.sh
#
# A full install is 15 to 20 GB, and the example uses a sliver of it, about 860 MB: the
# two programs QDECR calls (mris_fwhm for the smoothness of the residuals,
# mri_surfcluster for the cluster-wise correction), Freeview with the Qt and VTK it
# bundles for the snapshots (its RUNPATH looks for them in lib/qt/lib and lib/vtk), the
# fsaverage subject the maps are on (its surfaces and labels), fsaverage6 and fsaverage5
# for the site's downsampled viewer, and the precomputed simulations mri_surfcluster
# reads (average/mult-comp-cor, for fsaverage's cortex, absolute sign: the only sign
# QDECR uses). recon-all is taken too, a text script, as the record of the file names
# -qcache writes.
#
# The tarball itself is 9.5 GB and comes at a few MB/s from the FreeSurfer server, so it
# is kept in cache/ and the download resumes if it was cut off. Delete the cache once the
# install works. Safe to run again.
#
# FreeSurfer's licence is free but personal: it comes by email after registering at
# https://surfer.nmr.mgh.harvard.edu/registration.html and has to be saved by hand, as
# $FREESURFER_HOME/license.txt. This script copies one from ~/license.txt or from
# $FS_LICENSE if either exists, and says what to do if not.

set -euo pipefail
source "$(dirname "$0")/env.sh"

cache="$QDECR_EXAMPLE_ROOT/cache"
tarball="$cache/$FS_TARBALL"
mkdir -p "$cache" "$FREESURFER_HOME"

# ---------- download ----------

# The server's size for the file is the check that a download is complete; with -C -,
# curl carries on from where a partial file ends.
expected=$(curl -sI -L --max-time 60 "$FS_URL" | tr -d '\r' | awk 'tolower($1) == "content-length:" { size = $2 } END { print size }' || true)
if [ -z "$expected" ]; then
  echo "Could not read the size of $FS_URL; is the FreeSurfer server reachable?" >&2
  exit 1
fi
have=0
[ -f "$tarball" ] && have=$(stat -c %s "$tarball")
if [ "$have" -lt "$expected" ]; then
  echo "Fetching FreeSurfer $FS_VERSION ($((expected / 1024 / 1024)) MB; $((have / 1024 / 1024)) MB so far) into $cache"
  curl -L -C - --retry 10 --retry-delay 15 --retry-all-errors -o "$tarball" "$FS_URL"
fi
have=$(stat -c %s "$tarball")
if [ "$have" -ne "$expected" ]; then
  echo "$tarball is $have bytes, the server says $expected: the download is incomplete or the file changed." >&2
  exit 1
fi

# ---------- extract ----------

# Members are named ./freesurfer/...; the two leading parts are dropped so that
# $FREESURFER_HOME/bin/mris_fwhm is where FreeSurfer's own setup expects it. The list
# is kept in the install, in slim-install.txt, which also marks the extraction as done:
# when the list here differs from the recorded one, the extraction runs again, so a
# member added to this script reaches an existing install too.
members=(
  './freesurfer/SetUpFreeSurfer.sh'
  './freesurfer/FreeSurferEnv.sh'
  './freesurfer/build-stamp.txt'
  './freesurfer/FreeSurferColorLUT.txt'
  './freesurfer/bin/mris_fwhm'
  './freesurfer/bin/mri_surfcluster'
  './freesurfer/bin/freeview'
  # Tells Freeview's Qt where its plugins are (lib/qt/plugins); without it Qt cannot
  # start its X11 platform plugin.
  './freesurfer/bin/qt.conf'
  './freesurfer/bin/recon-all'
  './freesurfer/lib/qt/*'
  './freesurfer/lib/vtk/*'
  './freesurfer/subjects/fsaverage/surf/*'
  './freesurfer/subjects/fsaverage/label/*'
  # mri_surfcluster reads the subject's Talairach transform, for the MNI coordinates of
  # each cluster's peak in its summary, and exits without it.
  './freesurfer/subjects/fsaverage/mri/transforms/talairach.xfm'
  # The viewer's mesh (export.R uses fsaverage6) and the coarser fsaverage5, kept so
  # that the viewer can move to the smaller one, should its size budget call for it,
  # without another pass over the tarball.
  './freesurfer/subjects/fsaverage6/surf/*'
  './freesurfer/subjects/fsaverage6/label/*'
  './freesurfer/subjects/fsaverage5/surf/*'
  './freesurfer/subjects/fsaverage5/label/*'
  './freesurfer/average/mult-comp-cor/fsaverage/lh/cortex/*/abs/*'
  './freesurfer/average/mult-comp-cor/fsaverage/rh/cortex/*/abs/*'
)
record="$FREESURFER_HOME/slim-install.txt"
wanted=$(printf '  %s\n' "${members[@]}")
if [ ! -f "$record" ] || [ "$(tail -n +2 "$record")" != "$wanted" ]; then
  echo "Extracting the slim install into $FREESURFER_HOME (one pass over the tarball, a few minutes)"
  tar -xzf "$tarball" -C "$FREESURFER_HOME" --strip-components=2 --wildcards "${members[@]}"
  {
    echo "A slim FreeSurfer install made by tools/example/setup-freesurfer.sh on $(date -I), from $FS_TARBALL:"
    echo "$wanted"
  } > "$record"
fi
echo "FreeSurfer build: $(cat "$FREESURFER_HOME/build-stamp.txt")"
echo "Install size: $(du -sh "$FREESURFER_HOME" | cut -f1)"

# ---------- shared libraries ----------

# Anything "not found" here is a library setup-system.sh should add. The three programs
# find their bundled Qt and VTK through their RUNPATH; Qt's X11 plugin, which Freeview
# loads at run time, is checked with those on the library path.
missing=$(
  {
    ldd "$FREESURFER_HOME/bin/mris_fwhm" "$FREESURFER_HOME/bin/mri_surfcluster" "$FREESURFER_HOME/bin/freeview"
    LD_LIBRARY_PATH="$FREESURFER_HOME/lib/qt/lib:$FREESURFER_HOME/lib/vtk" ldd "$FREESURFER_HOME/lib/qt/plugins/platforms/libqxcb.so"
  } 2>&1 | grep "not found" | sort -u || true
)
if [ -n "$missing" ]; then
  echo "FreeSurfer's programs lack shared libraries (setup-system.sh installs them):" >&2
  echo "$missing" >&2
  exit 1
fi
echo "mris_fwhm, mri_surfcluster and Freeview have every library they need."

# ---------- licence ----------

licence="$FREESURFER_HOME/license.txt"
if [ ! -f "$licence" ]; then
  for candidate in "${FS_LICENSE:-}" "$HOME/license.txt"; do
    if [ -n "$candidate" ] && [ -f "$candidate" ]; then
      cp "$candidate" "$licence"
      echo "Copied the licence from $candidate"
      break
    fi
  done
fi
if [ -f "$licence" ]; then
  echo "Licence: $licence"
else
  cat <<MSG

FreeSurfer's tools will not run until its licence file is in place. It is free:
register at https://surfer.nmr.mgh.harvard.edu/registration.html, and save the
license.txt that arrives by email as

  $licence

then run tools/example/run.sh.
MSG
fi
