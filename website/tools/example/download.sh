#!/usr/bin/env bash
# Fetches the example's data from ABIDE: the phenotype file, and for each chosen subject
# the two files QDECR reads, lh and rh thickness smoothed at 10 mm on fsaverage, about
# 130 MB in all. Then links FreeSurfer's fsaverage into the subjects directory, as Get
# started tells a reader to do for subjects processed elsewhere.
#
#   bash tools/example/download.sh
#
# Needs R (setup-r.sh) for the subject selection and a FreeSurfer install
# (setup-freesurfer.sh) for the fsaverage link. Safe to run again: files already there
# are not fetched twice.

set -euo pipefail
source "$(dirname "$0")/env.sh"

data="$QDECR_EXAMPLE_ROOT/data"
mkdir -p "$data" "$SUBJECTS_DIR"

abide_pheno="$data/Phenotypic_V1_0b_preprocessed1.csv"
if [ ! -s "$abide_pheno" ]; then
  echo "Fetching the ABIDE phenotype file"
  curl -sS --fail --retry 5 -o "$abide_pheno" "$ABIDE_PHENOTYPES"
fi

Rscript "$EXAMPLE_DIR/select-subjects.R" "$abide_pheno" "$data/phenotypes.csv" "$data/excluded.csv"

# One line per file to fetch: its URL and where it goes. Then eight downloads at a time;
# the bucket is fast, and each file is under a megabyte.
list="$data/files.txt"
tail -n +2 "$data/phenotypes.csv" | cut -d, -f1 | while read -r id; do
  for hemi in lh rh; do
    file="$hemi.thickness.fwhm10.fsaverage.mgh"
    printf '%s %s\n' "$ABIDE_FREESURFER/$id/surf/$file" "$SUBJECTS_DIR/$id/surf/$file"
  done
done > "$list"

echo "Fetching $(wc -l < "$list") thickness files into $SUBJECTS_DIR"
xargs -P 8 -L 1 bash -c '
  url=$0; dest=$1
  [ -s "$dest" ] && exit 0
  mkdir -p "$(dirname "$dest")"
  curl -sS --fail --retry 5 --retry-delay 5 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
' < "$list"

# QDECR looks for the target in the subjects directory: fsaverage, from FreeSurfer.
if [ -d "$FREESURFER_HOME/subjects/fsaverage" ]; then
  ln -sfn "$FREESURFER_HOME/subjects/fsaverage" "$SUBJECTS_DIR/fsaverage"
else
  echo "No $FREESURFER_HOME/subjects/fsaverage yet: run setup-freesurfer.sh, then this script again for the link." >&2
fi

echo "Done: $(find "$SUBJECTS_DIR" -name '*.fsaverage.mgh' | wc -l) files, $(du -sh "$SUBJECTS_DIR" | cut -f1)"
