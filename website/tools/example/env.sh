# Shared settings for the example scripts in this directory, sourced by each of them:
#
#   source "$(dirname "$0")/env.sh"
#
# Everything the example needs sits in one directory inside the Linux file system, the
# FreeSurfer install, the subjects and the results alike: WSL reads Windows drives
# (/mnt/c) far more slowly than its own disk. Set QDECR_EXAMPLE_ROOT to put it elsewhere.
#
# Only the scripts source this; the site's Get started page tells readers to put the same
# two exports in their ~/.bashrc.

QDECR_EXAMPLE_ROOT="${QDECR_EXAMPLE_ROOT:-$HOME/qdecr-example}"
export QDECR_EXAMPLE_ROOT

# The slim FreeSurfer install setup-freesurfer.sh makes, and the subjects directory
# download.sh fills. QDECR reads both variables, as FreeSurfer's own tools do.
export FREESURFER_HOME="$QDECR_EXAMPLE_ROOT/freesurfer"
export SUBJECTS_DIR="$QDECR_EXAMPLE_ROOT/subjects"

# FreeSurfer's own setup, once it is installed: it adds bin/ to the PATH and sets the
# variables its tools expect. It keeps a SUBJECTS_DIR that is already set. Quiet, because
# it prints its banner to every script otherwise.
if [ -f "$FREESURFER_HOME/SetUpFreeSurfer.sh" ]; then
  # shellcheck disable=SC1091
  source "$FREESURFER_HOME/SetUpFreeSurfer.sh" > /dev/null
fi

# FreeSurfer 7.4.1, the release the example was run with, and where it comes from.
export FS_VERSION="7.4.1"
export FS_TARBALL="freesurfer-linux-ubuntu22_amd64-${FS_VERSION}.tar.gz"
export FS_URL="https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FS_VERSION}/${FS_TARBALL}"

# ABIDE I as the Preprocessed Connectomes Project shares it: a public S3
# bucket with the phenotype file and, per subject, the output of FreeSurfer 5.1 with
# -qcache. No account and no key are needed to read it.
export ABIDE_BASE="https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative"
export ABIDE_PHENOTYPES="$ABIDE_BASE/Phenotypic_V1_0b_preprocessed1.csv"
export ABIDE_FREESURFER="$ABIDE_BASE/Outputs/freesurfer/5.1"

# Where the scripts are, and the repository they are in, for export.sh.
EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export EXAMPLE_DIR
REPO_DIR="$(cd "$EXAMPLE_DIR/../../.." && pwd)"
export REPO_DIR
