#!/usr/bin/env bash
# The root part of the setup: the Ubuntu packages the example needs, on a fresh Ubuntu
# inside WSL2 or any Debian-based Linux.
#
#   sudo bash tools/example/setup-system.sh
#
# What and why:
# - r-base, r-base-dev: R and the compilers QDECR's dependencies build with (RcppEigen,
#   bigstatsr).
# - libcurl4-openssl-dev, libssl-dev, pkg-config: pak from CRAN is built from source
#   and embeds the curl package, which stops without libcurl's headers; magick needs
#   curl as well.
# - libmagick++-dev: for the magick package, which qdecr_snap() composes its
#   snapshots with.
# - xvfb: a virtual display, so Freeview can take snapshots in a shell without one.
# - The rest are shared libraries FreeSurfer 7.4.1 is built against and Ubuntu's WSL
#   image does not ship: OpenMP and OpenGL for its binaries, and X11, xcb, fontconfig and
#   dbus for the Qt that Freeview bundles. Found with ldd on the slim install
#   (setup-freesurfer.sh), which prints "not found" for anything still missing.
#
# Ubuntu renames libraries between releases (the t64 suffix since 24.04), so each name
# is checked against the archive first and only the ones that exist are installed; the
# script says which it skipped. Written for 26.04, the release the example ran on.

set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

if [ "$(id -u)" -ne 0 ]; then
  echo "Run this with sudo: it installs packages." >&2
  exit 1
fi

apt-get update -q

candidates="
  r-base r-base-dev
  libcurl4-openssl-dev libssl-dev pkg-config
  libmagick++-dev
  xvfb
  libgomp1 libglu1-mesa libgl1 libegl1 libopengl0
  libxmu6 libxt6 libxext6 libxrender1 libxi6 libxss1 libsm6 libice6
  libxrandr2 libxinerama1 libxcursor1 libxfixes3 libxcomposite1 libxdamage1 libxtst6
  libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0
  libxcb-shape0 libxcb-xinerama0 libxcb-xkb1 libxcb-xfixes0 libxcb-cursor0 libxcb-util1
  libxkbcommon-x11-0 libxkbcommon0 libfontconfig1 libfreetype6 libdbus-1-3
  libglib2.0-0t64 libasound2t64 libpulse0 libnss3 libnspr4
  libjpeg-turbo8 libtiff6 libpng16-16t64 libgfortran5
"
available=""
skipped=""
for p in $candidates; do
  if apt-cache show "$p" > /dev/null 2>&1; then
    available="$available $p"
  else
    skipped="$skipped $p"
  fi
done
[ -n "$skipped" ] && echo "Not in this Ubuntu's archive, skipped:$skipped"

# shellcheck disable=SC2086
apt-get install -y -q $available

echo
echo "Installed: $(R --version | head -1)"
