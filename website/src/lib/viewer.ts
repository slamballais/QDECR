// The home page viewer's settings and arithmetic, kept apart from the
// script that drives NiiVue (src/scripts/viewer-3d.ts) so they can be tested on Node:
// the files it loads, where the camera sits and how it turns between views, and the
// colours of the map.

import type { ExampleRun } from './example'

export const HEMISPHERES = ['lh', 'rh'] as const
export type Hemisphere = (typeof HEMISPHERES)[number]
export type Side = 'lateral' | 'medial'

/** What the viewer shows first: the poster, which is this hemisphere from this side. */
export const FIRST_VIEW: { hemi: Hemisphere; side: Side } = { hemi: 'lh', side: 'lateral' }

/**
 * The files tools/example/export.R writes for a hemisphere into public/viewer/, each a
 * gzipped MZ3: fsaverage6's inflated surface, 1 in its sulci and 0 on its gyri, and the
 * poster's map, the age stack's −log10(p) on its significant clusters and 0 off them.
 */
export function viewerFiles(hemi: Hemisphere) {
  return {
    mesh: `/viewer/${hemi}.inflated.mz3`,
    sulci: `/viewer/${hemi}.sulci.mz3`,
    map: `/viewer/${hemi}.age.p.mz3`,
  }
}

/**
 * NiiVue's camera azimuth, in degrees, for looking at a side of a hemisphere. At 90 the
 * camera is at the subject's left, which is the left hemisphere's lateral side and the
 * right hemisphere's medial one; at 270, the other way round.
 */
export function azimuthFor(hemi: Hemisphere, side: Side): number {
  return (hemi === 'lh') === (side === 'lateral') ? 90 : 270
}

/**
 * The azimuth that shows the other hemisphere as this one was shown: mirrored across the
 * midline, so a view from the side stays on the same side of the brain (lateral stays
 * lateral) and a view from the front or back stays put.
 */
export function mirrorAzimuth(azimuth: number): number {
  return (360 - azimuth) % 360
}

/** How long the camera takes to turn to another side, unless the reader asks for less
 * motion, when it cuts straight there. */
export const TURN_MS = 600

/** The camera's azimuth `progress` of the way (0 to 1) through a turn, eased in and out,
 * the short way round and within 0 to 360 degrees. */
export function turnAzimuth(from: number, to: number, progress: number): number {
  // The signed difference, folded into -180 to 180 so the turn is never the long way.
  const delta = ((((to - from) % 360) + 540) % 360) - 180
  // Smoothstep: slow off the mark, slow into place.
  const eased = progress * progress * (3 - 2 * progress)
  return (((from + delta * eased) % 360) + 360) % 360
}

/** The poster's colour scale, in −log10(p): where colour starts and where it saturates
 * (run.json's poster.scale). */
export type Scale = ExampleRun['poster']['scale']

/**
 * The folds, drawn as Freeview draws them, in two greys by the sign of the curvature: the
 * surface is the gyri's grey, and the sulci layer (1 in a sulcus) paints the darker one
 * over it, as a one-colour NiiVue colormap.
 */
export const GYRI: [number, number, number, number] = [150, 150, 150, 255]
export const SULCI = { R: [90, 90], G: [90, 90], B: [90, 90], A: [255, 255], I: [0, 255] }

/**
 * Freeview's heat colour scale, as the poster draws it, as a NiiVue colormap: red from the
 * low end of the scale (run.json's poster.scale.from) to its midpoint, which Freeview
 * puts halfway when given only the two ends, then red to yellow at the top (.to); and
 * opaque throughout, which is what Freeview's linearopaque overlay method means. NiiVue
 * spreads the index I (0 to 255) over the layer's cal_min to cal_max.
 */
export const HEAT = {
  R: [255, 255, 255],
  G: [0, 0, 255],
  B: [0, 0, 0],
  A: [255, 255, 255],
  I: [0, 128, 255],
}
