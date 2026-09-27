// Draws the QDECR logo from its parameters and writes every file made from it:
//
//   npm run make:logo
//
//   src/assets/logo/qdecr-logo.svg          the rich master: glow, gradients, wordmark
//   src/assets/logo/qdecr-mark.svg          the flat mark for small sizes (no filters)
//   src/assets/logo/qdecr-logo-240.png      README size (240px wide, shown at ~139px high)
//   src/assets/logo/qdecr-sticker.svg       print file: vector, 1.732" x 2", letters as outlines
//   src/assets/logo/qdecr-sticker-300dpi.png  the same at 300 dpi, for shops that want pixels
//   public/favicon.svg, favicon.ico, apple-touch-icon.png
//   ../man/figures/logo.png                 the README image again, where R packages keep it
//                                           (the usethis and pkgdown convention)
//
// The logo is a point-up hexagon (the hexb.in sticker standard) around a black triangle
// pointing down, over a flat lattice of equal equilateral triangles standing for a surface
// mesh. It comes from the cover of Sander's thesis; the parameters were picked by eye over
// four rounds of drafts, and the round codes are noted next to each one so a change can be
// traced back to the comparison that settled it.
//
// It is drawn from numbers rather than in an editor so the geometry stays exact: the
// lattice is anchored on the triangle, so the triangle's corners always sit on lattice
// points, whatever the size or density.
//
// The logo is the one place the site uses gradients and glow. It is brand artwork, and
// the site's own "no ornament" rule applies to the pages around it.

import sharp from 'sharp'
import { mkdir, writeFile } from 'node:fs/promises'
import { dirname, join, relative } from 'node:path'
import { fileURLToPath } from 'node:url'
import { loadFace, typeset } from './font-outline.mjs'

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..')
const ASSETS = join(ROOT, 'src', 'assets', 'logo')
const PUBLIC = join(ROOT, 'public')
const MAN_FIGURES = join(ROOT, '..', 'man', 'figures')
const NEWSREADER = join(ROOT, 'node_modules/@fontsource-variable/newsreader/files/newsreader-latin-opsz-normal.woff2')

/* ---------- parameters ---------- */

// Colours sampled from the thesis cover, which uses FreeSurfer's overlay colours: the cool
// end (blue to teal) for the lattice, the heat end (orange to amber) for the glow. These
// seed the site palette later.
const C = {
  bg: '#020305',     // the cover's black, a hair off pure black so it prints as rich black
  blue: '#2d6cc4',   // lattice, left
  mid: '#2a8fb0',    // lattice, middle; also the flat mark's single lattice colour
  teal: '#25b096',   // lattice, right
  orange: '#ff7a1c', // triangle edge, glow, hexagon border
  amber: '#ffb14e',  // inner glow, and the lattice where the glow falls on it
  core: '#fff4cf',   // the hot line along the triangle's edge
  ink: '#f3efe6',    // the wordmark
}

const RICH = {
  n: 4,          // lattice units along one side of the triangle (round 1, D4)
  size: 0.68,    // triangle side as a share of the hexagon's flat-to-flat width (S2)
  lift: -0.1,    // raises the triangle's bounding box, in half-heights of the hexagon (L3)
  glow: 1.6,     // spread and strength of the halo (G4)
  tint: 0.9,     // how strongly the lattice picks up the orange near the triangle (G4)
  fade: 0.6,     // how far the lattice darkens away from the triangle, as on the cover (V2)
  line: 0.9,     // lattice stroke width
  border: 5,     // hexagon outline, drawn on the edge and clipped to half this (B2)
}

// The wordmark runs down the triangle's left edge, from the top-left corner towards the
// bottom corner, letters standing on the edge (round 2, E1). Set along the top edge it
// looked off.
const WORDMARK = {
  text: 'QDECR',
  weight: 600,   // Newsreader SemiBold (round 1, W3)
  tracking: 0.04,
  pad: 0.06,     // gap from the edge to the baseline, as a share of the triangle's side
  maxLength: 0.8, // longest the word may run, as a share of the side
}

// The flat mark: no filters, one lattice colour, a heavy outline so the triangle holds at
// 16px (round 4, F11: the sticker's own density with thinner lines).
const FLAT = { n: 4, size: 0.72, lift: -0.1, stroke: 12, border: 16, line: 3.5, latticeOpacity: 0.9 }

/* ---------- geometry ---------- */

// The hexagon is 200 units point to point, so one unit is 0.01" on the 2" sticker, and
// 100·√3 ≈ 173.2 flat to flat. It is centred on the origin with a point at the top.
const SQ3 = Math.sqrt(3)
const HEX_W = 100 * SQ3
const HEX = [[0, -100], [HEX_W / 2, -50], [HEX_W / 2, 50], [0, 100], [-HEX_W / 2, 50], [-HEX_W / 2, -50]]
const HEX_VIEW = `${r(-HEX_W / 2)} -100 ${r(HEX_W)} 200`

function r(v) {
  return +v.toFixed(2)
}
const points = (list) => list.map(([x, y]) => `${r(x)},${r(y)}`).join(' ')

/** The triangle's side, height and corners. `lift` moves its bounding box, not its centroid. */
function triangle({ size, lift }) {
  const S = size * HEX_W
  const h = (S * SQ3) / 2
  const top = -lift * 100 - h / 2
  return { S, h, top, corners: [[-S / 2, top], [S / 2, top], [0, top + h]] }
}

/**
 * The lattice as one path. It is anchored on the triangle's top-left corner with spacing
 * S/n, so the top edge is a horizontal lattice line and the slanted edges are ±60° lattice
 * lines: the corners land on lattice points by construction. Each line is cut to the
 * hexagon's circumcircle (radius 100) to keep the file small; the hexagon's clip path does
 * the exact crop, which is what leaves the lattice cut at an angle along its edges.
 */
function lattice(n, t) {
  const a = t.S / n
  const R = 101
  const directions = [[1, 0], [0.5, SQ3 / 2], [-0.5, SQ3 / 2]]
  let d = ''
  for (const [ux, uy] of directions) {
    // Lines in this family pass through (x0 + i·a, top); their spacing, measured across
    // them, is a·√3/2 for every family.
    const [nx, ny] = [-uy, ux]
    const step = (a * SQ3) / 2
    const origin = (-t.S / 2) * nx + t.top * ny // signed distance of the i = 0 line
    const first = Math.ceil((-R - origin) / step)
    const last = Math.floor((R - origin) / step)
    for (let i = first; i <= last; i++) {
      const dist = origin + i * step
      const half = Math.sqrt(R * R - dist * dist)
      const [cx, cy] = [dist * nx, dist * ny]
      d += `M${r(cx - half * ux)} ${r(cy - half * uy)}L${r(cx + half * ux)} ${r(cy + half * uy)}`
    }
  }
  return d
}

/**
 * The wordmark, sized the way the approved draft was. Any line parallel to a side of the
 * triangle, at distance d inside it, is S − 2d/√3 long; the letters stand on the edge, so
 * their cap line (pad + cap height in) is the one that must fit, with 8% to spare. The
 * word is also held to maxLength of the side. Optical size follows the type size, so the
 * fit is repeated until the size settles.
 */
function wordmark(face, t) {
  const pad = WORDMARK.pad * t.S
  let size = 0.15 * t.S
  let run
  for (let i = 0; i < 4; i++) {
    run = typeset(face, WORDMARK.text, { weight: WORDMARK.weight, size, tracking: WORDMARK.tracking })
    // Per unit of type size, as the browser measured it: advances plus trailing tracking.
    const length = run.width / size + WORDMARK.tracking
    const cap = run.capHeight / size
    const byRoom = (0.92 * (t.S - (2 * pad) / SQ3)) / (length + (2 * cap) / SQ3)
    size = Math.min(byRoom, (WORDMARK.maxLength * t.S) / length)
  }
  // Turn the frame 60° about the left edge's midpoint: x now runs down the edge towards
  // the bottom corner, and "up" for the letters points into the triangle.
  const [mx, my] = [-t.S / 4, t.top + t.h / 2]
  return `<path d="${run.d}" fill="${C.ink}" transform="translate(${r(mx)} ${r(my)}) rotate(60) translate(${r(-run.width / 2)} ${r(-pad)})"/>`
}

/* ---------- drawings ---------- */

/** The rich version. `size` adds physical width and height for the print file. */
function rich(face, { size } = {}) {
  const p = RICH
  const t = triangle(p)
  const tri = points(t.corners)
  const hex = points(HEX)
  const lat = lattice(p.n, t)
  const region = 'filterUnits="userSpaceOnUse" x="-160" y="-160" width="320" height="320"'
  const dims = size ? ` width="${size.width}" height="${size.height}"` : ''
  // Ids are prefixed so the SVG can be inlined next to other SVGs on a page.
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${HEX_VIEW}"${dims} role="img" aria-label="QDECR">
  <defs>
    <clipPath id="qdecr-hex"><polygon points="${hex}"/></clipPath>
    <linearGradient id="qdecr-cool" gradientUnits="userSpaceOnUse" x1="${r(-HEX_W / 2)}" y1="0" x2="${r(HEX_W / 2)}" y2="0">
      <stop offset="0" stop-color="${C.blue}"/><stop offset=".5" stop-color="${C.mid}"/><stop offset="1" stop-color="${C.teal}"/>
    </linearGradient>
    <!-- The fade: the lattice at full strength around the triangle, dimming towards the hexagon's edge. -->
    <radialGradient id="qdecr-fade-ramp" gradientUnits="userSpaceOnUse" cx="0" cy="${r(t.top + t.h / 3)}" r="105">
      <stop offset=".35" stop-color="#fff"/><stop offset="1" stop-color="#fff" stop-opacity="${r(1 - p.fade)}"/>
    </radialGradient>
    <mask id="qdecr-fade" maskUnits="userSpaceOnUse" x="-160" y="-160" width="320" height="320">
      <rect x="-160" y="-160" width="320" height="320" fill="url(#qdecr-fade-ramp)"/>
    </mask>
    <!-- The tint: a blurred band along the triangle's edges, where the lattice turns amber. -->
    <filter id="qdecr-blur-band" ${region}><feGaussianBlur stdDeviation="${r(5 + 7 * p.glow)}"/></filter>
    <mask id="qdecr-tint" maskUnits="userSpaceOnUse" x="-160" y="-160" width="320" height="320">
      <polygon points="${tri}" fill="none" stroke="#fff" stroke-width="${r(10 + 16 * p.glow)}" filter="url(#qdecr-blur-band)"/>
    </mask>
    <filter id="qdecr-blur-wide" ${region}><feGaussianBlur stdDeviation="${r(7 * p.glow)}"/></filter>
    <filter id="qdecr-blur-near" ${region}><feGaussianBlur stdDeviation="${r(1.2 + 1.6 * p.glow)}"/></filter>
  </defs>
  <g clip-path="url(#qdecr-hex)">
    <rect x="-100" y="-100" width="200" height="200" fill="${C.bg}"/>
    <path d="${lat}" fill="none" stroke="url(#qdecr-cool)" stroke-width="${p.line}" stroke-opacity=".85" mask="url(#qdecr-fade)"/>
    <path d="${lat}" fill="none" stroke="${C.amber}" stroke-width="${r(p.line * 1.15)}" opacity="${p.tint}" mask="url(#qdecr-tint)"/>
    <!-- The glow: a wide orange halo and a tighter amber one. The black fill that follows
         covers their inner halves, so the light only spills outwards. -->
    <polygon points="${tri}" fill="none" stroke="${C.orange}" stroke-width="10" opacity="${r(Math.min(1, 0.6 * p.glow))}" filter="url(#qdecr-blur-wide)"/>
    <polygon points="${tri}" fill="none" stroke="${C.amber}" stroke-width="3" opacity=".9" filter="url(#qdecr-blur-near)"/>
    <polygon points="${tri}" fill="#000"/>
    <polygon points="${tri}" fill="none" stroke="${C.orange}" stroke-width="2.4"/>
    <polygon points="${tri}" fill="none" stroke="${C.core}" stroke-width=".9"/>
    ${wordmark(face, t)}
    <polygon points="${hex}" fill="none" stroke="${C.orange}" stroke-width="${p.border}"/>
  </g>
</svg>
`
}

/**
 * The flat mark: solid colours only, no filters or text, for favicons and anywhere the
 * logo is drawn smaller than about 64px. `square` widens the view to a square with the
 * hexagon centred, which is what favicon slots expect.
 */
function flat({ square = false, background } = {}) {
  const p = FLAT
  const t = triangle(p)
  const hex = points(HEX)
  const view = square ? '-100 -100 200 200' : HEX_VIEW
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${view}" role="img" aria-label="QDECR">
  <defs><clipPath id="qdecr-mark-hex"><polygon points="${hex}"/></clipPath></defs>
  ${background ? `<rect x="-100" y="-100" width="200" height="200" fill="${background}"/>` : ''}
  <g clip-path="url(#qdecr-mark-hex)">
    <rect x="-100" y="-100" width="200" height="200" fill="${C.bg}"/>
    <path d="${lattice(p.n, t)}" fill="none" stroke="${C.mid}" stroke-width="${p.line}" opacity="${p.latticeOpacity}"/>
    <polygon points="${points(t.corners)}" fill="#000" stroke="${C.orange}" stroke-width="${p.stroke}"/>
    <polygon points="${hex}" fill="none" stroke="${C.orange}" stroke-width="${p.border}"/>
  </g>
</svg>
`
}

/* ---------- rasters ---------- */

/** Rasterise at high density and downsample, rather than scaling a small render up. */
const png = (svg, width, height = width) =>
  sharp(Buffer.from(svg), { density: 1200 }).resize(width, height).png().toBuffer()

/**
 * An ICO is a 6-byte header, one 16-byte directory entry per image, then the images. PNG
 * data inside an ICO has been valid since Windows Vista and every browser accepts it.
 * (The same writer as sanderlamballais.com's tools/make-icons.mjs.)
 */
function ico(images) {
  const header = Buffer.alloc(6)
  header.writeUInt16LE(0, 0) // reserved
  header.writeUInt16LE(1, 2) // type: icon
  header.writeUInt16LE(images.length, 4)

  const entries = []
  let offset = 6 + 16 * images.length
  for (const { size, data } of images) {
    const entry = Buffer.alloc(16)
    entry.writeUInt8(size >= 256 ? 0 : size, 0) // width
    entry.writeUInt8(size >= 256 ? 0 : size, 1) // height
    entry.writeUInt8(0, 2) // palette size
    entry.writeUInt8(0, 3) // reserved
    entry.writeUInt16LE(1, 4) // colour planes
    entry.writeUInt16LE(32, 6) // bits per pixel
    entry.writeUInt32LE(data.length, 8)
    entry.writeUInt32LE(offset, 12)
    entries.push(entry)
    offset += data.length
  }
  return Buffer.concat([header, ...entries, ...images.map((image) => image.data)])
}

/* ---------- write ---------- */

const face = await loadFace(NEWSREADER)
await mkdir(ASSETS, { recursive: true })
await mkdir(PUBLIC, { recursive: true })
await mkdir(MAN_FIGURES, { recursive: true })

const logo = rich(face)
const mark = flat()
const favicon = flat({ square: true })

// The hexagon is 1.732 : 2, so a 240px-wide README image is 277px high.
const readmeHeight = Math.round((240 * 200) / HEX_W)
// hexb.in's sticker is 2" point to point; at 300 dpi that is 520 x 600 pixels.
const sticker = rich(face, { size: { width: `${r(HEX_W / 100)}in`, height: '2in' } })

const readme = await png(logo, 240, readmeHeight)
const icoSizes = [16, 32, 48]
const files = [
  [join(ASSETS, 'qdecr-logo.svg'), logo],
  [join(ASSETS, 'qdecr-mark.svg'), mark],
  [join(ASSETS, 'qdecr-logo-240.png'), readme],
  [join(MAN_FIGURES, 'logo.png'), readme],
  [join(ASSETS, 'qdecr-sticker.svg'), sticker],
  [join(ASSETS, 'qdecr-sticker-300dpi.png'), await png(logo, 520, 600)],
  [join(PUBLIC, 'favicon.svg'), favicon],
  [join(PUBLIC, 'favicon.ico'), ico(await Promise.all(icoSizes.map(async (size) => ({ size, data: await png(favicon, size) }))))],
  // iOS paints transparency black and rounds the corners itself, so the touch icon gets
  // the logo's own black behind the hexagon.
  [join(PUBLIC, 'apple-touch-icon.png'), await png(flat({ square: true, background: C.bg }), 180)],
]
for (const [path, data] of files) await writeFile(path, data)

console.log(files.map(([path]) => relative(ROOT, path)).join('\n'))
