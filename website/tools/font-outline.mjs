// Sets text as SVG outlines from the Fontsource font files, so the logo and any other
// generated artwork carry their letters as paths: librsvg (inside sharp) cannot see web
// fonts, and a print shop should not need the font either.
//
// The WOFF2 repacking below is the same as in sanderlamballais.com's tools/make-og.mjs,
// where it was first written; see sfnt() for why it is needed.

import { readFile } from 'node:fs/promises'
import { create } from 'fontkitten'

const ON_CURVE = 1
const ARG_WORDS = 1
const ARGS_ARE_XY = 2
const MORE_COMPONENTS = 32
const WE_HAVE_TWO_BY_TWO = 128

const pad4 = (buf) => (buf.length % 4 ? Buffer.concat([buf, Buffer.alloc(4 - (buf.length % 4))]) : buf)

/** A TrueType glyph: contours as absolute points, written back out uncompressed. */
function encodeSimpleGlyph(glyph) {
  const points = glyph.points ?? []
  if (points.length === 0) return Buffer.alloc(0) // a space, or an empty glyph

  const ends = points.flatMap((point, i) => (point.endContour ? [i] : []))
  if (ends.length === 0) ends.push(points.length - 1)
  const xs = points.map((point) => Math.round(point.x))
  const ys = points.map((point) => Math.round(point.y))

  const header = Buffer.alloc(12 + ends.length * 2)
  header.writeInt16BE(ends.length, 0)
  header.writeInt16BE(Math.min(...xs), 2)
  header.writeInt16BE(Math.min(...ys), 4)
  header.writeInt16BE(Math.max(...xs), 6)
  header.writeInt16BE(Math.max(...ys), 8)
  ends.forEach((end, i) => header.writeUInt16BE(end, 10 + i * 2))
  header.writeUInt16BE(0, 10 + ends.length * 2) // no hinting instructions

  // Long-form flags and 16-bit deltas throughout: larger than a compressed glyf table and
  // simpler to be sure of, and the file never leaves this process.
  const flags = Buffer.from(points.map((point) => (point.onCurve ? ON_CURVE : 0)))
  const coords = Buffer.alloc(points.length * 4)
  let at = 0
  let previous = 0
  for (const x of xs) {
    coords.writeInt16BE(x - previous, at)
    at += 2
    previous = x
  }
  previous = 0
  for (const y of ys) {
    coords.writeInt16BE(y - previous, at)
    at += 2
    previous = y
  }
  return Buffer.concat([header, flags, coords])
}

/** A composite glyph (an accented letter, say) as component references. */
function encodeCompositeGlyph(glyph) {
  const f2dot14 = (value) => Math.max(-32768, Math.min(32767, Math.round(value * 16384)))
  const parts = glyph.components.map((component, i) => {
    const scaled =
      component.scaleX !== 1 || component.scaleY !== 1 || component.scale01 !== 0 || component.scale10 !== 0
    let flags = ARG_WORDS | ARGS_ARE_XY
    if (i < glyph.components.length - 1) flags |= MORE_COMPONENTS
    if (scaled) flags |= WE_HAVE_TWO_BY_TWO

    const part = Buffer.alloc(scaled ? 16 : 8)
    part.writeUInt16BE(flags, 0)
    part.writeUInt16BE(component.glyphID, 2)
    part.writeInt16BE(Math.round(component.dx), 4)
    part.writeInt16BE(Math.round(component.dy), 6)
    if (scaled) {
      part.writeInt16BE(f2dot14(component.scaleX), 8)
      part.writeInt16BE(f2dot14(component.scale01), 10)
      part.writeInt16BE(f2dot14(component.scale10), 12)
      part.writeInt16BE(f2dot14(component.scaleY), 14)
    }
    return part
  })
  const header = Buffer.alloc(10)
  header.writeInt16BE(-1, 0) // numberOfContours < 0 marks a composite
  return Buffer.concat([header, ...parts])
}

/**
 * Repacks a WOFF2 file as a plain TrueType font.
 *
 * fontkitten reads WOFF2 happily, but it cannot apply variation axes to one: WOFF2 stores
 * glyphs in its own transformed form, and fontkitten's reader for that form never runs the
 * gvar deltas over the points. The logo needs a variation (Newsreader at weight 600 and a
 * set optical size), so the glyphs are decoded once here, written back as an ordinary
 * glyf/loca pair, and handed to fontkitten again as a TTF, where getVariation() works.
 * Every other table is copied across untouched: WOFF2 only transforms glyf and loca.
 */
function sfnt(woff2) {
  const font = create(woff2)
  if (font.type !== 'WOFF2') return woff2
  font._transformGlyfTable() // also brotli-decompresses the table data
  const source = font.stream.buffer

  const tables = new Map()
  for (const [tag, entry] of Object.entries(font.directory.tables)) {
    if (tag === 'glyf' || tag === 'loca') continue
    const length = entry.transformLength ?? entry.length
    tables.set(tag, Buffer.from(source.buffer, source.byteOffset + entry.offset, length))
  }

  const glyphs = font._transformedGlyphs.map((glyph) =>
    pad4(glyph.numberOfContours < 0 ? encodeCompositeGlyph(glyph) : encodeSimpleGlyph(glyph)),
  )
  const loca = Buffer.alloc((glyphs.length + 1) * 4)
  let offset = 0
  glyphs.forEach((glyph, i) => {
    loca.writeUInt32BE(offset, i * 4)
    offset += glyph.length
  })
  loca.writeUInt32BE(offset, glyphs.length * 4)
  tables.set('glyf', Buffer.concat(glyphs))
  tables.set('loca', loca)

  const head = Buffer.from(tables.get('head'))
  head.writeUInt32BE(0, 8) // checkSumAdjustment, meaningless for a font kept in memory
  head.writeInt16BE(1, 50) // indexToLocFormat: the long offsets written above
  tables.set('head', head)

  const tags = [...tables.keys()].sort()
  const entrySelector = Math.floor(Math.log2(tags.length))
  const searchRange = 2 ** entrySelector * 16
  const directory = Buffer.alloc(12 + tags.length * 16)
  directory.writeUInt32BE(0x00010000, 0) // sfntVersion: TrueType outlines
  directory.writeUInt16BE(tags.length, 4)
  directory.writeUInt16BE(searchRange, 6)
  directory.writeUInt16BE(entrySelector, 8)
  directory.writeUInt16BE(tags.length * 16 - searchRange, 10)

  const checksum = (buf) => {
    let sum = 0
    for (let i = 0; i < buf.length; i += 4) sum = (sum + buf.readUInt32BE(i)) >>> 0
    return sum
  }
  const body = []
  let position = directory.length
  tags.forEach((tag, i) => {
    const data = pad4(tables.get(tag))
    const record = 12 + i * 16
    directory.write(tag.padEnd(4, ' '), record, 4, 'ascii')
    directory.writeUInt32BE(checksum(data), record + 4)
    directory.writeUInt32BE(position, record + 8)
    directory.writeUInt32BE(tables.get(tag).length, record + 12)
    body.push(data)
    position += data.length
  })
  return Buffer.concat([directory, ...body])
}

/** Loads a Fontsource WOFF2 file as a variable font that getVariation() can instance. */
export async function loadFace(path) {
  return create(sfnt(await readFile(path)))
}

/**
 * Sets a string as a single SVG path, drawn from an origin on the baseline at its left
 * end, and reports its metrics.
 *
 * Optical size follows the type size, which is what `font-optical-sizing: auto` does in a
 * browser, so the outlines match what the same text looks like on a page.
 *
 * No shaping: fontkitten has no layout engine, so this is one glyph after another at its
 * own advance width, with `tracking` (in em, like CSS letter-spacing) added between them.
 */
export function typeset(face, text, { weight, size, tracking = 0 }) {
  const axes = { wght: weight }
  const opsz = face.variationAxes.opsz
  if (opsz) axes.opsz = Math.max(opsz.min, Math.min(opsz.max, size))
  const font = face.getVariation(axes)
  const scale = size / font.unitsPerEm
  const extra = tracking * size
  let x = 0
  let d = ''
  for (const glyph of font.glyphsForString(text)) {
    // The y axis points up in a font and down in SVG, hence the negative vertical scale.
    if (glyph.path.commands.length > 0) d += glyph.path.transform(scale, 0, 0, -scale, x, 0).toSVG()
    x += glyph.advanceWidth * scale + extra
  }
  return {
    d,
    width: Math.max(0, x - extra), // first glyph's origin to last glyph's advance
    capHeight: (font.capHeight / font.unitsPerEm) * size,
  }
}
