// Draws the share image for a page: the card LinkedIn, Slack, Mastodon and the rest show
// when someone posts a link to qdecr.com. 1200 x 630, the size Open Graph and Twitter's
// large card both use. The logo on the left, and on the right the section, the page's
// title and the address:
//
//   node tools/og.mjs "Performance and memory" --section "Tutorial 7" --out og.png
//   node tools/og.mjs qdecr_fastlm --section Reference --code --out og.png
//
// The build draws one per page (src/integrations/og-images.ts); the command above is for
// looking at a card while changing this file.
//
// It is the dark theme's paper and ink, read from tokens.css, so a card looks like the
// site does at night. Text is set as outlines (tools/font-outline.mjs), because librsvg,
// inside sharp, cannot see the web fonts.

import sharp from 'sharp'
import { readFile, writeFile } from 'node:fs/promises'
import { join } from 'node:path'
import { parseArgs } from 'node:util'
import { loadFace, typeset } from './font-outline.mjs'
import { themeTokens } from '../src/lib/palette.ts'

const ROOT = join(import.meta.dirname, '..')
const FONTS = join(ROOT, 'node_modules')

export const WIDTH = 1200
export const HEIGHT = 630

// The logo, and the column the text sits in to its right.
const LOGO = { x: 84, height: 450 }
const TEXT = { x: 530, right: 1116 }

const TITLE = {
  serif: { size: 76, min: 48, weight: 600, leading: 1.12 },
  code: { size: 64, min: 40, weight: 600, leading: 1.2 },
  maxLines: 3,
}

const escape = (text) => text.replace(/&/g, '&amp;').replace(/</g, '&lt;')

/**
 * Breaks a title into lines no wider than `width`, one word at a time. Words are split at
 * spaces, and a code title (a function name has none) also after its underscores and
 * dots, so qdecr_prep_mgh can break rather than overflow.
 */
function wrap(text, measure, width, code) {
  const words = code ? text.split(/(?<=[_.])/) : text.split(' ')
  const join = code ? '' : ' '
  const lines = []
  let line = ''
  for (const word of words) {
    const candidate = line ? line + join + word : word
    if (line && measure(candidate) > width) {
      lines.push(line)
      line = word
    } else {
      line = candidate
    }
  }
  if (line) lines.push(line)
  return lines
}

/** Loads the faces, the logo and the colours once, and returns a function that draws. */
export async function createCardRenderer() {
  const [serif, sans, mono, logoSvg, tokensCss] = await Promise.all([
    loadFace(join(FONTS, '@fontsource-variable/newsreader/files/newsreader-latin-opsz-normal.woff2')),
    loadFace(join(FONTS, '@fontsource-variable/ibm-plex-sans/files/ibm-plex-sans-latin-wght-normal.woff2')),
    loadFace(join(FONTS, '@fontsource/ibm-plex-mono/files/ibm-plex-mono-latin-600-normal.woff2')),
    readFile(join(ROOT, 'src/assets/logo/qdecr-logo.svg'), 'utf8'),
    readFile(join(ROOT, 'src/styles/tokens.css'), 'utf8'),
  ])
  const colour = themeTokens(tokensCss).dark

  // The logo's own <svg>, placed and sized; its viewBox keeps its proportions.
  const logoView = /viewBox="([^"]+)"/.exec(logoSvg)[1].split(/\s+/).map(Number)
  const logoWidth = (LOGO.height * logoView[2]) / logoView[3]
  const logo = logoSvg
    .replace(/<svg[^>]*>/, (open) =>
      open.replace('<svg', `<svg x="${LOGO.x}" y="${(HEIGHT - LOGO.height) / 2}" width="${logoWidth}" height="${LOGO.height}"`),
    )
    .replace(/\s*role="img"|\s*aria-label="[^"]*"/g, '')

  /** A card for one page. `section` is the small label above the title. */
  return async function render({ title, section, code = false }) {
    const face = code ? mono : serif
    const spec = code ? TITLE.code : TITLE.serif
    const width = TEXT.right - TEXT.x

    // The largest size, stepping down, at which the title fits in three lines. A function
    // name is read as one word, so it first shrinks to fit on one line, and only breaks
    // if even the smallest size will not hold it.
    const fit = (maxLines) => {
      for (let size = spec.size; ; size -= 4) {
        const measure = (text) => typeset(face, text, { weight: spec.weight, size }).width
        const lines = wrap(title, measure, width, code)
        const fits = lines.length <= maxLines && lines.every((line) => measure(line) <= width)
        if (fits || size - 4 < spec.min) return { size, lines, fits }
      }
    }
    const oneLine = code ? fit(1) : { fits: false }
    const { size, lines } = oneLine.fits ? oneLine : fit(TITLE.maxLines)

    const label = section ? typeset(sans, section.toUpperCase(), { weight: 600, size: 26, tracking: 0.08 }) : null
    const address = typeset(sans, 'qdecr.com', { weight: 500, size: 30 })
    const lineHeight = size * spec.leading

    // The block is centred on the logo: label, gap, title lines, gap, address.
    const labelBlock = label ? 26 + 30 : 0
    const titleBlock = lineHeight * (lines.length - 1) + size * 0.72
    const block = labelBlock + titleBlock + 56 + 30
    let y = (HEIGHT - block) / 2

    const parts = []
    if (label) {
      y += 26
      parts.push(`<path d="${label.d}" fill="${colour['--muted']}" transform="translate(${TEXT.x} ${y})"/>`)
      y += 30
    }
    y += size * 0.72
    for (const line of lines) {
      const run = typeset(face, line, { weight: spec.weight, size })
      parts.push(`<path d="${run.d}" fill="${colour['--ink']}" transform="translate(${TEXT.x} ${y})"/>`)
      y += lineHeight
    }
    y += 56 + 30 - lineHeight
    parts.push(`<path d="${address.d}" fill="${colour['--accent']}" transform="translate(${TEXT.x} ${y})"/>`)

    const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${WIDTH}" height="${HEIGHT}" viewBox="0 0 ${WIDTH} ${HEIGHT}">
  <title>${escape(title)}</title>
  <rect width="${WIDTH}" height="${HEIGHT}" fill="${colour['--paper']}"/>
  <rect y="${HEIGHT - 8}" width="${WIDTH}" height="8" fill="${colour['--brand']}"/>
  ${logo}
  ${parts.join('\n  ')}
</svg>`
    return sharp(Buffer.from(svg), { density: 144 }).resize(WIDTH, HEIGHT).png({ compressionLevel: 9 }).toBuffer()
  }
}

// Run directly: draw one card to a file.
if (import.meta.filename === process.argv[1]) {
  const { values, positionals } = parseArgs({
    allowPositionals: true,
    options: { section: { type: 'string' }, code: { type: 'boolean' }, out: { type: 'string', default: 'og.png' } },
  })
  const render = await createCardRenderer()
  await writeFile(values.out, await render({ title: positionals.join(' '), section: values.section, code: values.code }))
  console.log(values.out)
}
