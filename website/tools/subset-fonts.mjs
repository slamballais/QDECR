// Cuts the site's two variable typefaces down to what its CSS uses, and writes them with
// their @font-face rules:
//
//   npm run make:fonts
//
//   src/assets/fonts/*.woff2     Newsreader and IBM Plex Sans (upright and italic),
//                                Latin and Latin Extended
//   src/styles/fonts.css         the @font-face rules for them, imported by Base.astro
//
// Fontsource ships each typeface with every weight (and, for Newsreader, every optical
// size) the design allows. The site sets Newsreader only at 600, on headings from 20 to
// about 58 pixels, and Plex Sans from 400 to 600, so this keeps those ranges of the
// variable axes and drops the rest, which is where most of a variable font's size is:
// Newsreader's Latin file goes from 132 KB to 41 KB. Every glyph and every OpenType
// feature stays (the tables' tabular figures need theirs). The fonts are most of what the
// home page weighs, and its budget counts them.
//
// Latin and Latin Extended cover English and the names on the site. Text in any other
// script falls back to the system's fonts. Plex Mono is left to Fontsource: it is static,
// and its Latin files are already 15 KB.
//
// Needs subset-font (HarfBuzz's subsetter, compiled to WebAssembly). Run it after
// updating a Fontsource package, or when the CSS starts using a weight or size outside
// the ranges below.

import subsetFont from 'subset-font'
import { mkdir, readFile, writeFile } from 'node:fs/promises'
import { createRequire } from 'node:module'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'

const require = createRequire(import.meta.url)
const WEBSITE = fileURLToPath(new URL('../', import.meta.url))
const OUT_DIR = join(WEBSITE, 'src/assets/fonts')
const CSS_OUT = join(WEBSITE, 'src/styles/fonts.css')

// The scripts kept, by Fontsource's names for its files.
const SUBSETS = ['latin-ext', 'latin']

const FACES = [
  {
    // Headings and the wordmark: SemiBold only. The optical size follows the font size
    // (font-optical-sizing: auto), from h3 (1.25rem) to the display size (3.6rem).
    css: '@fontsource-variable/newsreader/opsz.css',
    axes: { wght: 600, opsz: { min: 18, max: 60 } },
    weight: '600',
  },
  {
    // Text: regular, medium (labels) and semibold (strong, table headings).
    css: '@fontsource-variable/ibm-plex-sans/wght.css',
    axes: { wght: { min: 400, max: 600 } },
    weight: '400 600',
  },
  {
    css: '@fontsource-variable/ibm-plex-sans/wght-italic.css',
    axes: { wght: { min: 400, max: 600 } },
    weight: '400 600',
  },
]

/** The @font-face rules of a Fontsource stylesheet, with what this script needs of each. */
function fontFaces(css) {
  return [...css.matchAll(/@font-face\s*{([^}]*)}/g)].map(([, body]) => {
    const get = (property) => new RegExp(`${property}:\\s*([^;]+);`).exec(body)?.[1]?.trim()
    const src = /url\(\.\/files\/([^)]+)\)\s*(format\([^)]*\))/.exec(get('src') ?? '')
    if (!src) throw new Error(`No src in @font-face {${body}}`)
    return {
      family: get('font-family'),
      style: get('font-style'),
      display: get('font-display'),
      unicodeRange: get('unicode-range'),
      file: src[1],
      format: src[2],
    }
  })
}

await mkdir(OUT_DIR, { recursive: true })
const rules = []
for (const face of FACES) {
  const cssPath = require.resolve(face.css)
  const faces = fontFaces(await readFile(cssPath, 'utf8'))
  for (const subset of SUBSETS) {
    // newsreader-latin-opsz-normal.woff2 is Latin; newsreader-latin-ext-... is not.
    const source = faces.find((candidate) => new RegExp(`-${subset}-(?!ext-)`).test(candidate.file))
    if (!source) throw new Error(`${face.css} has no ${subset} file`)
    const input = await readFile(join(dirname(cssPath), 'files', source.file))
    const output = await subsetFont(input, null, {
      targetFormat: 'woff2',
      // Instance the axes, keep the glyphs: Fontsource has already split the scripts.
      keepAllGlyphs: true,
      variationAxes: face.axes,
      // The licence's name records (description and URL), which the subsetter drops by
      // default. The OFL asks for the licence to travel with the font.
      preserveNameIds: [13, 14],
    })
    await writeFile(join(OUT_DIR, source.file), output)
    console.log(`${source.file}: ${(input.length / 1024).toFixed(1)} KB -> ${(output.length / 1024).toFixed(1)} KB`)
    rules.push(
      `/* ${source.file.replace(/\.woff2$/, '')} */\n` +
        `@font-face {\n` +
        `  font-family: ${source.family};\n` +
        `  font-style: ${source.style};\n` +
        `  font-display: ${source.display};\n` +
        `  font-weight: ${face.weight};\n` +
        `  src: url('../assets/fonts/${source.file}') ${source.format};\n` +
        `  unicode-range: ${source.unicodeRange};\n` +
        `}\n`,
    )
  }
}

await writeFile(
  CSS_OUT,
  `/*\n` +
    ` * Written by tools/subset-fonts.mjs (npm run make:fonts); edit that, not this.\n` +
    ` * Newsreader and IBM Plex Sans, cut to the weights and sizes the site uses, from\n` +
    ` * Fontsource's files and with Fontsource's unicode ranges, so a page fetches only the\n` +
    ` * scripts its text needs.\n` +
    ` */\n\n` +
    rules.join('\n'),
)
console.log(`Wrote ${rules.length} @font-face rules to src/styles/fonts.css`)
