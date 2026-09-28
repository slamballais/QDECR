/**
 * Contrast, per theme. Run with `npm test`.
 *
 * A palette fails contrast quietly: nothing breaks, the page may even look better as a
 * colour drifts towards its background, and the readers who can no longer read it are the
 * least likely to say so. So the tokens are read out of the stylesheet and checked, both
 * themes, every run. The pairs live in palette.ts.
 */

import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { contrast } from './contrast.ts'
import { PAIRS, THEME_BLOCKS, themeTokens } from './palette.ts'
import { readTokens } from './tokens.ts'

const read = (path: string) => readFileSync(new URL(path, import.meta.url), 'utf8')
const CSS = read('../styles/tokens.css')

const THEMES = themeTokens(CSS)

for (const [theme, tokens] of Object.entries(THEMES)) {
  for (const { fg, bg, min, where } of PAIRS) {
    test(`${theme}: ${fg} on ${bg}, ${where}`, () => {
      const front = tokens[fg]
      const back = tokens[bg]
      assert.ok(front, `${theme} defines no ${fg}`)
      assert.ok(back, `${theme} defines no ${bg}`)
      const ratio = contrast(front, back)
      assert.ok(ratio >= min, `${front} on ${back} is ${ratio.toFixed(2)}:1, below ${min}:1`)
    })
  }
}

test('both themes define every colour token the other does', () => {
  // A token missing from the dark block inherits its light value into the dark page, which
  // looks like a bug in whatever painted it rather than a hole in the palette.
  assert.deepEqual(Object.keys(THEMES.dark).sort(), Object.keys(THEMES.light).sort())
})

test('the two copies of the dark tokens are the same', () => {
  // One applies when the system is dark and no choice is stored, the other when the
  // toggle stored "dark". Edit one without the other and the page changes with the route
  // a reader took to dark mode.
  assert.deepEqual(readTokens(CSS, THEME_BLOCKS.darkBySystem), THEMES.dark)
})

test('the paper colours copied into public/theme.js are the stylesheet values', () => {
  // theme.js repaints the browser's theme-color with them when the toggle overrides the
  // system, and it runs before any stylesheet can be read.
  const script = read('../../public/theme.js')
  const copy = /PAPER = \{ light: '(#[0-9a-f]{6})', dark: '(#[0-9a-f]{6})' \}/.exec(script)
  assert.ok(copy, 'no PAPER = { light: ..., dark: ... } line in theme.js')
  assert.equal(copy[1], THEMES.light['--paper'])
  assert.equal(copy[2], THEMES.dark['--paper'])
})
