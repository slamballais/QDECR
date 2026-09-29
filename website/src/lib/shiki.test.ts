import { test } from 'node:test'
import assert from 'node:assert/strict'
import { codeToHtml } from 'shiki'
import { codeBlock, syntaxTheme } from './shiki.ts'

const highlight = (code: string, lang: string) =>
  codeToHtml(code, { lang, theme: syntaxTheme, transformers: [codeBlock()] })

/** The text of every span carrying `className`, in order. */
const tokens = (html: string, className: string) =>
  [...html.matchAll(new RegExp(`<span class="${className}">([^<]*)</span>`, 'g'))].map((m) => m[1])

const R = 'vw <- qdecr_fastlm(y ~ age, hemi = "lh", n_cores = 4, save = TRUE) # fit'

test('R tokens get the tok-* classes global.css colours', async () => {
  const html = await highlight(R, 'r')
  assert.deepEqual(tokens(html, 'tok-keyword'), ['&#x3C;-', '~', '=', '=', '='])
  assert.deepEqual(tokens(html, 'tok-function'), ['qdecr_fastlm'])
  assert.deepEqual(tokens(html, 'tok-string'), ['"lh"'])
  assert.deepEqual(tokens(html, 'tok-number'), ['4', 'TRUE'])
  assert.deepEqual(tokens(html, 'tok-comment'), ['# fit'])
})

test('no colour is written into the markup, so the CSP can forbid inline styles', async () => {
  assert.doesNotMatch(await highlight(R, 'r'), /style=/)
})

test('commas are punctuation, not operators, and stay in the text colour', async () => {
  const html = await highlight('c(1, 2)', 'r')
  assert.doesNotMatch(html, /<span class="tok-keyword">,<\/span>/)
})

test('R output lines (#>) are marked as output', async () => {
  const html = await highlight('mean(1:4)\n#> [1] 2.5', 'r')
  assert.match(html, /<span class="line output">.*#> \[1\] 2.5/)
  assert.doesNotMatch(html.split('\n')[0]!, /output/)
})

test('the block is wrapped with its language label and a copy button', async () => {
  const html = await highlight('x <- 1', 'r')
  assert.match(html, /^<div class="code-block">/)
  assert.match(html, /<span class="code-lang">R<\/span>/)
  assert.match(html, /<button type="button" class="code-copy" hidden>Copy<\/button>/)
  // Long lines scroll inside the block, so it must be reachable by keyboard.
  assert.match(html, /<pre tabindex="0"><code>/)
})

test('shell blocks are labelled Shell, and plain text gets no label', async () => {
  assert.match(await highlight('echo $SUBJECTS_DIR', 'bash'), /<span class="code-lang">Shell<\/span>/)
  assert.doesNotMatch(await highlight('just text', 'text'), /code-lang/)
})
