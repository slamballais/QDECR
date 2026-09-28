import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readTokens } from './tokens.ts'

const CSS = `
:root {
  color-scheme: light dark;
  --paper: #f7f5f0;
  --ink: #17191D;
  --serif: 'Newsreader Variable', serif;
  --shadow: rgb(0 0 0 / 0.1);
}

@media screen and (prefers-color-scheme: dark) {
  :root:not([data-theme='light']) {
    --paper: #0d0f12;
    --ink: #ebe7de;
  }
}

@media screen {
  :root[data-theme='dark'] {
    --paper: #0d0f12;
  }
}
`

test('reads the hex colour tokens of a block, and nothing else', () => {
  assert.deepEqual(readTokens(CSS, ':root'), { '--paper': '#f7f5f0', '--ink': '#17191d' })
})

test('finds a block nested in a media query, and stops at its own closing brace', () => {
  assert.deepEqual(readTokens(CSS, ":root:not([data-theme='light'])"), {
    '--paper': '#0d0f12',
    '--ink': '#ebe7de',
  })
  assert.deepEqual(readTokens(CSS, ":root[data-theme='dark']"), { '--paper': '#0d0f12' })
})

test('a selector that is only the start of another does not match it', () => {
  // ':root' must not pick up ':root:not(...)' or ':root[...]'.
  const css = `:root[data-theme='dark'] {\n  --a: #000000;\n}\n:root {\n  --a: #ffffff;\n}\n`
  assert.deepEqual(readTokens(css, ':root'), { '--a': '#ffffff' })
})

test('a missing block is an error, not an empty palette', () => {
  assert.throws(() => readTokens(CSS, '.nowhere'), /no \.nowhere block/)
})
