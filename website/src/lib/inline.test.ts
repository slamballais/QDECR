import { test } from 'node:test'
import assert from 'node:assert/strict'
import { findInline } from './inline.ts'

test('a page with only external scripts and stylesheets passes', () => {
  const html = `<!doctype html><html><head>
    <link rel="stylesheet" href="/_astro/index.css">
    <script src="/theme.js"></script>
    <script type="module" src="/_astro/page.js"></script>
  </head><body><p class="lead">Hi</p></body></html>`
  assert.deepEqual(findInline(html), [])
})

test('an inline script is refused', () => {
  const found = findInline('<body><script>alert(1)</script></body>')
  assert.equal(found.length, 1)
  assert.match(found[0]!, /inline <script>/)
})

test('an inline module script is refused too', () => {
  assert.equal(findInline('<script type="module">import "/x.js"</script>').length, 1)
})

test('JSON-LD is data, not a script, and passes', () => {
  assert.deepEqual(findInline('<script type="application/ld+json">{"@type":"Thing"}</script>'), [])
})

test('a <style> element is refused', () => {
  assert.match(findInline('<style>p{color:red}</style>')[0]!, /<style>/)
})

test('a style attribute is refused, whatever its quoting', () => {
  assert.equal(findInline('<pre style="background:#000">x</pre>').length, 1)
  assert.equal(findInline("<span style='color:red'>x</span>").length, 1)
  assert.equal(findInline('<span\n  style=color:red>x</span>').length, 1)
})

test('an event handler attribute is refused', () => {
  assert.match(findInline('<button onclick="go()">Go</button>')[0]!, /onclick/)
})

test('the words appearing in text or other attributes are not mistaken for markup', () => {
  const html = `<p>Set style="x" in <code>&lt;script&gt;</code>, or the online docs.</p>
    <a href="/help#style" data-onclick-note="no">Style guide</a>`
  assert.deepEqual(findInline(html), [])
})
