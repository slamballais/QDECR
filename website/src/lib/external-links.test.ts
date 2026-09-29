import { test } from 'node:test'
import assert from 'node:assert/strict'
import { markdownToHtml } from 'satteri'
import { externalLinks } from './external-links.ts'
import { OUT_ICON } from './out-icon.ts'

const render = (markdown: string) => markdownToHtml(markdown, { mdastPlugins: [externalLinks()] }).html.trim()

test('a link that leaves the site gets the out icon, inside the link after its text', () => {
  assert.equal(
    render('See [ABIDE](https://fcon_1000.projects.nitrc.org/indi/abide/).'),
    `<p>See <a href="https://fcon_1000.projects.nitrc.org/indi/abide/">ABIDE${OUT_ICON}</a>.</p>`,
  )
})

test('plain http links leave the site too', () => {
  assert.equal(render('[PCP](http://preprocessed-connectomes-project.org/abide/)'), `<p><a href="http://preprocessed-connectomes-project.org/abide/">PCP${OUT_ICON}</a></p>`)
})

test('links within the site, to an anchor, or to qdecr.com itself get none', () => {
  assert.equal(render('[Cite](/cite)'), '<p><a href="/cite">Cite</a></p>')
  assert.equal(render('[Below](#stacks)'), '<p><a href="#stacks">Below</a></p>')
  assert.equal(render('[Home](https://qdecr.com/cite)'), '<p><a href="https://qdecr.com/cite">Home</a></p>')
})

test('the icon is decorative, and holds no style the CSP would block', () => {
  assert.match(OUT_ICON, /^<svg class="out-icon"[^>]* aria-hidden="true"/)
  assert.doesNotMatch(OUT_ICON, /style=/)
})
