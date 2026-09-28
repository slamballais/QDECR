import { test } from 'node:test'
import assert from 'node:assert/strict'
import { markdownToHtml } from 'satteri'
import { callouts } from './callouts.ts'

/** The HTML, without the newlines between tags, which carry no meaning here. */
const render = (markdown: string) =>
  markdownToHtml(markdown, { mdastPlugins: [callouts()] }).html.replace(/>\s+</g, '><').trim()

test('a [!NOTE] blockquote becomes a labelled note callout, its content intact', () => {
  assert.equal(
    render('> [!NOTE]\n> Set `dir_tmp` to [shared memory](/tutorials/performance).'),
    '<aside class="callout" aria-label="Note"><span class="callout-label">Note</span>' +
      '<p>Set <code>dir_tmp</code> to <a href="/tutorials/performance">shared memory</a>.</p></aside>',
  )
})

test('tips and warnings get the classes that colour them', () => {
  assert.match(render('> [!TIP]\n> Use four cores.'), /^<aside class="callout callout-tip" aria-label="Tip"><span class="callout-label">Tip<\/span>/)
  assert.match(
    render('> [!WARNING]\n> Windows needs WSL2.'),
    /^<aside class="callout callout-warning" aria-label="Warning"><span class="callout-label">Warning<\/span>/,
  )
})

test('a marker without a style fails loudly, naming the ones there are', () => {
  assert.throws(() => render('> [!CAUTION]\n> Hot.'), /\[!CAUTION\].*NOTE, TIP, WARNING/)
})

test('a plain blockquote, or one that only mentions a marker later on, stays a blockquote', () => {
  assert.equal(render('> A quote.'), '<blockquote><p>A quote.</p></blockquote>')
  assert.equal(render('> See [!NOTE] below.'), '<blockquote><p>See [!NOTE] below.</p></blockquote>')
})

test('a callout keeps every paragraph, and a marker on a line of its own leaves no empty one', () => {
  assert.equal(
    render('> [!TIP]\n>\n> One.\n>\n> Two.'),
    '<aside class="callout callout-tip" aria-label="Tip"><span class="callout-label">Tip</span><p>One.</p><p>Two.</p></aside>',
  )
})
