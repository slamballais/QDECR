import { test } from 'node:test'
import assert from 'node:assert/strict'
import { markdownToHtml } from 'satteri'
import { figures } from './figures.ts'

/** The HTML, without the newlines between tags, which carry no meaning here. */
const render = (markdown: string) =>
  markdownToHtml(markdown, { mdastPlugins: [figures()] }).html.replace(/>\s+</g, '><').trim()

const CREDIT = 'Data: ABIDE I, preprocessed by the PCP, <a href="/colophon#the-example-data">CC BY-NC-SA 3.0</a>.'

test('an image alone in its paragraph, with a title, becomes a figure captioned by the title', () => {
  assert.equal(
    render('![A histogram](./hist.png "Thickness per vertex.")'),
    '<figure><img src="./hist.png" alt="A histogram"><figcaption>Thickness per vertex.</figcaption></figure>',
  )
})

test('a figure of the example data carries its credit and licence after the caption', () => {
  assert.equal(
    render('![A brain](../../../assets/example/lh.age.t.png "The t-statistic.")'),
    `<figure><img src="../../../assets/example/lh.age.t.png" alt="A brain"><figcaption>The t-statistic. ${CREDIT}</figcaption></figure>`,
  )
})

test('an example figure without a title still gets the credit as its caption', () => {
  assert.equal(
    render('![A brain](../../../assets/example/lh.age.t.png)'),
    `<figure><img src="../../../assets/example/lh.age.t.png" alt="A brain"><figcaption>${CREDIT}</figcaption></figure>`,
  )
})

test('a lone image with nothing to caption is a figure without a caption', () => {
  assert.equal(render('![The logo](./logo.svg)'), '<figure><img src="./logo.svg" alt="The logo"></figure>')
})

test('an image inside a sentence stays inline', () => {
  assert.equal(render('The ![logo](./logo.svg "Logo") mark.'), '<p>The <img src="./logo.svg" alt="logo" title="Logo"> mark.</p>')
})

test('an example image inside a sentence fails the build: it would lose its credit', () => {
  assert.throws(() => render('See ![a brain](../../../assets/example/lh.age.t.png) here.'), /lh\.age\.t\.png.*paragraph of its own/)
})
