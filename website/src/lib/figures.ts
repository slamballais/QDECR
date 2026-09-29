// Figures in Markdown. An image in a paragraph of its own becomes a <figure>, and its
// title, the quoted text after the path, becomes the caption:
//
//   ![Histogram of the mean thickness at each vertex](../../../assets/example/lh.hist-vertex.png "Mean thickness per vertex.")
//
// A figure drawn from the example data (anything under assets/example/) also carries the
// data's credit and licence after its caption, which ABIDE's licence asks of every one, so
// a page cannot show one without it. Such an image inside a sentence would have no
// caption to carry it, so that fails the build.

import { defineMdastPlugin } from 'satteri'
import type { Image, PhrasingContent } from 'mdast'
import { EXAMPLE_CREDIT, EXAMPLE_LICENCE } from './example.ts'

/** Where tools/example/export.R writes the figures. */
const EXAMPLE_ASSET = /(^|\/)assets\/example\//

/** The credit, as ExampleCredit.astro prints it on the .astro pages. */
const credit = (): PhrasingContent[] => [
  { type: 'text', value: `${EXAMPLE_CREDIT.text}, ` },
  { type: 'link', url: EXAMPLE_CREDIT.href, children: [{ type: 'text', value: EXAMPLE_LICENCE }] },
  { type: 'text', value: '.' },
]

export function figures() {
  return defineMdastPlugin({
    name: 'qdecr:figures',
    paragraph(node, ctx) {
      const loneImage = node.children.length === 1 ? node.children[0] : undefined
      if (loneImage?.type !== 'image') {
        const inline = node.children.find((child): child is Image => child.type === 'image' && EXAMPLE_ASSET.test(child.url))
        if (inline) throw new Error(`${inline.url} is from the example data and needs its credit: put it in a paragraph of its own.`)
        return
      }
      const caption: PhrasingContent[] = []
      if (loneImage.title) caption.push({ type: 'text', value: loneImage.title })
      if (EXAMPLE_ASSET.test(loneImage.url)) {
        if (caption.length) caption.push({ type: 'text', value: ' ' })
        caption.push(...credit())
      }
      // Still a paragraph in the Markdown tree, rendered as a <figure> (hName); the caption
      // is an emphasis rendered as <figcaption>. The title moves into the caption rather
      // than stay as a tooltip that says it again.
      const image: Image = { type: 'image', url: loneImage.url, alt: loneImage.alt ?? '' }
      ctx.replaceNode(node, {
        type: 'paragraph',
        data: { hName: 'figure' },
        children: [image, ...(caption.length ? [{ type: 'emphasis' as const, data: { hName: 'figcaption' }, children: caption }] : [])],
      })
    },
  })
}
