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
import { EXAMPLE_CREDIT_HREF, EXAMPLE_LICENCE } from './example.ts'

/** Where tools/example/export.R writes the figures. */
const EXAMPLE_ASSET = /(^|\/)assets\/example\//

/** The credit, linking to the Colophon's section on the example data. ExampleCredit.astro
 * says the same for the .astro pages. */
const credit = (): PhrasingContent[] => [
  { type: 'text', value: 'Data: ABIDE I, ' },
  { type: 'link', url: EXAMPLE_CREDIT_HREF, children: [{ type: 'text', value: EXAMPLE_LICENCE }] },
  { type: 'text', value: '.' },
]

export function figures() {
  return defineMdastPlugin({
    name: 'qdecr:figures',
    paragraph(node, ctx) {
      const only = node.children.length === 1 ? node.children[0] : undefined
      if (only?.type !== 'image') {
        const inline = node.children.find((child): child is Image => child.type === 'image' && EXAMPLE_ASSET.test(child.url))
        if (inline) throw new Error(`${inline.url} is from the example data and needs its credit: put it in a paragraph of its own.`)
        return
      }
      const caption: PhrasingContent[] = []
      if (only.title) caption.push({ type: 'text', value: only.title })
      if (EXAMPLE_ASSET.test(only.url)) {
        if (caption.length) caption.push({ type: 'text', value: ' ' })
        caption.push(...credit())
      }
      // Still a paragraph in the Markdown tree, rendered as a <figure> (hName); the caption
      // is an emphasis rendered as <figcaption>. The title moves into the caption rather
      // than stay as a tooltip that says it again.
      const image: Image = { type: 'image', url: only.url, alt: only.alt ?? '' }
      ctx.replaceNode(node, {
        type: 'paragraph',
        data: { hName: 'figure' },
        children: [image, ...(caption.length ? [{ type: 'emphasis' as const, data: { hName: 'figcaption' }, children: caption }] : [])],
      })
    },
  })
}
