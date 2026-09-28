// Callouts in Markdown, written the way GitHub writes its alerts, so the source reads the
// same on GitHub as on the site:
//
//   > [!TIP]
//   > Set `n_cores` to the number of physical cores.
//
// The blockquote becomes the <aside class="callout"> that global.css styles, with the kind
// in a label above the text.

import { defineMdastPlugin } from 'satteri'
import type { Paragraph, PhrasingContent } from 'mdast'

/** The marker on the first line of the blockquote: [!NOTE]. */
const MARKER = /^\[!([A-Za-z]+)\][ \t]*(?:\n|$)/

/** The three kinds global.css styles, by marker: the label, and the classes. */
const KINDS: Record<string, { label: string; className: string[] }> = {
  NOTE: { label: 'Note', className: ['callout'] },
  TIP: { label: 'Tip', className: ['callout', 'callout-tip'] },
  WARNING: { label: 'Warning', className: ['callout', 'callout-warning'] },
}

export function callouts() {
  return defineMdastPlugin({
    name: 'qdecr:callouts',
    blockquote(node, ctx) {
      const first = node.children[0]
      if (first?.type !== 'paragraph') return
      const lead = first.children[0]
      if (lead?.type !== 'text') return
      const marker = MARKER.exec(lead.value)
      if (!marker) return
      // GitHub also knows [!IMPORTANT] and [!CAUTION]; the site has no style for them, and
      // a typo would otherwise print as text, so either fails the build.
      const kind = KINDS[marker[1]!.toUpperCase()]
      if (!kind) {
        throw new Error(`Unknown callout [!${marker[1]}]. Use one of ${Object.keys(KINDS).join(', ')}.`)
      }
      const { label, className } = kind

      // The first paragraph, less the marker line. Nothing is left of it when the marker
      // stood in a paragraph of its own.
      const rest = lead.value.slice(marker[0].length)
      const children: PhrasingContent[] = rest ? [{ ...lead, value: rest }, ...first.children.slice(1)] : first.children.slice(1)
      const paragraph: Paragraph = { type: 'paragraph', children }

      ctx.replaceNode(node, {
        type: 'blockquote',
        data: { hName: 'aside', hProperties: { className, ariaLabel: label } },
        children: [
          { type: 'paragraph', data: { hName: 'span', hProperties: { className: ['callout-label'] } }, children: [{ type: 'text', value: label }] },
          ...(children.length ? [paragraph] : []),
          ...node.children.slice(1),
        ],
      })
    },
  })
}
