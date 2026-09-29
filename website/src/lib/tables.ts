// Tables in Markdown, wrapped the way the /design page wraps its own: a wide table
// scrolls inside <div class="table-wrap"> instead of pushing the page sideways at 320px.
// The wrapper can take focus so a keyboard can scroll it, and so it is a named region:
// named after the heading the table sits under, or its column names if there is none. A
// second table under the same heading is "<heading>, table 2", as two regions with one
// name cannot be told apart in a screen reader's list of them.

import { defineMdastPlugin } from 'satteri'

/** A table this plugin has already wrapped. */
const isWrap = (node: { data?: unknown }) =>
  (node.data as { hProperties?: { className?: string[] } } | undefined)?.hProperties?.className?.includes('table-wrap') ?? false

export function tables() {
  return defineMdastPlugin({
    name: 'qdecr:tables',
    table(node, ctx) {
      const parent = ctx.parent(node)
      const index = ctx.indexOf(node) ?? 0
      const before = parent?.children.slice(0, index) ?? []
      const headingAt = before.findLastIndex((child) => child.type === 'heading')
      const heading = before[headingAt]
      // The tables between that heading and this one, wrapped already or not yet.
      const earlier = before.slice(headingAt + 1).filter((child) => child.type === 'table' || isWrap(child)).length
      const columns = node.children[0]?.children.map((cell) => ctx.textContent(cell)).join(', ')
      const name = heading ? ctx.textContent(heading) : `Table: ${columns ?? ''}`
      const label = earlier ? `${name}, table ${earlier + 1}` : name
      // Markdown has no plain container node, so the wrapper is a blockquote that renders
      // as a <div>: hName and hProperties set the element the HTML gets.
      ctx.wrapNode(node, {
        type: 'blockquote',
        data: {
          hName: 'div',
          hProperties: { className: ['table-wrap'], tabIndex: 0, role: 'region', ariaLabel: label },
        },
        children: [],
      })
    },
  })
}
