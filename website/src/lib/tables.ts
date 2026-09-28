// Tables in Markdown, wrapped the way the /design page wraps its own: a wide table
// scrolls inside <div class="table-wrap"> instead of pushing the page sideways at 320px.
// The wrapper can take focus so a keyboard can scroll it, and so it is a named region:
// named after the heading the table sits under, or its column names if there is none.

import { defineMdastPlugin } from 'satteri'

export function tables() {
  return defineMdastPlugin({
    name: 'qdecr:tables',
    table(node, ctx) {
      const parent = ctx.parent(node)
      const index = ctx.indexOf(node) ?? 0
      const heading = parent?.children
        .slice(0, index)
        .reverse()
        .find((child) => child.type === 'heading')
      const columns = node.children[0]?.children.map((cell) => ctx.textContent(cell)).join(', ')
      const label = heading ? ctx.textContent(heading) : `Table: ${columns ?? ''}`
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
