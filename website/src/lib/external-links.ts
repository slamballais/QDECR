// Links that leave the site, in Markdown: each gets the out icon after its text, as the
// .astro pages add <OutIcon /> by hand (plan 3.1). A link counts as leaving when it is an
// absolute http(s) URL to any host but qdecr.com; paths, anchors and the site's own
// canonical URLs stay plain.

import { defineMdastPlugin } from 'satteri'
import { OUT_ICON } from './out-icon.ts'

const OUT = /^https?:\/\/(?!(www\.)?qdecr\.com(\/|$))/i

export function externalLinks() {
  return defineMdastPlugin({
    name: 'qdecr:external-links',
    link(node, ctx) {
      if (!OUT.test(node.url)) return
      ctx.appendChild(node, { type: 'html', value: OUT_ICON })
    },
  })
}
