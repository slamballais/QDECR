// Draws each page's share card (tools/og.mjs) once the build has written the page. Every
// indexable page names its card in <meta property="og:image"> (layouts/Base.astro), with
// its title and label in the metas beside it; this reads those back out of the finished
// HTML, so the card always says what the page says. Pages without a card (noindex) are
// skipped.

import type { AstroIntegration } from 'astro'
import { mkdir, writeFile } from 'node:fs/promises'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { createCardRenderer } from '../../tools/og.mjs'
import { builtPages } from './built-pages.ts'

/** The content of a <meta> named `key` (by property or name), entities decoded. */
function meta(html: string, key: string): string | undefined {
  const tag = new RegExp(`<meta (?:property|name)="${key}" content="([^"]*)"`).exec(html)
  return tag?.[1]
    ?.replace(/&quot;/g, '"')
    .replace(/&#39;|&#x27;/g, "'")
    .replace(/&lt;/g, '<')
    .replace(/&gt;/g, '>')
    .replace(/&amp;/g, '&')
}

export default function ogImages(): AstroIntegration {
  return {
    name: 'qdecr:og-images',
    hooks: {
      'astro:build:done': async ({ dir, logger }) => {
        const root = fileURLToPath(dir)
        const render = await createCardRenderer()
        let drawn = 0
        for (const { html } of await builtPages(dir)) {
          const image = meta(html, 'og:image')
          const title = meta(html, 'og:title')
          if (!image || !title) continue
          const out = join(root, new URL(image).pathname)
          await mkdir(dirname(out), { recursive: true })
          const section = meta(html, 'qdecr:card-section')
          await writeFile(out, await render({ title, section, code: meta(html, 'qdecr:card-code') === 'true' }))
          drawn += 1
        }
        logger.info(`${drawn} share cards drawn`)
      },
    },
  }
}
