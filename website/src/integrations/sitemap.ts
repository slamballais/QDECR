// Writes dist/sitemap.xml once the build has written the pages (src/lib/sitemap.ts). Built
// here rather than with @astrojs/sitemap: the pages already say which of them belong in it,
// through their canonical links, and reading that back is a dozen lines.

import type { AstroIntegration } from 'astro'
import { writeFile } from 'node:fs/promises'
import { sitemapXml } from '../lib/sitemap.ts'
import { builtPages } from './built-pages.ts'

export default function sitemap(): AstroIntegration {
  return {
    name: 'qdecr:sitemap',
    hooks: {
      'astro:build:done': async ({ dir, logger }) => {
        const xml = sitemapXml(await builtPages(dir))
        await writeFile(new URL('sitemap.xml', dir), xml)
        logger.info(`sitemap.xml lists ${xml.match(/<url>/g)?.length ?? 0} pages`)
      },
    },
  }
}
