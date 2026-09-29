// sitemap.xml, from the pages the build wrote (src/integrations/sitemap.ts). Each page is
// listed by its own canonical link, which layouts/Base.astro writes on every indexable page
// and leaves off noindex ones, so the sitemap and the pages cannot disagree about either.

import type { BuiltPage } from '../integrations/built-pages.ts'

const LINK = /<link\b[^>]*>/g

/** The href of the page's <link rel="canonical">, if it has one. */
function canonical(html: string): string | undefined {
  for (const [tag] of html.matchAll(LINK)) {
    if (/\brel="canonical"/.test(tag)) return /\bhref="([^"]*)"/.exec(tag)?.[1]
  }
  return undefined
}

/** The sitemap, sorted by URL so the file only changes when the pages do. */
export function sitemapXml(pages: readonly BuiltPage[]): string {
  const urls = pages
    .map((page) => canonical(page.html))
    .filter((url): url is string => url !== undefined)
    .sort()
  return (
    '<?xml version="1.0" encoding="UTF-8"?>\n' +
    '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n' +
    urls.map((url) => `  <url><loc>${url}</loc></url>\n`).join('') +
    '</urlset>\n'
  )
}
