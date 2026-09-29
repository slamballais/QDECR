// Every crawler may read everything; the sitemap (src/integrations/sitemap.ts) lists the
// pages worth indexing, and pages that are not carry their own noindex.

import type { APIRoute } from 'astro'
import { site } from '../data/site'

export const GET: APIRoute = () =>
  new Response(`User-agent: *\nAllow: /\n\nSitemap: ${new URL('/sitemap.xml', site.url).href}\n`, {
    headers: { 'Content-Type': 'text/plain; charset=utf-8' },
  })
