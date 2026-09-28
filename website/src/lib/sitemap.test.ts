import { test } from 'node:test'
import assert from 'node:assert/strict'
import { sitemapXml } from './sitemap.ts'

const page = (file: string, head: string) => ({ file, html: `<!doctype html><html><head>${head}</head><body></body></html>` })

test('every page with a canonical link is listed by that URL, in order; noindex pages are not', () => {
  const xml = sitemapXml([
    page('tutorials/plotting.html', '<link rel="canonical" href="https://qdecr.com/tutorials/plotting">'),
    page('404.html', '<meta name="robots" content="noindex">'),
    page('index.html', '<link rel="canonical" href="https://qdecr.com/">'),
    page('cite.html', '<meta name="description" content="x"><link rel="canonical" href="https://qdecr.com/cite">'),
  ])
  assert.equal(
    xml,
    '<?xml version="1.0" encoding="UTF-8"?>\n' +
      '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n' +
      '  <url><loc>https://qdecr.com/</loc></url>\n' +
      '  <url><loc>https://qdecr.com/cite</loc></url>\n' +
      '  <url><loc>https://qdecr.com/tutorials/plotting</loc></url>\n' +
      '</urlset>\n',
  )
})
