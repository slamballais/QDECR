// Runs Lighthouse on every page the build wrote, against the local server that serves dist/
// as Netlify will (tools/serve.ts), and holds them to the site's release checks:
//
//   npm run build && npm run serve        (in one terminal)
//   node tools/qa/lighthouse.ts [page...] (in another; all pages when none are named)
//
// It fails if a page scores under 95 in any category (bar SEO on the pages kept out of
// search engines, 404 and search, which score low there by design), requests anything
// from another origin, or, for the home page, goes over its budget before the reader asks
// for the 3D viewer: 100 KB for everything but the fonts, and 140 KB for the fonts, which
// the browser keeps for every page after (Sander's call, 29 Sep 2026). Lighthouse cannot
// open the viewer; tools/qa/browser.ts does.
//
// Lighthouse is fetched by npx at a pinned version rather than kept in package.json: it
// is a release check, and a large package with dependencies of its own that Netlify would
// otherwise install on every deploy.
// axe-core, which browser.ts reads into each page as a file, is a devDependency instead:
// npx runs a package's command, not a file inside it, and axe-core is 3 MB with no
// dependencies. The reports are written to tools/qa/reports/ (git-ignored).

import { exec } from 'node:child_process'
import { mkdir, readFile, stat } from 'node:fs/promises'
import { fileURLToPath } from 'node:url'
import { promisify } from 'node:util'
import { DIST, ORIGIN, pagesToCheck } from './pages.ts'

const LIGHTHOUSE = 'lighthouse@13.5.0'
const MIN_SCORE = 0.95
const HOME_BUDGET = { rest: 100 * 1024, fonts: 140 * 1024 }

const REPORTS = new URL('reports/', import.meta.url)
const run = promisify(exec)

interface Report {
  categories: Record<string, { title: string; score: number | null }>
  audits: Record<string, { details?: { items?: { url: string; transferSize?: number; resourceType?: string }[] } }>
}

const pages = await pagesToCheck()
await mkdir(REPORTS, { recursive: true })
const failures: string[] = []

for (const page of pages) {
  const out = new URL(`${page === '/' ? 'index' : page.slice(1).replaceAll('/', '_')}.json`, REPORTS)
  // One command line rather than an argument list: npx is a .cmd file on Windows, which
  // Node only runs through a shell. The page path is the one part from outside, the
  // command line, so it may hold nothing a shell would read.
  if (!/^\/[\w./-]*$/.test(page)) throw new Error(`Not a page path: ${page}`)
  const started = Date.now()
  await run(
    `npx --yes ${LIGHTHOUSE} "${ORIGIN}${page}" --output=json --output-path="${fileURLToPath(out)}" --chrome-flags=--headless=new --quiet`,
    { maxBuffer: 1 << 26 },
  ).catch(async (error: Error) => {
    // On Windows, Chrome can still hold its temporary profile when Lighthouse deletes it
    // at the end, and the run then fails with EPERM after the report is written. That
    // report is sound; any other failure is not.
    const written = await stat(out).then((file) => file.mtimeMs >= started, () => false)
    if (!(written && /EPERM/.test(error.message))) throw error
  })
  const report = JSON.parse(await readFile(out, 'utf8')) as Report
  const requests = report.audits['network-requests']?.details?.items ?? []

  const html = await readFile(new URL(page === '/' ? 'index.html' : `${page.slice(1)}.html`, DIST), 'utf8')
  const noindex = /<meta name="robots" content="noindex/.test(html)
  const scores = Object.entries(report.categories).map(([id, { score }]) => {
    const exempt = id === 'seo' && noindex
    if (!exempt && (score === null || score < MIN_SCORE)) failures.push(`${page}: ${id} scored ${score}`)
    return `${id} ${score === null ? '-' : Math.round(score * 100)}`
  })
  for (const { url } of requests) {
    if (!url.startsWith(ORIGIN) && !url.startsWith('data:')) failures.push(`${page}: requested ${url}`)
  }
  let weight = ''
  if (page === '/') {
    const sum = (font: boolean) =>
      requests
        .filter((request) => (request.resourceType === 'Font') === font)
        .reduce((total, request) => total + (request.transferSize ?? 0), 0)
    const kb = (bytes: number) => `${(bytes / 1024).toFixed(1)} KB`
    const fonts = sum(true)
    const rest = sum(false)
    weight = `, ${kb(rest)} + ${kb(fonts)} of fonts`
    if (rest >= HOME_BUDGET.rest) failures.push(`/: weighs ${kb(rest)} without its fonts`)
    if (fonts >= HOME_BUDGET.fonts) failures.push(`/: its fonts weigh ${kb(fonts)}`)
  }
  console.log(`${page}: ${scores.join(', ')}${weight}`)
}

if (failures.length) {
  console.error(`\n${failures.length} problem(s):\n${failures.join('\n')}`)
  process.exitCode = 1
}
