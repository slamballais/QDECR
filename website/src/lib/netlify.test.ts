/**
 * The promises netlify.toml makes, checked every run. Run with `npm test`.
 *
 * Nothing fails visibly when a header or a redirect goes missing: the page still loads,
 * and the link someone put in a paper in 2020 quietly lands on a 404. So the rules the
 * site depends on are held here.
 */

import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { parse } from 'smol-toml'

interface Redirect {
  from: string
  to: string
  status?: number
  force?: boolean
}

interface Headers {
  for: string
  values: Record<string, string>
}

const config = parse(readFileSync(new URL('../../netlify.toml', import.meta.url), 'utf8')) as {
  build?: { command?: string; publish?: string; environment?: Record<string, string> }
  redirects?: Redirect[]
  headers?: Headers[]
}
const redirects = config.redirects ?? []
const headers = config.headers ?? []

const redirectFrom = (from: string) => redirects.find((rule) => rule.from === from)
const headersFor = (path: string) => headers.find((rule) => rule.for === path)?.values ?? {}

test('builds with npm into dist/, on the Node version package.json asks for', () => {
  assert.equal(config.build?.command, 'npm run build')
  assert.equal(config.build?.publish, 'dist')
  const pkg = JSON.parse(readFileSync(new URL('../../package.json', import.meta.url), 'utf8'))
  const floor = /^>=(\d+\.\d+\.\d+)$/.exec(pkg.engines.node)?.[1]
  assert.equal(config.build?.environment?.['NODE_VERSION'], floor)
})

for (const from of ['https://www.qdecr.com/*', 'http://www.qdecr.com/*']) {
  test(`${from} redirects to the apex for good`, () => {
    const rule = redirectFrom(from)
    assert.ok(rule, `no redirect from ${from}`)
    assert.equal(rule.to, 'https://qdecr.com/:splat')
    assert.equal(rule.status, 301)
    assert.equal(rule.force, true)
  })
}

// Every URL of the 2021 site (the website branch), and where its content lives now.
const OLD_URLS: Record<string, string> = {
  '/01-getting-started.html': '/get-started',
  '/02-quick-start.html': '/tutorials/quick-start',
  '/03-using-qdecr.html': '/tutorials/formulas-and-design',
  '/03-post-processing.html': '/tutorials/inspecting-results',
  '/04-post-processing.html': '/tutorials/inspecting-results',
  '/about.html': '/about',
  '/contribution.html': '/help#contributing',
  '/code-of-conduct.html': '/help#code-of-conduct',
  '/ohbm2020.html': '/about#archive',
  '/data/qdecr_lmm.pdf': '/archive/qdecr_lmm.pdf',
  '/data/qdecr_ohbm2019_poster.pdf': '/archive/qdecr_ohbm2019_poster.pdf',
  '/data/qdecr_ohbm_brochure.pdf': '/archive/qdecr_ohbm_brochure.pdf',
  // Never existed, but the old site linked to it; the brochure is the file it meant.
  '/data/qdecr_ohbm2019_brochure.pdf': '/archive/qdecr_ohbm_brochure.pdf',
  // Netlify served the old pages without .html too, and the OHBM 2020 poster prints
  // qdecr.com/ohbm2020. (/about is a page of the new site, so needs no rule.)
  '/01-getting-started': '/get-started',
  '/02-quick-start': '/tutorials/quick-start',
  '/03-using-qdecr': '/tutorials/formulas-and-design',
  '/03-post-processing': '/tutorials/inspecting-results',
  '/04-post-processing': '/tutorials/inspecting-results',
  '/contribution': '/help#contributing',
  '/code-of-conduct': '/help#code-of-conduct',
  '/ohbm2020': '/about#archive',
}

for (const [from, to] of Object.entries(OLD_URLS)) {
  test(`${from} moves permanently to ${to}`, () => {
    const rule = redirectFrom(from)
    assert.ok(rule, `no redirect from ${from}`)
    assert.equal(rule.to, to)
    assert.equal(rule.status, 301)
    // Forced, because build.format 'file' writes /about as about.html: without force,
    // Netlify finds that file at /about.html and serves it instead of redirecting.
    assert.equal(rule.force, true, `${from} is not forced`)
  })
}

test('the CSP allows this origin only and nothing inline', () => {
  const csp = headersFor('/*')['Content-Security-Policy']
  assert.ok(csp, 'no Content-Security-Policy on /*')
  assert.doesNotMatch(csp, /unsafe-inline/)
  assert.doesNotMatch(csp, /'unsafe-eval'/, "only 'wasm-unsafe-eval' is allowed, for Pagefind")
  const directives = Object.fromEntries(
    csp.split(';').map((part) => {
      const [name, ...values] = part.trim().split(/\s+/)
      return [name, values]
    }),
  )
  assert.deepEqual(directives['default-src'], ["'self'"])
  assert.deepEqual(directives['object-src'], ["'none'"])
  assert.deepEqual(directives['frame-ancestors'], ["'none'"])
  assert.deepEqual(directives['base-uri'], ["'none'"])
  // Every source that is not a keyword must be this origin's: no third-party requests.
  for (const [name, values] of Object.entries(directives)) {
    for (const value of values) {
      assert.match(value, /^('[a-z-]+'|data:|blob:)$/, `${name} allows ${value}`)
    }
  }
})

test("images may be data: URLs, for the viewer, and nothing else is widened for it", () => {
  const csp = headersFor('/*')['Content-Security-Policy'] ?? ''
  const directive = (name: string) =>
    csp
      .split(';')
      .map((part) => part.trim().split(/\s+/))
      .find(([key]) => key === name)
      ?.slice(1)
  // NiiVue draws its font and its lighting from images it carries as data: URLs.
  assert.deepEqual(directive('img-src'), ["'self'", 'data:'])
  // Its WebGL needs nothing here, and it runs no workers: scripts stay as Pagefind needs.
  assert.deepEqual(directive('script-src'), ["'self'", "'wasm-unsafe-eval'"])
  assert.equal(directive('worker-src'), undefined)
})

test('the other security headers are set on every page', () => {
  const values = headersFor('/*')
  assert.equal(values['X-Frame-Options'], 'DENY')
  assert.equal(values['X-Content-Type-Options'], 'nosniff')
  assert.equal(values['Referrer-Policy'], 'strict-origin-when-cross-origin')
  assert.match(values['Permissions-Policy'] ?? '', /camera=\(\)/)
})

test('hashed build assets are cached for a year', () => {
  assert.equal(headersFor('/_astro/*')['Cache-Control'], 'public, max-age=31536000, immutable')
})
