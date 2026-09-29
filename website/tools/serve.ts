// Serves the built site the way Netlify will, for the checks before a release:
//
//   npm run build && npm run serve
//
// then http://localhost:8888 (the port `netlify dev` uses). Unlike `astro preview`, every
// response carries the headers in netlify.toml (so the CSP is enforced, and a page that
// breaks it breaks here too), the old URLs redirect, /about finds about.html, text is
// compressed as Netlify compresses it, and an unknown path gets dist/404.html with a 404.
// The rules are in src/lib/netlify-local.ts; this file only reads and sends.

import { createServer } from 'node:http'
import { readFile, stat } from 'node:fs/promises'
import { extname, join, sep } from 'node:path'
import { fileURLToPath } from 'node:url'
import { brotliCompressSync, gzipSync } from 'node:zlib'
import { parse } from 'smol-toml'
import {
  candidateFiles,
  headersFor,
  redirectFor,
  type HeaderRule,
  type RedirectRule,
} from '../src/lib/netlify-local.ts'

const ROOT = new URL('../', import.meta.url)
const DIST = fileURLToPath(new URL('dist/', ROOT))
const PORT = Number(process.env['PORT'] ?? 8888)

const config = parse(await readFile(new URL('netlify.toml', ROOT), 'utf8')) as {
  redirects?: RedirectRule[]
  headers?: HeaderRule[]
}
const redirects = config.redirects ?? []
const headerRules = config.headers ?? []

const TYPES: Record<string, string> = {
  '.html': 'text/html; charset=UTF-8',
  '.css': 'text/css; charset=UTF-8',
  '.js': 'application/javascript; charset=UTF-8',
  '.json': 'application/json',
  '.txt': 'text/plain; charset=UTF-8',
  '.xml': 'application/xml',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.webp': 'image/webp',
  '.avif': 'image/avif',
  '.ico': 'image/x-icon',
  '.woff2': 'font/woff2',
  '.woff': 'font/woff',
  '.pdf': 'application/pdf',
  '.wasm': 'application/wasm',
}
// Netlify compresses text; images, fonts and the gzipped surfaces are compressed already.
const COMPRESSIBLE = /^(text\/|application\/(javascript|json|xml)|image\/svg)/

async function findFile(path: string): Promise<string | undefined> {
  for (const candidate of candidateFiles(path)) {
    const file = join(DIST, candidate)
    // candidateFiles refuses paths that climb; this makes sure of it.
    if (!file.startsWith(DIST.endsWith(sep) ? DIST : DIST + sep)) continue
    const found = await stat(file).catch(() => undefined)
    if (found?.isFile()) return file
  }
  return undefined
}

// It listens on this machine only (127.0.0.1): a server for checks has no business being
// reachable from the network.
createServer(async (request, response) => {
  const url = new URL(request.url ?? '/', `http://localhost:${PORT}`)
  const path = url.pathname
  const redirect = redirectFor(path, redirects)
  let file = redirect?.force ? undefined : await findFile(path)

  const headers: Record<string, string> = {
    // Netlify's own default, which the rules in netlify.toml then override.
    'Cache-Control': 'public, max-age=0, must-revalidate',
    ...headersFor(path, headerRules),
  }
  if (redirect && !file) {
    response.writeHead(redirect.status, { ...headers, Location: redirect.to }).end()
    return
  }
  let status = 200
  if (!file) {
    status = 404
    file = join(DIST, '404.html')
  }

  const type = TYPES[extname(file)] ?? 'application/octet-stream'
  let body = await readFile(file)
  headers['Content-Type'] = type
  if (COMPRESSIBLE.test(type)) {
    headers['Vary'] = 'Accept-Encoding'
    const accepts = String(request.headers['accept-encoding'] ?? '')
    if (/\bbr\b/.test(accepts)) {
      body = brotliCompressSync(body)
      headers['Content-Encoding'] = 'br'
    } else if (/\bgzip\b/.test(accepts)) {
      body = gzipSync(body)
      headers['Content-Encoding'] = 'gzip'
    }
  }
  headers['Content-Length'] = String(body.length)
  response.writeHead(status, headers).end(request.method === 'HEAD' ? undefined : body)
}).listen(PORT, '127.0.0.1', () => console.log(`dist/ as Netlify serves it: http://localhost:${PORT}`))
