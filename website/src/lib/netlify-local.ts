// How Netlify answers a request, from the rules in netlify.toml, for the local server that
// serves dist/ as production will (tools/serve.ts). `astro preview` sends no CSP and
// compresses nothing, so Lighthouse, the link check and the browser checks before a
// release (tools/qa/) run against this instead. Only what this site's netlify.toml uses
// is covered: exact redirects on this host, header rules for a path or a folder, and
// pretty URLs.

import { posix } from 'node:path'

export interface RedirectRule {
  from: string
  to: string
  status?: number
  force?: boolean
}

export interface HeaderRule {
  for: string
  values: Record<string, string>
}

export interface Redirect {
  to: string
  status: number
  force: boolean
}

/**
 * The redirect for a path, if a rule names it. Rules for another host (the www ones) never
 * match: the local server has only the one.
 */
export function redirectFor(path: string, rules: readonly RedirectRule[]): Redirect | undefined {
  const rule = rules.find((candidate) => candidate.from.startsWith('/') && candidate.from === path)
  if (!rule) return undefined
  return { to: rule.to, status: rule.status ?? 301, force: rule.force ?? false }
}

/** The headers for a path: every rule that matches it, in order, a later value winning. */
export function headersFor(path: string, rules: readonly HeaderRule[]): Record<string, string> {
  const values: Record<string, string> = {}
  for (const rule of rules) {
    const folder = rule.for.endsWith('/*') ? rule.for.slice(0, -1) : undefined
    if (folder ? path.startsWith(folder) : path === rule.for) Object.assign(values, rule.values)
  }
  return values
}

/**
 * The files under dist/ that could answer a path, best first. A path with an extension is
 * that file; one without is a file, then the page with .html (build.format 'file' writes
 * /about as about.html), then a folder's index. Nothing for a path that climbs out of the
 * site or cannot be decoded.
 */
export function candidateFiles(path: string): string[] {
  let decoded: string
  try {
    decoded = decodeURIComponent(path)
  } catch {
    return []
  }
  // A backslash is a separator to Windows, so %5c..%5c would climb out of dist/ there. No
  // URL of this site has one.
  if (decoded.includes('\\') || decoded.split('/').includes('..')) return []
  const file = posix.normalize(decoded).replace(/^\/+/, '')
  if (file === '') return ['index.html']
  if (decoded.endsWith('/')) return [`${file.replace(/\/$/, '')}/index.html`, `${file.replace(/\/$/, '')}.html`]
  if (posix.extname(file)) return [file]
  return [file, `${file}.html`, `${file}/index.html`]
}
