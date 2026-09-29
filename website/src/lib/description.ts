// The package's DESCRIPTION file, read at build time. The site documents the release on
// master, so the version in the footer and anywhere else comes from here
// rather than from a number typed into the site.

import { readFileSync } from 'node:fs'
import { join } from 'node:path'

/**
 * Parses R's DESCRIPTION format (Debian control-file style): `Field: value` lines, with
 * any line that starts with whitespace continuing the field above it. Continuation lines
 * are folded into one line, a single space apart.
 */
export function parseDescription(text: string): Record<string, string> {
  const fields: Record<string, string> = {}
  let current: string | undefined
  for (const line of text.split(/\r?\n/)) {
    if (/^\s/.test(line)) {
      if (current && line.trim()) fields[current] = `${fields[current]} ${line.trim()}`.trim()
      continue
    }
    const match = /^([^:\s]+):\s*(.*)$/.exec(line)
    if (!match) continue
    current = match[1]!
    fields[current] = match[2]!.trim()
  }
  return fields
}

/**
 * The fields of the DESCRIPTION at the repo root, one level above website/. Resolved from
 * the working directory, not from this file: Astro bundles this module into dist/ before
 * running it, so import.meta.url no longer points into src/. npm scripts, and Netlify
 * with its base directory set to website/, both run from website/.
 */
export function packageDescription(): Record<string, string> {
  return parseDescription(readFileSync(join(process.cwd(), '..', 'DESCRIPTION'), 'utf8'))
}
