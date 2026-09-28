// Fails the build if any page it wrote carries markup the Content-Security-Policy blocks:
// an inline script, a <style> element, a style attribute or an event handler (see
// src/lib/inline.ts for why the pages are checked rather than the source).

import type { AstroIntegration } from 'astro'
import { readdir, readFile } from 'node:fs/promises'
import { fileURLToPath } from 'node:url'
import { join, relative } from 'node:path'
import { findInline } from '../lib/inline.ts'

export default function cspGuard(): AstroIntegration {
  return {
    name: 'qdecr:csp-guard',
    hooks: {
      'astro:build:done': async ({ dir, logger }) => {
        const root = fileURLToPath(dir)
        const files = (await readdir(root, { recursive: true })).filter((file) => file.endsWith('.html'))
        const failures: string[] = []
        for (const file of files) {
          const problems = findInline(await readFile(join(root, file), 'utf8'))
          for (const problem of problems) failures.push(`${relative(root, join(root, file))}: ${problem}`)
        }
        if (failures.length) {
          throw new Error(
            `${failures.length} thing(s) the CSP in netlify.toml would block:\n${failures.join('\n')}`,
          )
        }
        logger.info(`${files.length} pages checked: nothing inline for the CSP to block`)
      },
    },
  }
}
