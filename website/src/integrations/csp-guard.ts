// Fails the build if any page it wrote carries markup the Content-Security-Policy blocks:
// an inline script, a <style> element, a style attribute or an event handler (see
// src/lib/inline.ts for why the pages are checked rather than the source).

import type { AstroIntegration } from 'astro'
import { findInline } from '../lib/inline.ts'
import { builtPages } from './built-pages.ts'

export default function cspGuard(): AstroIntegration {
  return {
    name: 'qdecr:csp-guard',
    hooks: {
      'astro:build:done': async ({ dir, logger }) => {
        const pages = await builtPages(dir)
        const failures = pages.flatMap(({ file, html }) => findInline(html).map((problem) => `${file}: ${problem}`))
        if (failures.length) {
          throw new Error(
            `${failures.length} thing(s) the CSP in netlify.toml would block:\n${failures.join('\n')}`,
          )
        }
        logger.info(`${pages.length} pages checked: nothing inline for the CSP to block`)
      },
    },
  }
}
