// Fails the build if the code it sends to the browser comes from a package whose licence
// /licences does not print, or if /licences names a bundled package the bundle no longer
// carries. Most licences the site's dependencies use ask for their notice to travel with
// the code, and a bundle strips it, so the page is where it travels (src/lib/licences.ts,
// the notices in src/data/licences.ts).

import type { AstroIntegration } from 'astro'
import { notices } from '../data/licences.ts'
import { checkNotices } from '../lib/licences.ts'

export default function licenceGuard(): AstroIntegration {
  // Every module in the client bundle: the scripts pages load, and what they import.
  const moduleIds = new Set<string>()
  return {
    name: 'qdecr:licence-guard',
    hooks: {
      'astro:config:setup': ({ updateConfig }) => {
        updateConfig({
          vite: {
            plugins: [
              {
                name: 'qdecr:client-modules',
                generateBundle(_options, bundle) {
                  // Astro also bundles the pages themselves, for the build only; none of
                  // that reaches the browser.
                  if (this.environment.name !== 'client') return
                  for (const output of Object.values(bundle)) {
                    if (output.type === 'chunk') for (const id of output.moduleIds) moduleIds.add(id)
                  }
                },
              },
            ],
          },
        })
      },
      'astro:build:done': ({ logger }) => {
        const { missing, stale } = checkNotices(moduleIds, notices)
        const problems = [
          ...missing.map((name) => `${name} is in the bundle, but /licences has no notice for it`),
          ...stale.map((name) => `${name} has a notice on /licences, but is no longer in the bundle`),
        ]
        if (problems.length) {
          throw new Error(`Fix src/data/licences.ts:\n${problems.join('\n')}`)
        }
        logger.info(`${notices.length} licence notices cover every package in the bundle`)
      },
    },
  }
}
