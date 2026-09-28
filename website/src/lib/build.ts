// When the site was built, as the footer prints it.

import { execFileSync } from 'node:child_process'

/**
 * The date of the commit being built, not the time of the build: the site is rebuilt by
 * the monthly data refresh whether or not a page changed, and a build
 * time would claim every page changed each month. Works in Netlify's shallow clone; the
 * build time stands in where git is unavailable.
 */
export function commitDate(): Date {
  try {
    const iso = execFileSync('git', ['log', '-1', '--format=%cI'], {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
    })
    const date = new Date(iso.trim())
    if (!Number.isNaN(date.getTime())) return date
  } catch {
    // No git, or no repository: fall through.
  }
  return new Date()
}
