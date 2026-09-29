// Which packages' code the site sends to the browser, so that each one's licence is printed
// on /licences. The build lists the modules in its client bundles and fails when one comes
// from a package without a notice there, or when a notice names a package the bundle no
// longer carries (src/integrations/licence-guard.ts). The notices are src/data/licences.ts.

/** What the check needs from a notice: its packages, and whether the bundler ships them. */
export interface NoticeScope {
  packages: readonly string[]
  /** False for code and fonts the site ships as files (Pagefind, the fonts), which no bundle lists. */
  bundled: boolean
}

const NODE_MODULES = /.*node_modules\/((?:@[^/]+\/)?[^/]+)\//

/**
 * The package a module of the bundle comes from, or undefined for the site's own code. The
 * bundlers' own helpers have ids starting with a NUL byte: \0vite/preload-helper.js.
 */
export function packageOf(id: string): string | undefined {
  const path = id.replaceAll('\\', '/').replace(/\?.*$/, '')
  const helper = /^\0([^/]+)\//.exec(path)
  if (helper) return helper[1]
  return NODE_MODULES.exec(path)?.[1]
}

/**
 * The packages in the bundle without a notice (missing), and those with a notice for
 * bundled code that the bundle no longer carries (stale), each sorted.
 */
export function checkNotices(
  moduleIds: Iterable<string>,
  notices: readonly NoticeScope[],
): { missing: string[]; stale: string[] } {
  const shipped = new Set<string>()
  for (const id of moduleIds) {
    const name = packageOf(id)
    if (name) shipped.add(name)
  }
  const noticed = new Set(notices.flatMap((notice) => notice.packages))
  const bundled = notices.filter((notice) => notice.bundled).flatMap((notice) => notice.packages)
  return {
    missing: [...shipped].filter((name) => !noticed.has(name)).sort(),
    stale: bundled.filter((name) => !shipped.has(name)).sort(),
  }
}
