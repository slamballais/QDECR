/**
 * The path a page is served at, from the one Astro reports. With build.format 'file',
 * Astro.url.pathname is /cite.html at build time but /cite in dev; Netlify serves /cite.
 * Canonical links, current-page markers and share images all need the served form, or
 * they quietly point at the .html variant in production.
 */
export function servedPath(pathname: string): string {
  return pathname.replace(/(\/index)?\.html$/, '').replace(/\/$/, '') || '/'
}
