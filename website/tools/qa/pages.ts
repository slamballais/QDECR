// What the release checks share: where the built site is served, and which pages it has.

import { readdir } from 'node:fs/promises'

/** The local server that serves dist/ as Netlify will (tools/serve.ts). */
export const ORIGIN = process.env['QA_ORIGIN'] ?? 'http://localhost:8888'

export const DIST = new URL('../../dist/', import.meta.url)

/** Every page under dist/, as the path it is served at: about.html is /about. */
export async function allPages(): Promise<string[]> {
  const files = await readdir(DIST, { recursive: true })
  return files
    .map((file) => file.replaceAll('\\', '/'))
    .filter((file) => file.endsWith('.html') && !file.startsWith('pagefind/'))
    .map((file) => '/' + file.replace(/(^|\/)index\.html$/, '').replace(/\.html$/, ''))
    .sort()
}

/** The pages named on the command line, or all of them. */
export async function pagesToCheck(): Promise<string[]> {
  return process.argv.length > 2 ? process.argv.slice(2) : allPages()
}
