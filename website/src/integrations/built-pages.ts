// The pages a build wrote, for the integrations that check or draw from them.

import { readdir, readFile } from 'node:fs/promises'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'

export interface BuiltPage {
  /** Its path under dist/: tutorials/plotting.html. */
  file: string
  html: string
}

/** Every HTML page under the build's output directory. */
export async function builtPages(dir: URL): Promise<BuiltPage[]> {
  const root = fileURLToPath(dir)
  const files = (await readdir(root, { recursive: true })).filter((file) => file.endsWith('.html'))
  return Promise.all(files.map(async (file) => ({ file, html: await readFile(join(root, file), 'utf8') })))
}
