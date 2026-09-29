// Removes the parts of Pagefind's output the site never loads, after `pagefind --site dist`
// has written them (the last step of `npm run build`):
//
//   pagefind-ui.*, pagefind-component-ui.*, pagefind-modular-ui.*, pagefind-highlight.js
//
// Pagefind always writes its ready-made search interfaces, about 400 KB with Svelte inside
// one of them. The site queries the index through pagefind.js with a search page of its own
// (src/scripts/search.ts), because those interfaces set inline styles the CSP forbids. So
// they would only be deployed, never requested, and carry code whose licence /licences
// does not print. What stays: pagefind.js, its worker and WebAssembly, the entry file, and
// the index and fragments.

import { readdir, rm } from 'node:fs/promises'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'

const DIR = fileURLToPath(new URL('../dist/pagefind/', import.meta.url))
const UNUSED = /^pagefind-(ui|component-ui|modular-ui|highlight)\.(js|css)$/

const removed = (await readdir(DIR)).filter((name) => UNUSED.test(name))
await Promise.all(removed.map((name) => rm(join(DIR, name))))
console.log(`Removed Pagefind's unused interfaces: ${removed.join(', ') || 'none found'}`)
