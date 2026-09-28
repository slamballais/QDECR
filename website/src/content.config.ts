// The site's Markdown content. Pages written in .astro (the home page, Cite, the reference)
// are not here.

import { defineCollection } from 'astro:content'
import { glob } from 'astro/loaders'
import { z } from 'astro/zod'
import { readFile } from 'node:fs/promises'
import { join } from 'node:path'
import { parseNews } from './lib/news'

/**
 * The guide: Get started and the tutorials, one reading sequence (src/lib/sequence.ts).
 * A file's path is its URL: src/content/docs/tutorials/plotting.md is /tutorials/plotting.
 */
const docs = defineCollection({
  loader: glob({ base: './src/content/docs', pattern: '**/*.md' }),
  schema: z.object({
    title: z.string(),
    /** One or two sentences, for search results, link previews and the tutorials index. */
    description: z.string(),
    /** Its place in the reading sequence: Get started is 0, tutorial n is n. */
    order: z.number().int().min(0),
  }),
})

/**
 * The changelog: the package's NEWS.md, one level above website/, as a single entry so
 * the whole file renders in one pass. That keeps the heading anchors unique: every
 * release has a "Bug fixes", and rendered apart they would all be #bug-fixes. Each release
 * gets its own <h2 id="v0.9.0">, which the page's "On this page" links to.
 */
const changelog = defineCollection({
  loader: {
    name: 'qdecr:news',
    load: async ({ store, renderMarkdown, parseData, generateDigest, watcher }) => {
      const path = join(process.cwd(), '..', 'NEWS.md')
      const load = async () => {
        const releases = parseNews(await readFile(path, 'utf8'))
        const body = releases
          .map((release) => {
            const name = release.name ? ` <span class="release-name">${escape(release.name)}</span>` : ''
            return `<h2 id="${release.id}">${escape(release.version)}${name}</h2>\n\n${release.body}`
          })
          .join('\n\n')
        const data = await parseData({
          id: 'news',
          data: { releases: releases.map(({ version, name, id }) => ({ version, name, id })) },
        })
        store.clear()
        store.set({ id: 'news', data, body, digest: generateDigest(body), rendered: await renderMarkdown(body) })
      }
      await load()
      // In dev, an edit to NEWS.md shows up without a restart.
      watcher?.add(path)
      watcher?.on('change', (changed) => {
        if (changed === path) void load()
      })
    },
  },
  schema: z.object({
    releases: z.array(z.object({ version: z.string(), name: z.string().optional(), id: z.string() })),
  }),
})

const escape = (text: string) => text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')

export const collections = { docs, changelog }
