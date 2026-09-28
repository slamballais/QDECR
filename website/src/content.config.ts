// The site's Markdown content. Pages written in .astro (the home page, Cite, the reference)
// are not here.

import { defineCollection } from 'astro:content'
import { glob } from 'astro/loaders'
import { z } from 'astro/zod'

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

export const collections = { docs }
