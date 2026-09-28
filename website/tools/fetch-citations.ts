// Fetches the publications that cite the QDECR paper from OpenAlex and writes them to
// src/data/citations.json, for the Cite page:
//
//   npm run fetch:citations
//
// It runs monthly in GitHub Actions (.github/workflows/site-data.yml), which opens a pull
// request when the list changes; the file is committed so the build never calls OpenAlex.
// OpenAlex needs no key for this, and nothing about the reader is sent.

import { readFile, writeFile } from 'node:fs/promises'
import { site } from '../src/data/site.ts'
import { normaliseCitations, type OpenAlexWork } from '../src/lib/citations.ts'

const API = 'https://api.openalex.org'
const OUT = new URL('../src/data/citations.json', import.meta.url)
const FIELDS = 'id,doi,display_name,publication_year,authorships,primary_location,type'

async function get<T>(path: string): Promise<T> {
  const response = await fetch(`${API}${path}`, {
    headers: { 'User-Agent': 'qdecr.com citation sync (https://github.com/slamballais/QDECR)' },
  })
  if (!response.ok) throw new Error(`OpenAlex answered ${response.status} for ${path}`)
  return (await response.json()) as T
}

const doi = site.paper.url.replace('https://doi.org/', '')
const paper = await get<{ id: string }>(`/works/doi:${doi}?select=id`)
const paperId = paper.id.replace('https://openalex.org/', '')

// Cursor paging, 200 at a time: the most OpenAlex returns per page.
const results: OpenAlexWork[] = []
let cursor: string | null = '*'
while (cursor) {
  const page: { results: OpenAlexWork[]; meta: { next_cursor: string | null } } = await get(
    `/works?filter=cites:${paperId}&per-page=200&select=${FIELDS}&cursor=${encodeURIComponent(cursor)}`,
  )
  results.push(...page.results)
  cursor = page.results.length ? page.meta.next_cursor : null
}
const works = normaliseCitations(results)

// An empty answer where there used to be a list is far likelier to be an OpenAlex hiccup
// than every citing paper vanishing; refuse it rather than open a pull request that
// empties the Cite page.
const previous = await readFile(OUT, 'utf8').then(
  (text) => (JSON.parse(text) as { works: unknown[] }).works.length,
  () => 0,
)
if (works.length === 0 && previous > 0) {
  throw new Error(`OpenAlex returned no citing works; citations.json has ${previous}. Not overwriting it.`)
}

await writeFile(OUT, `${JSON.stringify({ paper: { openalex: paperId, doi }, works }, null, 2)}\n`)
console.log(`${works.length} works cite ${doi} (OpenAlex ${paperId}); ${previous} before.`)
