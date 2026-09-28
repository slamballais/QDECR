// The search page: queries the Pagefind index that `npm run build` writes to dist/pagefind/
// and lists the results. Results update as the reader types, and the query is kept in the
// address bar so a search can be linked to and survives going back.

import { servedPath } from '../lib/path'

/** The part of Pagefind's browser API used here. */
interface Pagefind {
  options(options: { excerptLength?: number }): Promise<void>
  debouncedSearch(query: string): Promise<{ results: { data(): Promise<ResultData> }[] } | null>
}

interface ResultData {
  url: string
  excerpt: string
  meta: { title?: string; section?: string }
}

const form = document.querySelector<HTMLFormElement>('.search-form')!
const input = document.querySelector<HTMLInputElement>('#search-query')!
const statusLine = document.querySelector<HTMLElement>('.search-status')!
const list = document.querySelector<HTMLOListElement>('.search-results')!

const MAX_RESULTS = 20

let pagefind: Pagefind | undefined

async function load(): Promise<Pagefind | undefined> {
  if (pagefind) return pagefind
  try {
    // A path, not a module Vite can see: the index only exists after the build.
    const url = '/pagefind/pagefind.js'
    pagefind = (await import(/* @vite-ignore */ url)) as Pagefind
    await pagefind.options({ excerptLength: 24 })
    return pagefind
  } catch {
    statusLine.textContent = 'Search is not available here: the index is built with the site (npm run build).'
    return undefined
  }
}


function render(results: ResultData[], total: number, query: string) {
  list.replaceChildren(
    ...results.map((result) => {
      const item = document.createElement('li')
      const link = document.createElement('a')
      link.href = servedPath(result.url)
      link.textContent = result.meta.title ?? servedPath(result.url)
      const where = document.createElement('p')
      where.className = 'search-section'
      where.textContent = result.meta.section ?? ''
      const excerpt = document.createElement('p')
      // Pagefind's excerpt is text from the site's own pages with the matches in <mark>.
      excerpt.innerHTML = result.excerpt
      item.append(link, ...(result.meta.section ? [where] : []), excerpt)
      return item
    }),
  )
  statusLine.textContent =
    total === 0
      ? `Nothing found for “${query}”.`
      : total > results.length
        ? `${total} pages match “${query}”; the first ${results.length} are shown.`
        : `${total} ${total === 1 ? 'page matches' : 'pages match'} “${query}”.`
}

async function run(query: string) {
  const trimmed = query.trim()
  const url = new URL(window.location.href)
  if (trimmed) url.searchParams.set('q', trimmed)
  else url.searchParams.delete('q')
  history.replaceState(null, '', url)

  if (!trimmed) {
    list.replaceChildren()
    statusLine.textContent = ''
    return
  }
  const search = await (await load())?.debouncedSearch(trimmed)
  // null: a newer query has superseded this one.
  if (!search) return
  const results = await Promise.all(search.results.slice(0, MAX_RESULTS).map((result) => result.data()))
  render(results, search.results.length, trimmed)
}

form.addEventListener('submit', (event) => {
  event.preventDefault()
  void run(input.value)
})
input.addEventListener('input', () => void run(input.value))

const initial = new URLSearchParams(window.location.search).get('q') ?? ''
input.value = initial
if (initial) void run(initial)
input.focus()
