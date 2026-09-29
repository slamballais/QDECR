// The publications that cite the QDECR paper, as tools/fetch-citations.ts stores them in
// src/data/citations.json from OpenAlex. The Cite page reads that file.

export interface CitingWork {
  /** OpenAlex's id, without the https://openalex.org/ prefix: W4289317516. */
  id: string
  /** Without the https://doi.org/ prefix, or null where OpenAlex has none. */
  doi: string | null
  title: string
  year: number
  authors: string[]
  /** The journal, repository or proceedings, where OpenAlex knows it. */
  venue: string | null
  /** OpenAlex's type: article, preprint, dissertation, ... */
  type: string
}

/** The fields of an OpenAlex work that are read here. */
export interface OpenAlexWork {
  id: string
  doi: string | null
  display_name: string | null
  publication_year: number
  authorships: { author: { display_name: string } }[]
  primary_location: { source: { display_name: string } | null } | null
  type: string
}

const tidy = (text: string) => text.replace(/\s+/g, ' ').trim()
const bare = (url: string, prefix: string) => (url.startsWith(prefix) ? url.slice(prefix.length) : url)

/**
 * OpenAlex results as the Cite page needs them: once each, newest first and then by
 * title, so the committed file changes only when the list does. Titles lose any markup
 * (OpenAlex keeps the <i> of species names and the like).
 */
export function normaliseCitations(results: readonly OpenAlexWork[]): CitingWork[] {
  const byId = new Map<string, CitingWork>()
  for (const work of results) {
    const id = bare(work.id, 'https://openalex.org/')
    if (byId.has(id) || !work.display_name) continue
    byId.set(id, {
      id,
      doi: work.doi ? bare(work.doi, 'https://doi.org/') : null,
      title: tidy(work.display_name.replace(/<[^>]+>/g, '')),
      year: work.publication_year,
      authors: work.authorships.map((authorship) => tidy(authorship.author.display_name)),
      venue: work.primary_location?.source?.display_name ?? null,
      type: work.type,
    })
  }
  return [...byId.values()].sort((a, b) => b.year - a.year || a.title.localeCompare(b.title, 'en'))
}

/** A citing work kept off the site, and why, so the list can be reviewed later. */
export interface HiddenCitation {
  id: string
  reason: string
}

/**
 * The works to list: OpenAlex's, less the ones hidden by hand (src/data/citations-hidden.ts).
 * OpenAlex counts some works twice (a preprint and its published version) and counts the
 * paper as citing itself; it has no way to say so, so the fix lives here.
 */
export function hideCitations(works: readonly CitingWork[], hidden: readonly HiddenCitation[]): CitingWork[] {
  const ids = new Set(hidden.map((entry) => entry.id))
  return works.filter((work) => !ids.has(work.id))
}
