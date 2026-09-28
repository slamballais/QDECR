// The function reference: the help pages as tools/rd-to-json.R writes them to
// src/data/reference.json, and the groups the reference is read in.

/** One help page (an Rd file), which may document several functions (its aliases). */
export interface Topic {
  name: string
  /** Its URL: /reference/<slug>. The name with dots as hyphens. */
  slug: string
  aliases: string[]
  /** False for a help page whose functions the package does not export. */
  exported: boolean
  title: string
  /** The R file whose roxygen comments the page is generated from. */
  source: string | null
  /** The text sections are HTML (a small subset); usage and examples are R code. */
  description: string | null
  usage: string | null
  arguments: { name: string; html: string | null }[]
  value: string | null
  details: string | null
  sections: { title: string; html: string | null }[]
  note: string | null
  author: string | null
  references: string | null
  seealso: string | null
  examples: string | null
  keywords: string[]
}

export interface Group {
  label: string
  /** Help pages by name, in reading order. */
  topics: string[]
}

/**
 * How the reference is grouped. Every help page must be in exactly one
 * group: a new one fails the tests and the build until it is placed here.
 */
export const GROUPS: Group[] = [
  { label: 'Run an analysis', topics: ['qdecr_fastlm', 'qdecr'] },
  {
    label: 'Inspect results',
    topics: ['print.vw', 'summary.vw_fastlm', 'hist.vw', 'stacks', 'formula.vw_fastlm', 'nobs.vw_fastlm', 'qdecr_fwhm'],
  },
  { label: 'Plot', topics: ['qdecr_snap', 'freeview'] },
  { label: 'Save and load', topics: ['qdecr_save', 'qdecr_load', 'unload', 'reload'] },
  { label: 'MGH and annotation I/O', topics: ['qdecr_read', 'load.mgh', 'save.mgh', 'as_mgh', 'load.annot'] },
  { label: 'Imputation helpers', topics: ['imp2list'] },
  { label: 'Low-level FBM helpers', topics: ['qdecr_prep_mgh', 'bsfbm2mgh', 'fbm_functions_from_qdecr'] },
]

/** The topics sorted into their groups. Throws on a topic in no group, or in two, or a
 * group naming a topic that does not exist. */
export function groupTopics<T extends Pick<Topic, 'name'>>(topics: readonly T[], groups: Group[] = GROUPS) {
  const byName = new Map(topics.map((topic) => [topic.name, topic]))
  const placed = new Set<string>()
  const result = groups.map((group) => ({
    label: group.label,
    topics: group.topics.map((name) => {
      const topic = byName.get(name)
      if (!topic) throw new Error(`reference group "${group.label}" lists ${name}, which has no help page`)
      if (placed.has(name)) throw new Error(`${name} is listed twice in the reference groups`)
      placed.add(name)
      return topic
    }),
  }))
  const unplaced = topics.filter((topic) => !placed.has(topic.name)).map((topic) => topic.name)
  if (unplaced.length) {
    throw new Error(`help pages in no reference group (add them to GROUPS in src/lib/reference.ts): ${unplaced.join(', ')}`)
  }
  return result
}
