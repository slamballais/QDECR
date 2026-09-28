// The Glossary page's terms, read from src/data/glossary.md: the project's glossary,
// which the code and the docs both follow. The file is written for people first:
//
//   ### Group
//
//   **Term**:
//   The definition, in Markdown, over one or more lines.
//   _Avoid_: other words, for the same thing
//
// so the page takes what it needs from it rather than the other way round.

import { readFileSync } from 'node:fs'
import { join } from 'node:path'

export interface GlossaryTerm {
  term: string
  /** Markdown, on one line. */
  definition: string
  /** Words used elsewhere for the same thing, which the docs do not use. */
  avoid: string[]
}

export interface GlossaryGroup {
  label: string
  terms: GlossaryTerm[]
}

export interface Glossary {
  groups: GlossaryGroup[]
}

const SECTION = /^##\s+(.+)$/
const GROUP = /^###\s+(.+)$/
const TERM = /^\*\*(.+)\*\*:\s*$/
const AVOID = /^_Avoid_:\s*(.*)$/

/** The terms under the file's "## Language" heading; its other sections are not glossary. */
export function parseGlossary(markdown: string): Glossary {
  const groups: GlossaryGroup[] = []
  let inLanguage = false
  let term: GlossaryTerm | undefined
  for (const line of markdown.split(/\r?\n/)) {
    const section = SECTION.exec(line)
    if (section) {
      inLanguage = section[1]!.trim() === 'Language'
      term = undefined
      continue
    }
    if (!inLanguage) continue
    const group = GROUP.exec(line)
    if (group) {
      groups.push({ label: group[1]!.trim(), terms: [] })
      term = undefined
      continue
    }
    const start = TERM.exec(line)
    if (start) {
      term = { term: start[1]!.trim(), definition: '', avoid: [] }
      if (!groups.length) groups.push({ label: '', terms: [] })
      groups.at(-1)!.terms.push(term)
      continue
    }
    if (!term || !line.trim()) {
      term = undefined
      continue
    }
    const avoid = AVOID.exec(line)
    if (avoid) {
      term.avoid = avoid[1]!.split(',').map((word) => word.trim()).filter(Boolean)
      continue
    }
    term.definition = `${term.definition} ${line.trim()}`.trim()
  }
  return { groups }
}

/**
 * The glossary in src/data/glossary.md. Resolved from the
 * working directory for the reason packageDescription() gives (src/lib/description.ts).
 */
export function projectGlossary(): Glossary {
  return parseGlossary(readFileSync(join(process.cwd(), 'src', 'data', 'glossary.md'), 'utf8'))
}
