// The package's NEWS.md, split into releases for the changelog page. NEWS.md stays the
// one place release notes are written; the site only reads it.

export interface Release {
  /** 0.9.0 */
  version: string
  /** The release's name, where it has one: Lausanne, Momo. */
  name: string | undefined
  /** An anchor for the release's heading: v0.9.0. */
  id: string
  /** The release's notes as Markdown, their headings one level down (## becomes ###). */
  body: string
}

const RELEASE = /^# QDECR (\d+\.\d+\.\d+)(?::\s*(.+?))?\s*$/

/**
 * Releases in the order NEWS.md lists them, which is newest first. Every first-level
 * heading must be a release ("# QDECR 0.9.0" or "# QDECR 0.9.0: Lausanne"); anything
 * else would otherwise be folded into the release above it without a word.
 */
export function parseNews(text: string): Release[] {
  const releases: Release[] = []
  let lines: string[] | undefined
  // Inside a code fence a leading # is an R or shell comment, not a heading.
  let fenced = false
  const finish = () => {
    if (lines && releases.length) releases.at(-1)!.body = lines.join('\n').trim()
  }
  for (const line of text.split(/\r?\n/)) {
    if (/^(```|~~~)/.test(line)) fenced = !fenced
    if (!fenced && line.startsWith('# ')) {
      const match = RELEASE.exec(line)
      if (!match) throw new Error(`NEWS.md: "${line}" is a first-level heading but not a release`)
      finish()
      const version = match[1]!
      releases.push({ version, name: match[2], id: `v${version}`, body: '' })
      lines = []
      continue
    }
    lines?.push(!fenced && line.startsWith('#') ? `#${line}` : line)
  }
  finish()
  return releases
}
