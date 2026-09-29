// What R printed in the example analysis, pulled into the tutorials' code blocks rather
// than copied into them, so a page cannot show output the run did not produce:
//
//   ```r output=lh.stacks.txt
//   stacks(out)
//   ```
//
// appends the file from src/data/example/output/ (tools/example/export.R writes them)
// after the code, each line marked #> the way the site shows R's output (src/lib/shiki.ts
// mutes those lines, and the copy button leaves them out). lines=12-30 keeps a range of
// the file, counted from 1, both ends included, for a long log. Anything wrong, a missing
// file or a range the file does not have, fails the build.

import { defineMdastPlugin } from 'satteri'

interface Options {
  /** The text of an output file, by name; throws if there is none. */
  read: (name: string) => string
}

/** A file in the output directory itself: no path, so nothing outside it. */
const NAME = /^[\w.-]+$/

/** The fence's meta, word by word: output=lh.print.txt lines=1-20. */
function parseMeta(meta: string): Map<string, string> {
  const fields = new Map<string, string>()
  for (const word of meta.split(/\s+/)) {
    const eq = word.indexOf('=')
    if (eq > 0) fields.set(word.slice(0, eq), word.slice(eq + 1))
  }
  return fields
}

/** The output's lines, less the file's final newline and each line's trailing spaces. */
function outputLines(text: string, name: string, range: string | undefined): string[] {
  const lines = text.replace(/\n$/, '').split('\n').map((line) => line.trimEnd())
  if (range === undefined) return lines
  const match = /^(\d+)-(\d+)$/.exec(range)
  const from = Number(match?.[1])
  const to = Number(match?.[2])
  if (!match || from < 1 || to < from || to > lines.length) {
    throw new Error(`${name}: lines=${range} is not a range of its ${lines.length} lines (write it as from-to, counted from 1)`)
  }
  return lines.slice(from - 1, to)
}

export function exampleOutput({ read }: Options) {
  return defineMdastPlugin({
    name: 'qdecr:example-output',
    code(node, ctx) {
      const meta = parseMeta(node.meta ?? '')
      const output = meta.get('output')
      const lines = meta.get('lines')
      if (output === undefined) {
        if (lines !== undefined) throw new Error(`A code block has lines=${lines}, but lines= needs output= to say which file.`)
        return
      }
      if (!NAME.test(output)) throw new Error(`output=${output}: name a file in the example's output directory, without a path.`)
      const printed = outputLines(read(output), output, lines).map((line) => (line ? `#> ${line}` : '#>'))
      ctx.setProperty(node, 'value', [...(node.value ? [node.value] : []), ...printed].join('\n'))
    },
  })
}
