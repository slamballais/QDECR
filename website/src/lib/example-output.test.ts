import { test } from 'node:test'
import assert from 'node:assert/strict'
import { markdownToHtml } from 'satteri'
import { exampleOutput } from './example-output.ts'

const FILES: Record<string, string> = {
  'lh.stacks.txt': '[1] "(Intercept)" "age"         "sexmale"    \n',
  'lh.log.txt': 'one\ntwo\n\nfour\nfive\n',
}

const read = (name: string) => {
  const text = FILES[name]
  if (text === undefined) throw new Error(`ENOENT: ${name}`)
  return text
}

/** The text of the rendered code block, entities decoded; the fence has no highlighter here. */
const code = (markdown: string) =>
  markdownToHtml(markdown, { mdastPlugins: [exampleOutput({ read })] })
    .html.replace(/^<pre><code[^>]*>|\n?<\/code><\/pre>\s*$/g, '')
    .replace(/&gt;/g, '>')
    .replace(/&quot;/g, '"')

test('a fence naming an output file gets that output after its code, as #> lines', () => {
  assert.equal(code('```r output=lh.stacks.txt\nstacks(out)\n```'), 'stacks(out)\n#> [1] "(Intercept)" "age"         "sexmale"')
})

test('an empty line of output stays a line, with the marker alone', () => {
  assert.equal(code('```r output=lh.log.txt\nx\n```'), 'x\n#> one\n#> two\n#>\n#> four\n#> five')
})

test('lines= keeps a range of the output, counted from 1, both ends included', () => {
  assert.equal(code('```r output=lh.log.txt lines=2-4\nx\n```'), 'x\n#> two\n#>\n#> four')
})

test('a fence with no code shows the output alone', () => {
  assert.equal(code('```r output=lh.log.txt lines=4-5\n```'), '#> four\n#> five')
})

test('a fence without output= is left as it is', () => {
  assert.equal(code('```r\nout\n```'), 'out')
  assert.equal(code('```r title=x\nout\n```'), 'out')
})

test('a file that is not there fails the build, naming it', () => {
  assert.throws(() => code('```r output=lh.nothing.txt\nx\n```'), /lh\.nothing\.txt/)
})

test('a range outside the output, or back to front, fails the build', () => {
  assert.throws(() => code('```r output=lh.log.txt lines=4-9\nx\n```'), /lines=4-9.*5 lines/)
  assert.throws(() => code('```r output=lh.log.txt lines=3-2\nx\n```'), /lines=3-2/)
  assert.throws(() => code('```r output=lh.log.txt lines=0-2\nx\n```'), /lines=0-2/)
})

test('lines= without output=, or a name that leaves the directory, fails the build', () => {
  assert.throws(() => code('```r lines=1-2\nx\n```'), /lines= needs output=/)
  assert.throws(() => code('```r output=../secret.txt\nx\n```'), /\.\.\/secret\.txt/)
})
