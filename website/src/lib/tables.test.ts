import { test } from 'node:test'
import assert from 'node:assert/strict'
import { markdownToHtml } from 'satteri'
import { tables } from './tables.ts'

const labels = (markdown: string) =>
  [...markdownToHtml(markdown, { mdastPlugins: [tables()] }).html.matchAll(/aria-label="([^"]*)"/g)].map((m) => m[1])

const table = '| File | What it holds |\n|---|---|\n| a | b |\n'

test('a table is a focusable region named after the heading above it', () => {
  const html = markdownToHtml(`## Output\n\n${table}`, { mdastPlugins: [tables()] }).html
  assert.match(html, /<div class="table-wrap" tabindex="0" role="region" aria-label="Output">\s*<table>/)
})

test('without a heading above it, a table is named by its columns', () => {
  assert.deepEqual(labels(table), ['Table: File, What it holds'])
})

test('tables under the same heading are numbered, so no two regions share a name', () => {
  assert.deepEqual(labels(`## Output\n\nPer stack:\n\n${table}\nOverall:\n\n${table}\n## Next\n\n${table}`), [
    'Output',
    'Output, table 2',
    'Next',
  ])
})
