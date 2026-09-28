import { test } from 'node:test'
import assert from 'node:assert/strict'
import { buildToc } from './toc.ts'

const h = (depth: number, text: string) => ({ depth, text, slug: text.toLowerCase().replaceAll(' ', '-') })

test('second-level headings make the list, third-level ones nest under them', () => {
  const toc = buildToc([h(2, 'Install'), h(3, 'Linux'), h(3, 'macOS'), h(2, 'Check it works')])
  assert.deepEqual(toc, [
    {
      slug: 'install',
      text: 'Install',
      children: [
        { slug: 'linux', text: 'Linux', children: [] },
        { slug: 'macos', text: 'macOS', children: [] },
      ],
    },
    { slug: 'check-it-works', text: 'Check it works', children: [] },
  ])
})

test('the page title and headings below the third level are left out', () => {
  const toc = buildToc([h(1, 'Title'), h(2, 'Section'), h(4, 'Detail')])
  assert.deepEqual(toc, [{ slug: 'section', text: 'Section', children: [] }])
})

test('a third-level heading before any second-level one stands on its own', () => {
  assert.deepEqual(buildToc([h(3, 'Aside'), h(2, 'Section')]).map((item) => item.slug), ['aside', 'section'])
})
