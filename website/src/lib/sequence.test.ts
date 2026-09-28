import { test } from 'node:test'
import assert from 'node:assert/strict'
import { neighbours, readingOrder } from './sequence.ts'

const PAGES = [
  { id: 'tutorials/plotting', title: 'Plotting', order: 5 },
  { id: 'get-started', title: 'Get started', order: 0 },
  { id: 'tutorials/quick-start', title: 'Quick start', order: 1 },
  { id: 'tutorials/formulas-and-design', title: 'Formulas and design', order: 2 },
]

test('pages are read in their order, whatever order they were found in', () => {
  assert.deepEqual(
    readingOrder(PAGES).map((page) => page.id),
    ['get-started', 'tutorials/quick-start', 'tutorials/formulas-and-design', 'tutorials/plotting'],
  )
})

test('a page in the middle links back and forward', () => {
  const { prev, next } = neighbours(PAGES, 'tutorials/quick-start')
  assert.equal(prev?.id, 'get-started')
  assert.equal(next?.id, 'tutorials/formulas-and-design')
})

test('the first page has no previous one, the last no next', () => {
  assert.equal(neighbours(PAGES, 'get-started').prev, undefined)
  assert.equal(neighbours(PAGES, 'tutorials/plotting').next, undefined)
})

test('two pages claiming the same place is an error, not a coin toss', () => {
  const clash = [...PAGES, { id: 'tutorials/saving', title: 'Saving', order: 5 }]
  assert.throws(() => readingOrder(clash), /tutorials\/plotting.*tutorials\/saving|tutorials\/saving.*tutorials\/plotting/)
})

test('asking for a page that is not in the sequence is an error', () => {
  assert.throws(() => neighbours(PAGES, 'tutorials/nope'), /tutorials\/nope/)
})
