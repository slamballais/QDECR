import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { GROUPS, groupTopics, type Topic } from './reference.ts'

const topic = (name: string) => ({ name, title: name }) as Topic

const SMALL = [
  { label: 'Run an analysis', topics: ['qdecr_fastlm', 'qdecr'] },
  { label: 'Plot', topics: ['qdecr_snap'] },
]

test('topics come out in their groups, in the order the groups list them', () => {
  const grouped = groupTopics([topic('qdecr_snap'), topic('qdecr'), topic('qdecr_fastlm')], SMALL)
  assert.deepEqual(
    grouped.map((group) => [group.label, group.topics.map((t) => t.name)]),
    [
      ['Run an analysis', ['qdecr_fastlm', 'qdecr']],
      ['Plot', ['qdecr_snap']],
    ],
  )
})

test('a help page no group lists fails, so a new function cannot go missing from the reference', () => {
  assert.throws(
    () => groupTopics([topic('qdecr_fastlm'), topic('qdecr'), topic('qdecr_snap'), topic('qdecr_new')], SMALL),
    /qdecr_new/,
  )
})

test('a group listing a help page that does not exist fails, so a rename is noticed', () => {
  assert.throws(() => groupTopics([topic('qdecr_fastlm'), topic('qdecr')], SMALL), /qdecr_snap/)
})

test('a help page listed in two groups fails', () => {
  const twice = [...SMALL, { label: 'Again', topics: ['qdecr'] }]
  assert.throws(() => groupTopics([topic('qdecr_fastlm'), topic('qdecr'), topic('qdecr_snap')], twice), /qdecr.*twice|twice.*qdecr/)
})

test('every help page in reference.json is in exactly one of the seven groups', () => {
  const { topics } = JSON.parse(readFileSync(new URL('../data/reference.json', import.meta.url), 'utf8')) as {
    topics: Topic[]
  }
  const grouped = groupTopics(topics)
  assert.deepEqual(
    grouped.map((group) => group.label),
    [
      'Run an analysis',
      'Inspect results',
      'Plot',
      'Save and load',
      'MGH and annotation I/O',
      'Imputation helpers',
      'Low-level FBM helpers',
    ],
  )
  assert.equal(grouped.flatMap((group) => group.topics).length, topics.length)
  assert.equal(GROUPS.length, 7)
})
