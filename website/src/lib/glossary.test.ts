import { test } from 'node:test'
import assert from 'node:assert/strict'
import { parseGlossary, projectGlossary } from './glossary.ts'
import { headingId } from './toc.ts'

const SAMPLE = `# QDECR

Vertex-wise statistics on FreeSurfer surfaces, in R.

## Language

### Surfaces

**Vertex**:
A point of the surface mesh, where every measure is sampled.
_Avoid_: voxel, node

**Vertex measure**:
The surface map analysed at every vertex, such as cortical
thickness, named in the formula as \`qdecr_thickness\`.

### Results

**Stack**:
One column of the design matrix, with the maps that belong to it.
_Avoid_: contrast
`

test('terms are read with their definitions, under the group they sit in', () => {
  const glossary = parseGlossary(SAMPLE)
  assert.deepEqual(
    glossary.groups.map((group) => [group.label, group.terms.map((term) => term.term)]),
    [
      ['Surfaces', ['Vertex', 'Vertex measure']],
      ['Results', ['Stack']],
    ],
  )
  assert.equal(glossary.groups[0]!.terms[0]!.definition, 'A point of the surface mesh, where every measure is sampled.')
})

test('a definition over several lines is one line, and the words to avoid are a list of their own', () => {
  const [surfaces, results] = parseGlossary(SAMPLE).groups
  assert.deepEqual(surfaces!.terms[1], {
    term: 'Vertex measure',
    definition: 'The surface map analysed at every vertex, such as cortical thickness, named in the formula as `qdecr_thickness`.',
    avoid: [],
  })
  assert.deepEqual(surfaces!.terms[0]!.avoid, ['voxel', 'node'])
  assert.deepEqual(results!.terms[0], {
    term: 'Stack',
    definition: 'One column of the design matrix, with the maps that belong to it.',
    avoid: ['contrast'],
  })
})

test('only the Language section is read; a term outside a group gets a group with no label', () => {
  const glossary = parseGlossary(`# QDECR

## Language

**Hemisphere**:
The left or right half of the cortex, analysed one at a time.

## Flagged ambiguities

**Target**:
Used for both the template subject and its directory.
`)
  assert.deepEqual(glossary.groups, [
    { label: '', terms: [{ term: 'Hemisphere', definition: 'The left or right half of the cortex, analysed one at a time.', avoid: [] }] },
  ])
})

test('src/data/glossary.md parses: every term has a definition, and every name a unique anchor', () => {
  const { groups } = projectGlossary()
  const terms = groups.flatMap((group) => group.terms)
  assert.ok(terms.length > 20, `only ${terms.length} terms read`)
  assert.deepEqual(terms.filter((term) => !term.definition).map((term) => term.term), [])
  const anchors = [...groups.map((group) => group.label), ...terms.map((term) => term.term)].map(headingId)
  assert.deepEqual(anchors.filter((anchor, i) => anchors.indexOf(anchor) !== i), [])
})
