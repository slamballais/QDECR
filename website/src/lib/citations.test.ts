import { test } from 'node:test'
import assert from 'node:assert/strict'
import { hideCitations, normaliseCitations, type CitingWork } from './citations.ts'

// Trimmed from what https://api.openalex.org/works?filter=cites:W3158050040 returns.
const author = (name: string) => ({ author: { id: 'https://openalex.org/A1', display_name: name } })
const RESULTS = [
  {
    id: 'https://openalex.org/W4289317516',
    doi: 'https://doi.org/10.1001/jamanetworkopen.2022.24701',
    display_name: 'Association of Maternal Tobacco Use During Pregnancy With Preadolescent Brain Morphology',
    publication_year: 2022,
    authorships: [author('Runyu Zou'), author('Ryan L. Muetzel')],
    primary_location: { source: { display_name: 'JAMA Network Open' } },
    type: 'article',
  },
  {
    id: 'https://openalex.org/W4000000001',
    doi: null,
    display_name: 'Cortical thickness in <i>adolescents</i>',
    publication_year: 2024,
    authorships: [author('A. Author')],
    primary_location: null,
    type: 'preprint',
  },
  {
    id: 'https://openalex.org/W4000000002',
    doi: 'https://doi.org/10.1000/abc',
    display_name: 'Another 2024 paper',
    publication_year: 2024,
    authorships: [],
    primary_location: { source: null },
    type: 'article',
  },
]

test('each work keeps what the Cite page prints, with bare ids and DOIs', () => {
  const [, , jama] = normaliseCitations(RESULTS)
  assert.deepEqual(jama, {
    id: 'W4289317516',
    doi: '10.1001/jamanetworkopen.2022.24701',
    title: 'Association of Maternal Tobacco Use During Pregnancy With Preadolescent Brain Morphology',
    year: 2022,
    authors: ['Runyu Zou', 'Ryan L. Muetzel'],
    venue: 'JAMA Network Open',
    type: 'article',
  })
})

test('newest first, then by title, so the file only changes when the list does', () => {
  assert.deepEqual(
    normaliseCitations(RESULTS).map((work) => work.id),
    ['W4000000002', 'W4000000001', 'W4289317516'],
  )
})

test('markup in a title is dropped, and a missing DOI or venue is null', () => {
  const work = normaliseCitations(RESULTS).find((w) => w.id === 'W4000000001')!
  assert.equal(work.title, 'Cortical thickness in adolescents')
  assert.equal(work.doi, null)
  assert.equal(work.venue, null)
})

test('a work listed twice (two pages of results) is kept once', () => {
  assert.equal(normaliseCitations([...RESULTS, RESULTS[0]!]).length, 3)
})

test("stray spaces inside an author's name are closed up", () => {
  const [work] = normaliseCitations([{ ...RESULTS[0]!, authorships: [author('Charlotte AM  Cecil ')] }])
  assert.deepEqual(work!.authors, ['Charlotte AM Cecil'])
})

const work = (id: string, title: string, type = 'article'): CitingWork => ({
  id,
  doi: null,
  title,
  year: 2022,
  authors: ['A. Author'],
  venue: null,
  type,
})

test('hidden works are left out, and the rest keep their order', () => {
  const works = [
    work('W4308637228', 'ADHD symptoms and brain morphology'),
    work('W4225394107', 'ADHD symptoms and brain morphology', 'preprint'),
    work('W4310281878', 'Poverty, cortical structure'),
    work('W3158050040', 'QDECR: A Flexible, Extensible Vertex-Wise Analysis Framework in R'),
  ]
  const hidden = [
    { id: 'W4225394107', reason: 'Preprint of W4308637228.' },
    { id: 'W3158050040', reason: 'The QDECR paper itself.' },
  ]
  assert.deepEqual(
    hideCitations(works, hidden).map((w) => w.id),
    ['W4308637228', 'W4310281878'],
  )
})
