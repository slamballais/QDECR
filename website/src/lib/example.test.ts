import { test } from 'node:test'
import assert from 'node:assert/strict'
import { existsSync, readFileSync } from 'node:fs'
import { parseExampleRun } from './example.ts'

const hemisphere = (hemi: 'lh' | 'rh') => ({
  project: `${hemi}.age_sex.thickness`,
  vertices: { loaded: 163842, analysed: 149955 },
  fwhmEstimate: 16,
  seconds: 90.5,
  stacks: [
    { number: 1, name: '(Intercept)' },
    { number: 2, name: 'age' },
    { number: 3, name: 'sexmale' },
  ],
  clusters: [
    {
      stack: 'age',
      cluster: 1,
      nVertices: 5000,
      sizeMm2: 3200.5,
      cwp: 0.0001,
      peak: { value: -5.2, vertex: 12345, region: 'superiorfrontal' },
      meanThickness: 2.71,
      meanCoefficient: -0.021,
      meanSe: 0.0035,
      regions: [
        { name: 'superiorfrontal', ofCluster: 45.2, ofRegion: 12.1 },
        { name: 'rostralmiddlefrontal', ofCluster: 30.1, ofRegion: 9.8 },
      ],
    },
  ],
})

const run = {
  date: '2026-09-28',
  dataset: {
    name: 'ABIDE I',
    site: 'NYU',
    n: 99,
    sex: { female: 26, male: 73 },
    age: { min: 6.5, max: 31.8, mean: 15.6, median: 14.4 },
    excluded: [{ id: 'NYU_0051059', reason: 'anatomical scan failed the PCP quality check' }],
  },
  software: {
    qdecr: '0.9.0',
    r: 'R version 4.5.2 (2025-10-31)',
    freesurfer: 'freesurfer-linux-ubuntu22_amd64-7.4.1-20230614-7eb8460',
    os: 'Ubuntu 26.04 LTS',
    platform: 'WSL2 on Windows 11',
  },
  model: { formula: 'qdecr_thickness ~ age + sex', measure: 'thickness', fwhm: 10, mczThr: 30, cwpThr: 0.025, nCores: 4 },
  hemispheres: { lh: hemisphere('lh'), rh: hemisphere('rh') },
}

test('a complete run parses as it is, with both hemispheres', () => {
  assert.deepEqual(parseExampleRun(run), run)
})

test('a run needs both hemispheres: the default cwp_thr already splits 0.05 over the two', () => {
  assert.throws(() => parseExampleRun({ ...run, hemispheres: { lh: hemisphere('lh') } }), /rh/)
})

test('the sample has to add up: the sexes sum to n', () => {
  const dataset = { ...run.dataset, sex: { female: 26, male: 74 } }
  assert.throws(() => parseExampleRun({ ...run, dataset }), /add up/)
})

test('the model is the one the example fixes: thickness on age and sex, nothing else', () => {
  const model = { ...run.model, formula: 'qdecr_thickness ~ age + sex + site' }
  assert.throws(() => parseExampleRun({ ...run, model }), /the formula must be/)
})

test('clusters are numbered from 1 within their stack, in order', () => {
  const lh = hemisphere('lh')
  lh.clusters = [{ ...lh.clusters[0]!, cluster: 2 }]
  assert.throws(() => parseExampleRun({ ...run, hemispheres: { lh, rh: hemisphere('rh') } }), /numbered/)
})

const committed = new URL('../data/example/run.json', import.meta.url)
test('the committed run.json is a valid run', { skip: !existsSync(committed) && 'no run.json yet: step 5 has not run' }, () => {
  const parsed = parseExampleRun(JSON.parse(readFileSync(committed, 'utf8')))
  assert.equal(parsed.dataset.site, 'NYU')
})
