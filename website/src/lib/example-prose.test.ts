// The tutorials quote the example run's numbers in their prose, where the Markdown cannot
// read run.json the way the output blocks and the home page do. These tests hold each
// quoted number to the committed run, so an export that changes one fails here, naming
// what to rewrite, rather than leave a page describing a different analysis.

import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { parseExampleRun } from './example.ts'
import { count, signed } from './format.ts'

const run = parseExampleRun(JSON.parse(readFileSync(new URL('../data/example/run.json', import.meta.url), 'utf8')))
const { dataset, hemispheres, model, poster } = run
const age = (hemi: 'lh' | 'rh') => {
  const cluster = hemispheres[hemi].clusters.find((c) => c.stack === 'age')
  assert.ok(cluster, `the ${hemi} run has no age cluster, which the tutorials describe`)
  return cluster
}

const tutorial = (name: string) => readFileSync(new URL(`../content/docs/tutorials/${name}.md`, import.meta.url), 'utf8')

/** Every phrase must appear in the page as written. */
function quotes(name: string, phrases: string[]) {
  const text = tutorial(name)
  const missing = phrases.filter((phrase) => !text.includes(phrase))
  assert.deepEqual(missing, [], `${name}.md no longer matches run.json; rewrite: ${missing.join(' | ')}`)
}

test('the quick start quotes the sample, the clusters and the run as run.json has them', () => {
  const lh = age('lh')
  const rh = age('rh')
  quotes('quick-start', [
    `${dataset.n} people aged ${dataset.age.min} to ${dataset.age.max}, ${dataset.sex.female} of them female`,
    `${count(lh.nVertices)} vertices and ${count(lh.sizeMm2)} mm²`,
    `${count(rh.nVertices)} vertices and ${count(rh.sizeMm2)} mm²`,
    `mean coefficient is ${signed(lh.meanCoefficient, 3)}`,
    `mean coefficient of ${signed(rh.meanCoefficient, 3)}`,
    `came out at ${hemispheres.lh.smoothness} mm`,
    `took ${count(hemispheres.lh.seconds)} seconds`,
    `n_cores = ${model.nCores}`,
  ])
})

test('Inspecting results quotes the left hemisphere as run.json has it', () => {
  const lh = age('lh')
  const top = lh.regions[0]!
  quotes('inspecting-results', [
    `${dataset.n} people`,
    `${count(lh.nVertices)} of the ${count(hemispheres.lh.vertices.analysed)} vertices`,
    `thickness drops ${signed(-lh.meanCoefficient, 3)} mm a year`,
    `${top.ofCluster.toFixed(1)}%, lies in`,
    `covers ${Math.round(top.ofRegion)}% of that region`,
    `${hemispheres.lh.smoothness} mm`,
  ])
})

test("Plotting quotes the poster's colour scale as run.json has it", () => {
  quotes('plotting', [`${dataset.n} people`, `overlay_threshold = c(${poster.scale.from}, ${poster.scale.to})`])
})
