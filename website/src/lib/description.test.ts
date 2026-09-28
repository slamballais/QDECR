import { test } from 'node:test'
import assert from 'node:assert/strict'
import { packageDescription, parseDescription } from './description.ts'

const SAMPLE = `Package: QDECR
Type: Package
Version: 0.9.0
Authors@R: as.person(c(
    "Ryan Muetzel [aut, cre]",
    "Sander Lamballais [aut]"
    ))
License: GPL-3
Imports:
    bigstatsr (>= 1.5.1),
    methods
`

test('reads single-line fields', () => {
  const fields = parseDescription(SAMPLE)
  assert.equal(fields['Package'], 'QDECR')
  assert.equal(fields['Version'], '0.9.0')
  assert.equal(fields['License'], 'GPL-3')
})

test('folds continuation lines into their field, one space apart', () => {
  const fields = parseDescription(SAMPLE)
  assert.equal(fields['Imports'], 'bigstatsr (>= 1.5.1), methods')
  assert.equal(
    fields['Authors@R'],
    'as.person(c( "Ryan Muetzel [aut, cre]", "Sander Lamballais [aut]" ))',
  )
})

test('copes with Windows line endings', () => {
  assert.equal(parseDescription(SAMPLE.replaceAll('\n', '\r\n'))['Version'], '0.9.0')
})

test("the package's own DESCRIPTION has a version the site can print", () => {
  // Run from website/, as npm test and the build both are.
  assert.match(packageDescription()['Version'] ?? '', /^\d+\.\d+\.\d+$/)
})
