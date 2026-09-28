import { test } from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { parseNews } from './news.ts'

const SAMPLE = `# QDECR 0.9.0: Lausanne

Version 0.9.0 is the first update after publication.

## New features
* Weights.

## Bug fixes
* Cores.

# QDECR 0.8.5

## Bug fixes
* A fix.
`

test('each first-level heading starts a release, newest first, as NEWS.md lists them', () => {
  assert.deepEqual(
    parseNews(SAMPLE).map(({ version, name }) => ({ version, name })),
    [
      { version: '0.9.0', name: 'Lausanne' },
      { version: '0.8.5', name: undefined },
    ],
  )
})

test("a release's body is its text, with its headings moved down a level under the release", () => {
  const [latest] = parseNews(SAMPLE)
  assert.equal(
    latest!.body,
    'Version 0.9.0 is the first update after publication.\n\n### New features\n* Weights.\n\n### Bug fixes\n* Cores.',
  )
})

test('an anchor per release, safe to use as an id', () => {
  assert.deepEqual(
    parseNews(SAMPLE).map((release) => release.id),
    ['v0.9.0', 'v0.8.5'],
  )
})

test('R comments inside a code fence are code, not headings', () => {
  const news = '# QDECR 1.0.0\n\n```r\n# fit the model\nqdecr_fastlm()\n```\n'
  assert.equal(parseNews(news)[0]!.body, '```r\n# fit the model\nqdecr_fastlm()\n```')
})

test('Windows line endings make no difference', () => {
  assert.deepEqual(parseNews(SAMPLE.replaceAll('\n', '\r\n')), parseNews(SAMPLE))
})

test('a first-level heading that is not a release is an error, not a silent merge', () => {
  assert.throws(() => parseNews('# Changes\n\n* Something.\n'), /# Changes/)
})

test("the package's own NEWS.md parses, with the latest release first", () => {
  const releases = parseNews(readFileSync(new URL('../../../NEWS.md', import.meta.url), 'utf8'))
  assert.ok(releases.length >= 8)
  assert.match(releases[0]!.version, /^\d+\.\d+\.\d+$/)
  assert.ok(releases.some((release) => release.name === 'Momo'))
})
