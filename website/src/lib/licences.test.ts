import { test } from 'node:test'
import assert from 'node:assert/strict'
import { checkNotices, packageOf } from './licences.ts'

test('a module under node_modules belongs to its package, scoped or not', () => {
  assert.equal(packageOf('C:/site/node_modules/fflate/esm/browser.js'), 'fflate')
  assert.equal(packageOf('/site/node_modules/@niivue/niivue/dist/index.js'), '@niivue/niivue')
  assert.equal(packageOf('C:\\site\\node_modules\\gl-matrix\\esm\\vec3.js'), 'gl-matrix')
})

test('a package nested in another belongs to the innermost one', () => {
  assert.equal(packageOf('/site/node_modules/@niivue/niivue/node_modules/zarrita/dist/index.js'), 'zarrita')
})

test("a query on the module's id does not change its package", () => {
  assert.equal(packageOf('/site/node_modules/fflate/esm/browser.js?commonjs-es-import'), 'fflate')
})

test("the bundlers' own helpers belong to the bundler", () => {
  assert.equal(packageOf('\0vite/preload-helper.js'), 'vite')
  assert.equal(packageOf('\0rolldown/runtime.js'), 'rolldown')
})

test("the site's own modules belong to no package", () => {
  assert.equal(packageOf('/site/src/scripts/viewer.ts'), undefined)
  assert.equal(packageOf('/site/src/components/Viewer.astro?astro&type=script&index=0&lang.ts'), undefined)
})

const notices = [
  { packages: ['@niivue/niivue'], bundled: true },
  { packages: ['fflate', 'gl-matrix'], bundled: true },
  { packages: ['pagefind'], bundled: false },
]

test('every package in the bundle needs a notice', () => {
  const ids = ['/s/node_modules/fflate/a.js', '/s/node_modules/zarrita/b.js', '/s/node_modules/zarrita/c.js', '/s/src/x.ts']
  assert.deepEqual(checkNotices(ids, notices).missing, ['zarrita'])
})

test('a notice for a bundled package that is no longer in the bundle is stale', () => {
  const ids = ['/s/node_modules/fflate/a.js', '/s/node_modules/@niivue/niivue/b.js']
  assert.deepEqual(checkNotices(ids, notices).stale, ['gl-matrix'])
})

test('a notice for something shipped outside the bundle is never stale', () => {
  assert.deepEqual(checkNotices([], notices).stale, ['@niivue/niivue', 'fflate', 'gl-matrix'])
})
