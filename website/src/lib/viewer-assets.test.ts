// The home page viewer's files, as tools/example/export.R writes them to public/viewer/:
// per hemisphere, fsaverage6's inflated surface, where its sulci are, and the poster's map
// (the age stack's −log10(p) on its significant clusters), each a gzipped MZ3. Nothing
// checks them in the browser until someone opens the viewer, so they are held here to
// what the viewer asks for and to the run they come from.

import { test } from 'node:test'
import assert from 'node:assert/strict'
import { existsSync, readFileSync, statSync } from 'node:fs'
import { gunzipSync } from 'node:zlib'
import { parseExampleRun } from './example.ts'
import { HEMISPHERES, viewerFiles } from './viewer.ts'

const run = parseExampleRun(JSON.parse(readFileSync(new URL('../data/example/run.json', import.meta.url), 'utf8')))
const publicDir = new URL('../../public/', import.meta.url)
/** Where a file the site serves at `url` lives on disk, in public/. */
const inPublic = (url: string) => new URL(url.replace(/^\//, ''), publicDir)

/** fsaverage6: an icosahedron divided six times. */
const VERTICES = 40962
const TRIANGLES = 81920
/** fsaverage, whose maps export.R cuts down to fsaverage6's first 40,962 vertices. */
const FSAVERAGE_VERTICES = 163842

interface Mz3 {
  faces: Uint32Array | null
  vertices: Float32Array | null
  values: Float32Array | null
  nVertices: number
}

// The reading half of export.R's write_mz3: a 16-byte header, then each part the bit
// field names, little-endian, the whole gzipped. NiiVue reads the same.
function readMz3(url: string): Mz3 {
  const file = readFileSync(inPublic(url))
  assert.deepEqual([...file.subarray(0, 2)], [0x1f, 0x8b], `${url} is not gzipped`)
  const bytes = gunzipSync(file)
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength)
  assert.equal(view.getUint16(0, true), 23117, `${url} is not MZ3`)
  const attributes = view.getUint16(2, true)
  const nFaces = view.getUint32(4, true)
  const nVertices = view.getUint32(8, true)
  let at = 16 + view.getUint32(12, true)
  // Each part is copied out, since it need not be aligned for a typed array over the
  // buffer. Typed arrays read in the machine's byte order, little-endian on anything
  // Node runs on.
  function part<T>(bit: number, length: number, Type: new (buffer: ArrayBuffer) => T): T | null {
    if ((attributes & bit) === 0) return null
    const start = bytes.byteOffset + at
    at += length * 4
    return new Type(bytes.buffer.slice(start, start + length * 4))
  }
  const faces = part(1, nFaces * 3, Uint32Array)
  const vertices = part(2, nVertices * 3, Float32Array)
  assert.equal(attributes & 4, 0, `${url} has colours, which the viewer does not use`)
  const values = part(8, nVertices, Float32Array)
  assert.equal(at, bytes.byteLength, `${url} is longer than its header says`)
  return { faces, vertices, values, nVertices }
}

for (const hemi of HEMISPHERES) {
  const files = viewerFiles(hemi)

  test(`${hemi}: every file the viewer asks for is there`, () => {
    for (const url of Object.values(files)) {
      assert.ok(existsSync(inPublic(url)), `no ${url} in public/`)
    }
  })

  // The viewer's size budget. Opening the viewer downloads one
  // hemisphere's files, about 0.9 MB as gzipped MZ3; a re-export that doubled that
  // would go unnoticed on a fast connection.
  test(`${hemi}: opening the viewer costs under 1 MB of surfaces and maps`, () => {
    const bytes = Object.values(files).reduce((sum, url) => sum + statSync(inPublic(url)).size, 0)
    assert.ok(bytes < 1_000_000, `${hemi}'s files come to ${(bytes / 1e6).toFixed(2)} MB`)
  })

  test(`${hemi}: the surface is fsaverage6's, whole`, () => {
    const mesh = readMz3(files.mesh)
    assert.equal(mesh.nVertices, VERTICES)
    assert.equal(mesh.faces?.length, TRIANGLES * 3)
    assert.equal(mesh.vertices?.length, VERTICES * 3)
    assert.equal(mesh.values, null)
    // Every triangle names real vertices, and every vertex is in a triangle.
    const used = new Uint8Array(VERTICES)
    for (const v of mesh.faces ?? []) {
      assert.ok(v < VERTICES, `a triangle names vertex ${v}`)
      used[v] = 1
    }
    assert.equal(used.indexOf(0), -1, 'a vertex is in no triangle')
    assert.ok(mesh.vertices?.every(Number.isFinite))
  })

  test(`${hemi}: the sulci are a 0 or a 1 for every vertex, about half of each`, () => {
    const sulci = readMz3(files.sulci)
    assert.equal(sulci.faces, null)
    assert.equal(sulci.vertices, null)
    assert.equal(sulci.values?.length, VERTICES)
    const ones = sulci.values?.filter((v) => v === 1).length ?? 0
    const zeros = sulci.values?.filter((v) => v === 0).length ?? 0
    assert.equal(ones + zeros, VERTICES, 'a value is neither 0 nor 1')
    assert.ok(ones > VERTICES * 0.3 && ones < VERTICES * 0.7, `${ones} of ${VERTICES} vertices in a sulcus`)
  })

  test(`${hemi}: the map is the age clusters' −log10(p), and 0 off them`, () => {
    const map = readMz3(files.map)
    assert.equal(map.faces, null)
    assert.equal(map.vertices, null)
    assert.equal(map.values?.length, VERTICES)
    const values = map.values ?? new Float32Array()
    assert.ok(values.every(Number.isFinite), 'a value is not finite')
    // Inside a cluster every vertex passed the cluster-forming threshold.
    const threshold = -Math.log10(run.model.clusterFormingThreshold)
    const inside = values.filter((v) => v !== 0)
    assert.ok(
      inside.every((v) => v >= threshold - 1e-3),
      `a vertex is on a cluster below −log10(p) = ${threshold}`,
    )
    // fsaverage6's vertices are spread over the surface like fsaverage's, so the share on
    // a cluster is about the same as in run.json. A map of another stack would not be.
    const clustered = run.hemispheres[hemi].clusters
      .filter((c) => c.stack === 'age')
      .reduce((sum, c) => sum + c.nVertices, 0)
    const expected = clustered / FSAVERAGE_VERTICES
    const share = inside.length / VERTICES
    assert.ok(Math.abs(share - expected) < 0.02, `${(share * 100).toFixed(1)}% on a cluster, run.json says ${(expected * 100).toFixed(1)}%`)
  })
}

test('the licence of the data sits beside the files', () => {
  const licence = readFileSync(inPublic('/viewer/LICENCE.txt'), 'utf8')
  assert.match(licence, new RegExp(run.credit.licence.replace(/\./g, '\\.')))
})
