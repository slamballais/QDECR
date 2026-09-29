/**
 * The home page viewer's settings and arithmetic: where the camera sits and how it turns
 * between views, the colours, and the files. Run with `npm test`.
 */

import { test } from 'node:test'
import assert from 'node:assert/strict'
import { FIRST_VIEW, GYRI, HEAT, HEMISPHERES, SULCI, azimuthFor, mirrorAzimuth, turnAzimuth, viewerFiles } from './viewer.ts'

test('a turn starts where the camera is and ends where it is going', () => {
  assert.equal(turnAzimuth(90, 270, 0), 90)
  assert.equal(turnAzimuth(90, 200, 1), 200)
})

test('a turn goes the short way round, across 0 if that is shorter', () => {
  // From 350 to 10 is 20 degrees forward, not 340 back.
  const halfway = turnAzimuth(350, 10, 0.5)
  assert.ok(halfway > 350 || halfway < 10, `went the long way: ${halfway}`)
  assert.equal(turnAzimuth(350, 10, 1), 10)
  // And the other way.
  const back = turnAzimuth(10, 350, 0.5)
  assert.ok(back < 10 || back > 350, `went the long way: ${back}`)
})

test('the azimuth stays within 0 to 360', () => {
  for (const t of [0, 0.25, 0.5, 0.75, 1]) {
    const a = turnAzimuth(350, 10, t)
    assert.ok(a >= 0 && a < 360, `${a} at ${t}`)
  }
})

test('a turn eases in and out: slow at the ends, fastest in the middle', () => {
  const at = (t: number) => turnAzimuth(0, 100, t)
  assert.ok(at(0.1) - at(0) < at(0.55) - at(0.45))
  assert.ok(at(1) - at(0.9) < at(0.55) - at(0.45))
  assert.equal(at(0.5), 50)
})

test('the viewer opens on the poster: the left hemisphere, seen from the left', () => {
  // NiiVue's azimuth 90 puts the camera at the subject's left.
  assert.deepEqual(FIRST_VIEW, { hemi: 'lh', side: 'lateral' })
  assert.equal(azimuthFor('lh', 'lateral'), 90)
})

test("each hemisphere's two sides are half a turn apart, and mirror the other's", () => {
  for (const hemi of HEMISPHERES) {
    assert.equal(Math.abs(azimuthFor(hemi, 'lateral') - azimuthFor(hemi, 'medial')), 180)
  }
  assert.equal(azimuthFor('lh', 'lateral'), azimuthFor('rh', 'medial'))
  assert.equal(azimuthFor('lh', 'medial'), azimuthFor('rh', 'lateral'))
})

test("the colours are Freeview's heat, the poster's: red to halfway up the scale, then to yellow", () => {
  const { R, G, B, A, I } = HEAT
  for (const channel of [R, G, B, A]) assert.equal(channel.length, I.length)
  const at = (index: number) => {
    const i = I.indexOf(index)
    assert.notEqual(i, -1, `no stop at ${index}`)
    return [R[i], G[i], B[i]]
  }
  // NiiVue spreads the index over cal_min to cal_max; Freeview, given only those two,
  // puts its midpoint halfway between them.
  assert.equal(I[0], 0)
  assert.equal(I[I.length - 1], 255)
  assert.deepEqual(at(0), [255, 0, 0])
  assert.deepEqual(at(128), [255, 0, 0])
  assert.deepEqual(at(255), [255, 255, 0])
  assert.ok(A.every((a) => a === 255), 'the map is opaque, as linearopaque draws it')
})

test('the files for a hemisphere are its own', () => {
  for (const hemi of HEMISPHERES) {
    for (const url of Object.values(viewerFiles(hemi))) {
      assert.match(url, new RegExp(String.raw`^/viewer/${hemi}\.[a-z.]+\.mz3$`))
    }
  }
})

test('switching hemispheres mirrors the camera, so the reader sees the same side of the other', () => {
  for (const side of ['lateral', 'medial'] as const) {
    assert.equal(mirrorAzimuth(azimuthFor('lh', side)), azimuthFor('rh', side))
    assert.equal(mirrorAzimuth(azimuthFor('rh', side)), azimuthFor('lh', side))
  }
  // From the front or the back, the view stays where it is.
  assert.equal(mirrorAzimuth(0), 0)
  assert.equal(mirrorAzimuth(180), 180)
  // Anything between turns the other way, and stays within 0 to 360.
  assert.equal(mirrorAzimuth(60), 300)
  assert.equal(mirrorAzimuth(300), 60)
})

test('the sulci are the darker grey, as Freeview draws them, and opaque', () => {
  const sulcus = SULCI.R[0] ?? 0
  assert.ok(sulcus < GYRI[0], 'a sulcus is no darker than a gyrus')
  for (const channel of [SULCI.R, SULCI.G, SULCI.B]) assert.ok(channel.every((v) => v === sulcus), 'the sulci are not one grey')
  assert.ok(SULCI.A.every((a) => a === 255))
  assert.equal(GYRI[3], 255)
})
