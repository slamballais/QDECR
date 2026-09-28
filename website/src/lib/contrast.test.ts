import { test } from 'node:test'
import assert from 'node:assert/strict'
import { contrast, luminance } from './contrast.ts'

test('black and white are the extremes of relative luminance', () => {
  assert.equal(luminance('#000000'), 0)
  assert.equal(luminance('#ffffff'), 1)
})

test('black on white is 21:1, whichever way round', () => {
  assert.equal(contrast('#000000', '#ffffff'), 21)
  assert.equal(contrast('#ffffff', '#000000'), 21)
})

test('a colour on itself is 1:1', () => {
  assert.equal(contrast('#2d6cc4', '#2d6cc4'), 1)
})

test('#767676 on white is the familiar 4.54:1, just over AA', () => {
  assert.equal(contrast('#767676', '#ffffff').toFixed(2), '4.54')
})

test('hex case does not matter', () => {
  assert.equal(contrast('#FF7A1C', '#020305'), contrast('#ff7a1c', '#020305'))
})

test('anything but a six-digit hex colour is refused', () => {
  assert.throws(() => luminance('#fff'), /six-digit/)
  assert.throws(() => luminance('rgb(0 0 0)'), /six-digit/)
})
