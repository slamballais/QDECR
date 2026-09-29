// WCAG 2.x relative luminance and contrast ratio, as the spec defines them. Used by the
// palette test, which holds every text-on-background pair the site paints to AA, and by
// the design page, which prints the same ratios next to the swatches.
//
// Only six-digit hex is accepted: the tokens are written that way in tokens.css, and a
// token that is not (a three-digit shorthand, an rgb() with alpha) should fail loudly
// rather than be half-understood.

const HEX = /^#[0-9a-f]{6}$/i

/** One sRGB channel (0 to 255), linearised. */
function channel(value: number): number {
  const c = value / 255
  return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4
}

/** Relative luminance of a colour, from 0 (black) to 1 (white). */
export function luminance(hex: string): number {
  if (!HEX.test(hex)) throw new Error(`expected a six-digit hex colour, got ${hex}`)
  const r = channel(parseInt(hex.slice(1, 3), 16))
  const g = channel(parseInt(hex.slice(3, 5), 16))
  const b = channel(parseInt(hex.slice(5, 7), 16))
  return 0.2126 * r + 0.7152 * g + 0.0722 * b
}

/** Contrast ratio between two colours, from 1 to 21. The order does not matter. */
export function contrast(a: string, b: string): number {
  const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x)
  return (hi! + 0.05) / (lo! + 0.05)
}
