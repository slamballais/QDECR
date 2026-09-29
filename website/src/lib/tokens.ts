// Reads colour tokens out of a stylesheet, so the palette test and the design page check
// the values the browser actually gets rather than a copy of them that could drift.

/**
 * The hex colour custom properties (`--name: #rrggbb;`) declared directly in the block
 * that opens with `selector {`, lower-cased. Other declarations (fonts, sizes, colours
 * with alpha) are skipped.
 *
 * Braces are counted rather than matched with a regex: a block sits inside a media
 * query, and a lazy `\{([^}]*)\}` stops at whichever closing brace comes first, which is
 * wrong in a way that still parses.
 */
export function readTokens(css: string, selector: string): Record<string, string> {
  // The selector must start a line (after indentation) and be followed by its brace, so
  // ':root' does not match the start of ':root[data-theme=...]'.
  const escaped = selector.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')
  const opening = new RegExp(`^[ \\t]*${escaped} \\{`, 'm').exec(css)
  if (!opening) throw new Error(`no ${selector} block in the stylesheet`)

  const start = opening.index + opening[0].length
  let depth = 1
  let end = start
  for (; end < css.length && depth > 0; end++) {
    if (css[end] === '{') depth += 1
    else if (css[end] === '}') depth -= 1
  }

  const tokens: Record<string, string> = {}
  for (const match of css.slice(start, end).matchAll(/^\s*(--[a-z0-9-]+):\s*(#[0-9a-f]{6})\s*;/gim)) {
    tokens[match[1]!] = match[2]!.toLowerCase()
  }
  return tokens
}
