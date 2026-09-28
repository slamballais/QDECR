// Which colour is painted on which, and how much contrast each pairing needs. The palette
// test holds src/styles/tokens.css to this list in both themes, and the design page prints
// the same list with the measured ratios, so the two cannot disagree.
//
// A pair belongs here the day something paints it. A token that is never set on a
// background it has to be legible against is not listed.

import { readTokens } from './tokens.ts'

/** The selectors of the token blocks in tokens.css. */
export const THEME_BLOCKS = {
  light: ':root',
  // The dark tokens are written twice, because CSS cannot say "the system is dark and no
  // light choice is stored, or a dark choice is stored" in one selector. The test checks
  // that the two copies match.
  dark: ":root[data-theme='dark']",
  darkBySystem: ":root:not([data-theme='light'])",
} as const

/** Both themes' colour tokens, read out of tokens.css. */
export function themeTokens(css: string): Record<'light' | 'dark', Record<string, string>> {
  return { light: readTokens(css, THEME_BLOCKS.light), dark: readTokens(css, THEME_BLOCKS.dark) }
}

export interface Pair {
  fg: string
  bg: string
  /**
   * WCAG AA: 4.5 for text, 3 for large text and for the parts of a control a reader has
   * to see to use it, such as a focus ring (1.4.11, non-text contrast).
   */
  min: 4.5 | 3
  /** What paints it this way. */
  where: string
}

export const PAIRS: Pair[] = [
  { fg: '--ink', bg: '--paper', min: 4.5, where: 'body text' },
  { fg: '--ink', bg: '--surface', min: 4.5, where: 'text on a raised panel' },
  { fg: '--ink', bg: '--sunken', min: 4.5, where: 'inline code' },
  { fg: '--muted', bg: '--paper', min: 4.5, where: 'captions, the footer, secondary text' },
  { fg: '--muted', bg: '--surface', min: 4.5, where: 'secondary text on a raised panel' },
  { fg: '--muted', bg: '--sunken', min: 4.5, where: 'secondary text on a sunken panel' },
  { fg: '--accent', bg: '--paper', min: 4.5, where: 'links and secondary buttons' },
  { fg: '--accent', bg: '--surface', min: 4.5, where: 'links on a raised panel' },
  { fg: '--accent', bg: '--sunken', min: 4.5, where: 'links on a sunken panel' },
  { fg: '--accent-strong', bg: '--paper', min: 4.5, where: 'a link on hover' },
  { fg: '--accent-strong', bg: '--accent-soft', min: 4.5, where: 'a secondary button on hover' },
  { fg: '--accent', bg: '--accent-soft', min: 4.5, where: 'the current page in navigation, links in a note' },
  { fg: '--on-accent', bg: '--accent', min: 4.5, where: 'primary button text' },
  { fg: '--on-accent', bg: '--accent-strong', min: 4.5, where: 'primary button text on hover' },
  { fg: '--ink', bg: '--accent-soft', min: 4.5, where: 'the body of a note' },
  { fg: '--tip', bg: '--tip-bg', min: 4.5, where: 'the label of a tip' },
  { fg: '--ink', bg: '--tip-bg', min: 4.5, where: 'the body of a tip' },
  { fg: '--accent', bg: '--tip-bg', min: 4.5, where: 'links in a tip' },
  { fg: '--warn', bg: '--warn-bg', min: 4.5, where: 'the label of a warning' },
  { fg: '--ink', bg: '--warn-bg', min: 4.5, where: 'the body of a warning' },
  { fg: '--accent', bg: '--warn-bg', min: 4.5, where: 'links in a warning' },
  { fg: '--code-ink', bg: '--code-bg', min: 4.5, where: 'code' },
  { fg: '--code-muted', bg: '--code-bg', min: 4.5, where: 'comments, R output, the language label' },
  { fg: '--syntax-keyword', bg: '--code-bg', min: 4.5, where: 'keywords and operators in code' },
  { fg: '--syntax-string', bg: '--code-bg', min: 4.5, where: 'strings in code' },
  { fg: '--syntax-number', bg: '--code-bg', min: 4.5, where: 'numbers and constants in code' },
  { fg: '--syntax-function', bg: '--code-bg', min: 4.5, where: 'function names in code' },
  { fg: '--heat', bg: '--paper', min: 3, where: 'the focus ring' },
  { fg: '--heat', bg: '--surface', min: 3, where: 'the focus ring on a raised panel' },
  { fg: '--heat', bg: '--sunken', min: 3, where: 'the focus ring on a sunken panel' },
  { fg: '--heat', bg: '--code-bg', min: 3, where: 'the focus ring inside a code block' },
]
