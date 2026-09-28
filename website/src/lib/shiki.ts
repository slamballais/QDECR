// Syntax highlighting for every code block on the site: Markdown fences (astro.config.mjs)
// and the CodeBlock component alike.
//
// Shiki colours code by writing style="color:..." on every token, which the CSP forbids.
// So the theme below gives each kind of token a placeholder colour that stands for a
// class, and the transformer swaps each placeholder for its tok-* class and drops every
// style. The real colours live in tokens.css (--syntax-*), are painted by global.css, and
// are held to AA on the code panel by the palette test.

import type { ShikiTransformer, ThemeRegistration } from 'shiki'

/** Placeholder colour → the class global.css paints. Plain text gets no class. */
const CLASSES: Record<string, string> = {
  '#000001': 'tok-keyword',
  '#000002': 'tok-string',
  '#000003': 'tok-number',
  '#000004': 'tok-function',
  '#000005': 'tok-comment',
}
const PLAIN = '#000000'

const colourOf = (className: string) => Object.keys(CLASSES).find((key) => CLASSES[key] === className)!

/** TextMate scopes, grouped the way the logo's colours are: heat for operators and
 * values, the lattice's teal and blue for strings and functions. */
export const syntaxTheme: ThemeRegistration = {
  name: 'qdecr',
  type: 'dark',
  fg: PLAIN,
  bg: PLAIN,
  settings: [
    { settings: { foreground: PLAIN, background: PLAIN } },
    { scope: ['comment', 'punctuation.definition.comment'], settings: { foreground: colourOf('tok-comment') } },
    { scope: ['keyword', 'storage'], settings: { foreground: colourOf('tok-keyword') } },
    { scope: ['string', 'punctuation.definition.string'], settings: { foreground: colourOf('tok-string') } },
    { scope: ['constant.numeric', 'constant.language'], settings: { foreground: colourOf('tok-number') } },
    {
      scope: ['entity.name.function', 'support.function', 'entity.name.command'],
      settings: { foreground: colourOf('tok-function') },
    },
  ],
}

/** How a fence's language is labelled in the block's corner; '' for no label. */
const LABELS: Record<string, string> = {
  r: 'R',
  bash: 'Shell',
  sh: 'Shell',
  shell: 'Shell',
  console: 'Shell',
  text: '',
  plaintext: '',
  txt: '',
}

/** Text inside a hast node, however deep. */
type Node = { type: string; value?: string; children?: Node[] }
const textOf = (node: Node): string =>
  node.type === 'text' ? (node.value ?? '') : (node.children ?? []).map(textOf).join('')

/**
 * Turns Shiki's output into the site's code block: classes instead of colours, R output
 * lines (#>) marked, the language label and a copy button (hidden until the script in
 * src/scripts/code-copy.ts wires it up, so a reader without JavaScript never sees a dead
 * control) around a keyboard-focusable <pre>.
 */
export function codeBlock(): ShikiTransformer {
  return {
    name: 'qdecr:code-block',
    preprocess(_code, options) {
      // Keep spaces out of the tokens beside them, so a function's class covers its name
      // and nothing else.
      options.mergeWhitespaces = false
    },
    pre(node) {
      node.properties = { tabindex: '0' }
    },
    code(node) {
      node.properties = {}
    },
    line(node) {
      if (textOf(node).startsWith('#>')) this.addClassToHast(node, 'output')
    },
    span(node) {
      const style = String(node.properties['style'] ?? '')
      delete node.properties['style']
      const colour = /color:\s*(#[0-9a-fA-F]{6})/.exec(style)?.[1]?.toLowerCase()
      const className = colour ? CLASSES[colour] : undefined
      // The R grammar calls a comma an operator; painted as one, every argument list
      // would be spotted with orange.
      if (className && !(className === 'tok-keyword' && /^[,;]$/.test(textOf(node)))) {
        this.addClassToHast(node, className)
      }
    },
    root(root) {
      const lang = this.options.lang
      const label = LABELS[lang] ?? lang
      root.children = [
        {
          type: 'element',
          tagName: 'div',
          properties: { class: 'code-block' },
          children: [
            ...(label
              ? [{ type: 'element' as const, tagName: 'span', properties: { class: 'code-lang' }, children: [{ type: 'text' as const, value: label }] }]
              : []),
            {
              type: 'element',
              tagName: 'button',
              properties: { type: 'button', class: 'code-copy', hidden: true },
              children: [{ type: 'text', value: 'Copy' }],
            },
            ...root.children.filter((child) => child.type === 'element'),
          ],
        },
      ]
    },
  }
}
