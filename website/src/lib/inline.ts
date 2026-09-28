// Finds the markup the Content-Security-Policy would block. netlify.toml allows scripts
// and styles from this origin only, never 'unsafe-inline', so an inline script, a <style>
// element, a style attribute or an onclick handler on a page is dead on arrival in
// production, and silently: the browser only says so in the console. Astro and its
// integrations add these on their own (Shiki's code colours, small scripts it inlines),
// so the build checks every page it writes (src/integrations/csp-guard.ts) rather than
// trusting the source.

/** One opening tag, with its attributes, however they are quoted. */
const TAG = /<([a-zA-Z][\w-]*)((?:\s+[^\s"'>/=]+(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s"'=<>`]+))?)*)\s*\/?>/g
const ATTRIBUTE = /([^\s"'>/=]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'=<>`]+)))?/g

function attributes(source: string): Map<string, string> {
  const found = new Map<string, string>()
  for (const match of source.matchAll(ATTRIBUTE)) {
    found.set(match[1]!.toLowerCase(), match[2] ?? match[3] ?? match[4] ?? '')
  }
  return found
}

/**
 * What in `html` a strict CSP would block, one line per finding, or an empty list. A
 * script with a src is fine, and so is a JSON-LD block, which is data the browser never
 * runs.
 */
export function findInline(html: string): string[] {
  const problems: string[] = []
  for (const match of html.matchAll(TAG)) {
    const tag = match[1]!.toLowerCase()
    const attrs = attributes(match[2] ?? '')
    if (tag === 'script' && !attrs.has('src') && attrs.get('type') !== 'application/ld+json') {
      problems.push(`an inline <script>: ${match[0]}`)
    }
    if (tag === 'style') problems.push('a <style> element')
    for (const name of attrs.keys()) {
      if (name === 'style') problems.push(`a style attribute on <${tag}>: ${match[0].slice(0, 120)}`)
      if (name.startsWith('on')) problems.push(`an ${name} handler on <${tag}>`)
    }
  }
  return problems
}
