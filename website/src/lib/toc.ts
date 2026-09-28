// The "On this page" list: a page's sections and their subsections, from the headings
// Astro collects while rendering it.

export interface Heading {
  depth: number
  slug: string
  text: string
}

export interface TocItem {
  slug: string
  text: string
  children: TocItem[]
}

/**
 * h2s as the list and h3s nested under the h2 before them. The h1 is the page's title and
 * anything below h3 is too fine to be worth a link.
 */
export function buildToc(headings: readonly Heading[]): TocItem[] {
  const toc: TocItem[] = []
  for (const { depth, slug, text } of headings) {
    const item = { slug, text, children: [] }
    const parent = toc.at(-1)
    if (depth === 2 || (depth === 3 && !parent)) toc.push(item)
    else if (depth === 3) parent!.children.push(item)
  }
  return toc
}
