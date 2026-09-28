// The guide (Get started and the tutorials) as the docs layout needs it: the pages in
// reading order and the sidebar that lists them.

import { getCollection, type CollectionEntry } from 'astro:content'
import { readingOrder, type PageLink } from './sequence'
import type { SidebarGroup } from './sidebar'

export interface GuidePage extends PageLink {
  id: string
  order: number
  /** A numbered tutorial, rather than a page like Get started. */
  tutorial: boolean
  /** Its label above the title and on its share card: "Tutorial 7", or "Guide". */
  section: string
  entry: CollectionEntry<'docs'>
}

/** Every page of the guide, in reading order. */
export async function guidePages(): Promise<GuidePage[]> {
  const entries = await getCollection('docs')
  return readingOrder(
    entries.map((entry) => {
      const tutorial = entry.id.startsWith('tutorials/')
      return {
        id: entry.id,
        href: `/${entry.id}`,
        title: entry.data.title,
        order: entry.data.order,
        tutorial,
        section: tutorial ? `Tutorial ${entry.data.order}` : 'Guide',
        entry,
      }
    }),
  )
}

/** The guide's sidebar: where to start, then the numbered tutorials. */
export function guideSidebar(pages: GuidePage[]): SidebarGroup[] {
  const tutorials = pages.filter((page) => page.tutorial)
  return [
    {
      label: 'Start here',
      items: pages.filter((page) => !page.tutorial).map(({ href, title }) => ({ href, label: title })),
    },
    {
      label: 'Tutorials',
      href: '/tutorials',
      numbered: true,
      items: tutorials.map(({ href, title }) => ({ href, label: title })),
    },
  ]
}
