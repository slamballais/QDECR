// The guide (Get started and the tutorials) as the docs layout needs it: the pages in
// reading order and the sidebar that lists them.

import { getCollection, type CollectionEntry } from 'astro:content'
import { readingOrder } from './sequence'
import type { SidebarGroup } from './sidebar'

export interface GuidePage {
  id: string
  href: string
  title: string
  order: number
  entry: CollectionEntry<'docs'>
}

/** Every page of the guide, in reading order. */
export async function guidePages(): Promise<GuidePage[]> {
  const entries = await getCollection('docs')
  return readingOrder(
    entries.map((entry) => ({
      id: entry.id,
      href: `/${entry.id}`,
      title: entry.data.title,
      order: entry.data.order,
      entry,
    })),
  )
}

/** The guide's sidebar: where to start, then the numbered tutorials. */
export function guideSidebar(pages: GuidePage[]): SidebarGroup[] {
  const tutorials = pages.filter((page) => page.id.startsWith('tutorials/'))
  return [
    {
      label: 'Start here',
      items: pages.filter((page) => !page.id.includes('/')).map(({ href, title }) => ({ href, label: title })),
    },
    {
      label: 'Tutorials',
      href: '/tutorials',
      numbered: true,
      items: tutorials.map(({ href, title }) => ({ href, label: title })),
    },
  ]
}
