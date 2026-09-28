// The reference as the pages use it: the topics from reference.json in their groups, and
// the sidebar that lists them. Kept apart from reference.ts, which the tests import
// without the JSON.

import data from '../data/reference.json'
import { groupTopics, type Topic } from './reference'
import type { SidebarGroup } from './sidebar'

export const topics = data.topics as Topic[]
export const groups = groupTopics(topics)

export const topicHref = (topic: Pick<Topic, 'slug'>) => `/reference/${topic.slug}`

export const referenceSidebar: SidebarGroup[] = groups.map((group) => ({
  label: group.label,
  items: group.topics.map((topic) => ({ href: topicHref(topic), label: topic.name, code: true })),
}))
