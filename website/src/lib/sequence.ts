// The guide as one reading sequence: Get started, then the tutorials in their numbered
// order. It sets the order of the sidebar and the previous/next links at the foot of each
// page, so both come from the `order` in each page's front matter.

export interface SequencePage {
  /** The page's id in the docs collection, which is also its path: 'tutorials/plotting'. */
  id: string
  title: string
  /** Its place in the sequence; unique. Get started is 0, tutorial n is n. */
  order: number
}

/** The pages in reading order. Two pages with the same order are refused. */
export function readingOrder<T extends SequencePage>(pages: readonly T[]): T[] {
  const sorted = [...pages].sort((a, b) => a.order - b.order)
  for (let i = 1; i < sorted.length; i++) {
    const [before, page] = [sorted[i - 1]!, sorted[i]!]
    if (before.order === page.order) {
      throw new Error(`${before.id} and ${page.id} both have order ${page.order}`)
    }
  }
  return sorted
}

/** The pages before and after `id` in reading order. */
export function neighbours<T extends SequencePage>(pages: readonly T[], id: string): { prev?: T; next?: T } {
  const sorted = readingOrder(pages)
  const at = sorted.findIndex((page) => page.id === id)
  if (at === -1) throw new Error(`${id} is not in the reading sequence`)
  const result: { prev?: T; next?: T } = {}
  if (at > 0) result.prev = sorted[at - 1]!
  if (at < sorted.length - 1) result.next = sorted[at + 1]!
  return result
}
