// The citing works as the pages list them: citations.json less the ones hidden by hand.
// Kept apart from citations.ts, which the tests and the fetch script import without the
// data (like reference-data.ts beside reference.ts).

import citations from '../data/citations.json'
import { hiddenCitations } from '../data/citations-hidden'
import { hideCitations } from './citations'

/** Every work that cites the paper and is shown, newest first. */
export const citingWorks = hideCitations(citations.works, hiddenCitations)
