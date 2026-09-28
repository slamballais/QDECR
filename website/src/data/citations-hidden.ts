// Works in citations.json that the publications page leaves out, each with its reason.
// citations.json is rewritten from OpenAlex every month, so corrections to it live here.
// The rule: the paper itself, and a preprint whose published version is also listed
// (or, for two copies of one preprint, the second copy). An id that OpenAlex stops
// returning can stay; it hides nothing.

import type { HiddenCitation } from '../lib/citations'

export const hiddenCitations: HiddenCitation[] = [
  { id: 'W3158050040', reason: 'The QDECR paper itself, which OpenAlex counts as citing itself.' },
  { id: 'W4405641072', reason: 'bioRxiv preprint of W4415724738 (Translational Psychiatry, 2025).' },
  { id: 'W4281677564', reason: 'Research Square preprint of W4375954885 (European Journal of Epidemiology, 2023).' },
  { id: 'W4225394107', reason: 'medRxiv preprint of W4308637228 (eLife, 2022).' },
  { id: 'W3161864457', reason: 'Preprints.org preprint of W4384068781 (Translational Psychiatry, 2023).' },
  { id: 'W4285129253', reason: 'SSRN copy of the bioRxiv preprint W4281923042.' },
]
