// The example run as the pages use it: src/data/example/run.json, checked against the
// schema as it is imported, so that a bad export fails the build. Kept apart from
// example.ts, which the tests import without the JSON (like reference-data.ts beside
// reference.ts).

import data from '../data/example/run.json'
import { parseExampleRun } from './example'

/** The one analysis every output on the site comes from. */
export const exampleRun = parseExampleRun(data)
