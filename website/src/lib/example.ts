// The example analysis: the results of the one real run every output on the site comes
// from, as tools/example/export.R writes them to src/data/example/run.json.
// The shape is checked when the JSON is read, strictly (a key the schema does not know
// is an error too), so that a re-export from a changed script fails the build here
// rather than leave a page showing the wrong numbers. Kept apart from example-data.ts,
// which imports the JSON (like reference.ts beside reference-data.ts). The names follow
// the glossary (src/data/glossary.md): smoothness, cluster-forming threshold,
// cluster-wise p-value, cluster map.

import { z } from 'astro/zod'

/** The example's model: cortical thickness on age and sex, on controls only. */
export const EXAMPLE_FORMULA = 'qdecr_thickness ~ age + sex'

/** ABIDE I is shared under this licence, which every derived figure carries. */
export const EXAMPLE_LICENCE = 'CC BY-NC-SA 3.0'

const stack = z.strictObject({
  /** Its number in stacks(out), which names the files: stack2.coef.mgh. */
  number: z.number().int().min(1),
  /** The column of the design matrix, as R names it: "(Intercept)", "age", "sexmale". */
  name: z.string(),
})

const region = z.strictObject({
  /** A region of the Desikan-Killiany atlas (aparc), as FreeSurfer names it. */
  name: z.string(),
  /** How much of the cluster lies in the region, in percent. */
  ofCluster: z.number().min(0).max(100),
  /** How much of the region the cluster covers, in percent. */
  ofRegion: z.number().min(0).max(100),
})

const cluster = z.strictObject({
  /** The stack the cluster belongs to, by name. */
  stack: z.string(),
  /** Its number within the stack, from 1, as mri_surfcluster numbers them. */
  cluster: z.number().int().min(1),
  nVertices: z.number().int().min(1),
  /** Its area on the white surface. */
  sizeMm2: z.number().min(0),
  /** The cluster-wise p-value, from FreeSurfer's simulations. */
  clusterwiseP: z.number().min(0).max(1),
  /**
   * The vertex with the strongest signal: the −log10(p) there, its number, and its
   * region. The value is null where p underflowed to zero, as it does for the intercept.
   */
  peak: z.strictObject({ value: z.number().nullable(), vertex: z.number().int().min(0), region: z.string() }),
  /** The mean of the measure over the cluster, and of the model's coefficient and its SE. */
  meanThickness: z.number(),
  meanCoefficient: z.number(),
  meanSe: z.number(),
  /** The regions the cluster covers most, largest share first. */
  regions: z.array(region),
})

const hemisphere = z
  .strictObject({
    /** The project's full name, which is its output directory: lh.age_sex.thickness. */
    project: z.string(),
    vertices: z.strictObject({
      /** Vertices per hemisphere of fsaverage. */
      loaded: z.number().int().min(1),
      /** Those inside the mask, where the model was fitted. */
      analysed: z.number().int().min(1),
    }),
    /** The smoothness of the residuals, as a FWHM in mm, which picks the simulation. */
    smoothness: z.number().min(1).max(30),
    /** How long the analysis took, from the call to the result. */
    seconds: z.number().min(0),
    stacks: z.array(stack).min(1),
    clusters: z.array(cluster),
  })
  .superRefine((hemi, ctx) => {
    // mri_surfcluster numbers a stack's clusters 1, 2, 3 from the largest; the export
    // keeps that order, and the pages count on it.
    const seen = new Map<string, number>()
    for (const c of hemi.clusters) {
      const expected = (seen.get(c.stack) ?? 0) + 1
      if (c.cluster !== expected) {
        ctx.addIssue({
          code: 'custom',
          message: `clusters of ${c.stack} are not numbered from 1 in order: got ${c.cluster}, expected ${expected}`,
          path: ['clusters'],
        })
      }
      seen.set(c.stack, expected)
    }
  })

export const exampleRunSchema = z.strictObject({
  /** The day the analysis ran, YYYY-MM-DD. */
  date: z.string().regex(/^\d{4}-\d{2}-\d{2}$/),
  dataset: z
    .strictObject({
      name: z.string(),
      /** The one ABIDE site the controls come from. */
      site: z.string(),
      /** Subjects in the analysis. */
      n: z.number().int().min(1),
      sex: z.strictObject({ female: z.number().int().min(0), male: z.number().int().min(0) }),
      /** Age at scan, in years. */
      age: z.strictObject({ min: z.number(), max: z.number(), mean: z.number(), median: z.number() }),
      /** Subjects left out, with the reason, so a page can say so. */
      excluded: z.array(z.strictObject({ id: z.string(), reason: z.string() })),
    })
    .refine((d) => d.sex.female + d.sex.male === d.n, {
      message: 'the sexes do not add up to n',
      path: ['sex'],
    }),
  software: z.strictObject({
    qdecr: z.string(),
    r: z.string(),
    /** FreeSurfer's build stamp, which names the release and the build. */
    freesurfer: z.string(),
    os: z.string(),
    platform: z.string(),
  }),
  model: z.strictObject({
    formula: z.literal(EXAMPLE_FORMULA, { message: `the formula must be ${EXAMPLE_FORMULA}` }),
    measure: z.string(),
    /** The smoothing of the maps read, as a FWHM in mm. */
    fwhm: z.number().min(0),
    /** The cluster-forming threshold, a vertex-wise p-value: mcz_thr, 0.001 by default. */
    clusterFormingThreshold: z.number().min(0).max(1),
    /** The cluster-wise threshold: cwp_thr, 0.025 by default. */
    clusterwiseThreshold: z.number().min(0).max(1),
    nCores: z.number().int().min(1),
  }),
  /** Both hemispheres: the analysis is whole-brain, and cwp_thr splits 0.05 over the two. */
  hemispheres: z.strictObject({ lh: hemisphere, rh: hemisphere }),
  /**
   * What every page showing the run has to print: the data's licence,
   * which is not the site's, the two projects that collected and preprocessed them with
   * the papers they ask to be cited, and the funding ABIDE asks to be acknowledged.
   */
  credit: z.strictObject({
    licence: z.literal(EXAMPLE_LICENCE, { message: `ABIDE is shared under ${EXAMPLE_LICENCE}` }),
    abide: z.strictObject({ url: z.url(), cite: z.string() }),
    pcp: z.strictObject({ url: z.url(), cite: z.string() }),
    funding: z.string(),
  }),
})

export type ExampleRun = z.infer<typeof exampleRunSchema>
export type ExampleHemisphere = ExampleRun['hemispheres']['lh']
export type ExampleCluster = ExampleHemisphere['clusters'][number]

/** Reads a run, or throws with the first thing wrong with it. */
export function parseExampleRun(data: unknown): ExampleRun {
  const result = exampleRunSchema.safeParse(data)
  if (result.success) return result.data
  const issue = result.error.issues[0]!
  throw new Error(`run.json is not a valid example run: ${issue.path.join('.')}: ${issue.message}`)
}
