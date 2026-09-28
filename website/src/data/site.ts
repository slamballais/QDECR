// Facts about the site and the package that more than one page prints. The version is not
// here: it comes from the package's DESCRIPTION (src/lib/description.ts).

export const site = {
  name: 'QDECR',
  url: 'https://qdecr.com',
  description:
    'QDECR is an R package for vertex-wise statistical analysis of FreeSurfer surface data, ' +
    "with FreeSurfer's own correction for multiple testing.",
  repo: 'https://github.com/slamballais/QDECR',
  /** Where the site's own source lives in the repo, for the "Edit this page" link. */
  editBase: 'https://github.com/slamballais/QDECR/edit/master/website/',
  /**
   * The release in development, for the home page's status line. DESCRIPTION only knows
   * the released version. If 0.10.0 is not underway within a few months,
   * set this to null and the line says only which release is current.
   */
  nextRelease: '0.10.0' as string | null,
  /** The paper to cite, as Crossref has it for its DOI. */
  paper: {
    url: 'https://doi.org/10.3389/fninf.2021.561689',
    doi: '10.3389/fninf.2021.561689',
    title: 'QDECR: A Flexible, Extensible Vertex-Wise Analysis Framework in R',
    authors: [
      { given: 'Sander', family: 'Lamballais', orcid: '0000-0003-3118-6330' },
      { given: 'Ryan L.', family: 'Muetzel', orcid: '0000-0003-3215-1287' },
    ],
    journal: 'Frontiers in Neuroinformatics',
    volume: '15',
    /** Frontiers numbers articles rather than pages. */
    article: '561689',
    published: '2021-04-22',
  },
} as const

export interface NavItem {
  label: string
  href: string
}

/** The masthead: the site's own sections, in reading order. */
export const nav: NavItem[] = [
  { label: 'Get started', href: '/get-started' },
  { label: 'Tutorials', href: '/tutorials' },
  { label: 'Reference', href: '/reference' },
  { label: 'Changelog', href: '/changelog' },
  { label: 'Cite', href: '/cite' },
  { label: 'Help', href: '/help' },
]

/** The footer's pages: about the project and the site rather than about using it. */
export const footerNav: NavItem[] = [
  { label: 'About', href: '/about' },
  { label: 'Glossary', href: '/glossary' },
  { label: 'Colophon', href: '/colophon' },
]
