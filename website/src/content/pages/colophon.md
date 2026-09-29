---
title: Colophon
description: 'How qdecr.com is made: what it is built with, its type and colours, and what it does not do with your visit.'
---

## Built with

The site is a static site made with [Astro](https://astro.build/): every page is plain HTML by the time it reaches you, with a little JavaScript for search, the copy buttons and the theme switch, and more for the 3D viewer, which loads only when you ask for it. The layout and styles are written by hand, as CSS custom properties, without a framework.

- The guide and the project pages are Markdown, in the same repository as the package, under [`website/`](https://github.com/slamballais/QDECR/tree/master/website).
- The [reference](/reference) is generated from the package's own help pages, and the [changelog](/changelog) from its `NEWS.md`, so neither can drift from the release they describe.
- The list of [publications](/cite/publications) comes from [OpenAlex](https://openalex.org/), refreshed monthly.
- Code is highlighted with [Shiki](https://shiki.style/) when the site is built, and search runs on [Pagefind](https://pagefind.app/), whose index is built at the same time.
- The 3D viewer on the home page is [NiiVue](https://niivue.com/), served from this site with the rest of its code (its licence is below). It draws fsaverage6, a brain surface that comes with [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/), with the same map as the picture it replaces.
- [Netlify](https://www.netlify.com/) serves the site.

## Type and colour

Headings are set in [Newsreader](https://github.com/productiontype/Newsreader), the text in [IBM Plex Sans](https://github.com/IBM/plex) and code in IBM Plex Mono. All three are served from this site.

The colours come from FreeSurfer's overlays, the maps drawn on a brain surface: orange for the hot end, for the logo and for emphasis, and blue-teal for the cold end, for links and buttons. The logo is a surface mesh around a glowing triangle, drawn by a script in the repository.

Every pairing of text and background meets WCAG AA contrast in both the light and the dark theme, and a test holds them to it. The [design system](/design) page shows the palette with its measured contrast, the type and the components.

## The example data

Every result on this site, the printed output, the cluster tables, the histograms, the snapshots and the viewer, comes from one real analysis: cortical thickness against age and sex in 99 typically developing participants of the NYU site of [ABIDE I](https://fcon_1000.projects.nitrc.org/indi/abide/), the Autism Brain Imaging Data Exchange, as preprocessed with FreeSurfer by the [Preprocessed Connectomes Project](http://preprocessed-connectomes-project.org/abide/). The analysis is [reproducible from the repository](https://github.com/slamballais/QDECR/tree/master/website/tools/example).

Those data are shared for non-commercial research under the [Creative Commons Attribution-NonCommercial-ShareAlike 3.0](https://creativecommons.org/licenses/by-nc-sa/3.0/) licence, and so is everything on this site derived from them: the figures, the numbers and the viewer's files carry that licence, whatever the licence of the package or of the site's own text. If you use the data, cite Di Martino et al. (2014, [10.1038/mp.2013.78](https://doi.org/10.1038/mp.2013.78)) for ABIDE and Craddock et al. (2013, [10.3389/conf.fninf.2013.09.00041](https://doi.org/10.3389/conf.fninf.2013.09.00041)) for the preprocessing, as we do.

ABIDE asks that its funding be acknowledged: primary support for the work by Adriana Di Martino was provided by the NIMH (K23MH087770) and the Leon Levy Foundation; primary support for the work by Michael P. Milham and the INDI team was provided by gifts from Joseph P. Healy and the Stavros Niarchos Foundation to the Child Mind Institute, as well as by an NIMH award to MPM (R03MH096321).

## Your visit

No analytics, no cookies, and no requests to any other site: fonts, scripts, the search index and the viewer's surfaces all come from qdecr.com. Search runs in your browser. If you pick a light or dark theme, that choice is kept in your browser's local storage, and nowhere else.

The site works at any width from 320 pixels up and at 400% zoom, with a keyboard, and without motion when your system asks for less.

## NiiVue's licence

The viewer's code includes NiiVue, which is shared under this licence:

> BSD 2-Clause License
>
> Copyright (c) 2021, Niivue
>
> Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
>
> 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
> 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
