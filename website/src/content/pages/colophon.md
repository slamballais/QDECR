---
title: Colophon
description: 'How qdecr.com is made: what it is built with, its type and colours, and what it does not do with your visit.'
---

## Built with

The site is a static site made with [Astro](https://astro.build/): every page is plain HTML by the time it reaches you, with a little JavaScript for search, the copy buttons and the theme switch. The layout and styles are written by hand, as CSS custom properties, without a framework.

- The guide and the project pages are Markdown, in the same repository as the package, under [`website/`](https://github.com/slamballais/QDECR/tree/master/website).
- The [reference](/reference) is generated from the package's own help pages, and the [changelog](/changelog) from its `NEWS.md`, so neither can drift from the release they describe.
- The list of [publications](/cite/publications) comes from [OpenAlex](https://openalex.org/), refreshed monthly.
- Code is highlighted with [Shiki](https://shiki.style/) when the site is built, and search runs on [Pagefind](https://pagefind.app/), whose index is built at the same time.
- [Netlify](https://www.netlify.com/) serves the site.

## Type and colour

Headings are set in [Newsreader](https://github.com/productiontype/Newsreader), the text in [IBM Plex Sans](https://github.com/IBM/plex) and code in IBM Plex Mono. All three are served from this site.

The colours come from FreeSurfer's overlays, the maps drawn on a brain surface: orange for the hot end, for the logo and for emphasis, and blue-teal for the cold end, for links and buttons. The logo is a surface mesh around a glowing triangle, drawn by a script in the repository.

Every pairing of text and background meets WCAG AA contrast in both the light and the dark theme, and a test holds them to it. The [design system](/design) page shows the palette with its measured contrast, the type and the components.

## Your visit

No analytics, no cookies, and no requests to any other site: fonts, scripts and the search index all come from qdecr.com. Search runs in your browser. If you pick a light or dark theme, that choice is kept in your browser's local storage, and nowhere else.

The site works at any width from 320 pixels up and at 400% zoom, with a keyboard, and without motion when your system asks for less.
