// @ts-check
import { defineConfig } from 'astro/config'
import cspGuard from './src/integrations/csp-guard.ts'
import ogImages from './src/integrations/og-images.ts'
import sitemap from './src/integrations/sitemap.ts'
import { satteri } from '@astrojs/markdown-satteri'
import { callouts } from './src/lib/callouts.ts'
import { tables } from './src/lib/tables.ts'
import { figures } from './src/lib/figures.ts'
import { exampleOutput } from './src/lib/example-output.ts'
import { readFileSync } from 'node:fs'
import { codeBlock, syntaxTheme } from './src/lib/shiki.ts'

// What R printed in the example analysis (tools/example/export.R writes it).
const exampleOutputDir = new URL('./src/data/example/output/', import.meta.url)

// Fully static: no adapter and no server runtime. Netlify serves dist/ as it is.
export default defineConfig({
  integrations: [
    // Fails the build on any inline script or style the CSP in netlify.toml would block.
    cspGuard(),
    // Draws each page's share card from its title (tools/og.mjs).
    ogImages(),
    // Lists every indexable page in sitemap.xml, which robots.txt points to.
    sitemap(),
  ],
  // The canonical host is the apex domain.
  site: 'https://qdecr.com',
  build: {
    // Emit /cite.html rather than /cite/index.html, so URLs have no trailing slash.
    format: 'file',
    // Astro inlines small stylesheets into <style> tags by default. The CSP (netlify.toml)
    // allows styles from this origin only and never 'unsafe-inline',
    // so every stylesheet is a file.
    inlineStylesheets: 'never',
  },
  trailingSlash: 'never',
  markdown: {
    // Sätteri is Astro's default Markdown processor; it is named here only to add four
    // plugins: GitHub's "> [!NOTE]" blockquotes rendered as the site's note, tip and warning
    // boxes (src/lib/callouts.ts); tables wrapped so a wide one scrolls on its own
    // (src/lib/tables.ts); a lone image set as a figure with its caption and, for the
    // example data, its credit (src/lib/figures.ts); and code blocks marked output=<file>
    // filled with what R printed in the example analysis (src/lib/example-output.ts).
    processor: satteri({
      mdastPlugins: [
        callouts(),
        tables(),
        figures(),
        exampleOutput({ read: (name) => readFileSync(new URL(name, exampleOutputDir), 'utf8') }),
      ],
    }),
    // Code in Markdown is highlighted into the site's tok-* classes rather than Shiki's
    // inline colours, which the CSP forbids (src/lib/shiki.ts).
    shikiConfig: { theme: syntaxTheme, transformers: [codeBlock()] },
  },
  vite: {
    build: {
      // Vite inlines small scripts and assets into the page. An inline script is blocked
      // by the CSP, so nothing is inlined; the build check (csp-guard) would catch it.
      assetsInlineLimit: 0,
    },
  },
})
