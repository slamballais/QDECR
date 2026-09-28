// @ts-check
import { defineConfig } from 'astro/config'
import cspGuard from './src/integrations/csp-guard.ts'
import { codeBlock, syntaxTheme } from './src/lib/shiki.ts'

// Fully static: no adapter and no server runtime. Netlify serves dist/ as it is.
export default defineConfig({
  integrations: [
    // Fails the build on any inline script or style the CSP in netlify.toml would block.
    cspGuard(),
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
