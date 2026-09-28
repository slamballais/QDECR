// @ts-check
import { defineConfig } from 'astro/config'
import cspGuard from './src/integrations/csp-guard.ts'

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
})
