// @ts-check
import { defineConfig } from 'astro/config'

// Fully static: no adapter and no server runtime. Netlify serves dist/ as it is.
export default defineConfig({
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
