// NiiVue reads Zarr volumes through zarrita, which imports the Blosc, LZ4 and Zstandard
// decoders (1.4 MB of WebAssembly, with C libraries' licences of their own) on demand. The
// home page's viewer only ever opens MZ3 surfaces, so they are never asked for, but the
// bundler still writes them into dist/ and the site would deploy them. This hands zarrita
// a stand-in instead, which says why if it is ever called.

import type { AstroIntegration } from 'astro'

const STAND_IN = '\0qdecr:no-zarr-codecs'

export default function noZarrCodecs(): AstroIntegration {
  return {
    name: 'qdecr:no-zarr-codecs',
    hooks: {
      'astro:config:setup': ({ updateConfig }) => {
        updateConfig({
          vite: {
            plugins: [
              {
                name: 'qdecr:no-zarr-codecs',
                enforce: 'pre',
                resolveId: (id) => (/^numcodecs\/(blosc|lz4|zstd)$/.test(id) ? STAND_IN : undefined),
                load: (id) =>
                  id === STAND_IN
                    ? "export default { fromConfig() { throw new Error('qdecr.com ships no Zarr codecs: the viewer reads MZ3 surfaces only.') } }"
                    : undefined,
              },
            ],
          },
        })
      },
    },
  }
}
