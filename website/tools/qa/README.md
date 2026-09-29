# Checks before a release

Every build already refuses a page with anything inline the CSP would block
(`src/integrations/csp-guard.ts`) and a bundle carrying code whose licence `/licences`
does not print (`src/integrations/licence-guard.ts`), and `npm test` holds the redirects,
headers, contrast pairs and the rest. What is here is slower, or needs a browser or the
network, so it is run by hand before a release.

All of it runs against the built site served the way Netlify serves it, with the headers
of `netlify.toml` (the CSP is enforced), its redirects, its pretty URLs and compression:

```sh
npm run build
npm run serve      # http://localhost:8888, leave it running
```

Then, from `website/`, in another terminal:

| Check | Command | Passes when |
| --- | --- | --- |
| Lighthouse | `node tools/qa/lighthouse.ts` | Every page scores 95 or more in every category (SEO aside on the two noindex pages), nothing is requested from another origin, and the home page stays under 100 KB without its fonts and 140 KB of fonts. Reports land in `tools/qa/reports/`. |
| Reflow, themes and axe | `node tools/qa/browser.ts` | At 320 px, which is also a 1280 px window at 400% zoom, no page scrolls sideways, axe finds nothing and nothing is logged as an error, in the light theme, the dark theme the system asks for, and the dark theme the toggle picks; on the home page, with the 3D viewer open as well. |
| Links | `lychee --config lychee.toml http://localhost:8888/...` (see below) | No link inside the site is broken, anchors included, and none outside it is dead. |
| HTML | `vnu --skip-non-html --errors-only dist` | The Nu HTML Checker, the W3C's, reports no errors. |
| Spelling | `npx cspell@10.3.5 lint --no-progress` | No unknown words in the prose, the reference and the glossary (`cspell.config.yaml` holds the words the dictionaries lack). |

Both scripts take page paths to check only those: `node tools/qa/browser.ts / /licences`.
In Git Bash on Windows, set `MSYS_NO_PATHCONV=1` first, or it turns `/` into a Windows
path.

## Tools that are not npm packages

Lighthouse and cspell come through `npx` at pinned versions, so Netlify does not install
them on every deploy. Two checks need programs of their own:

- **lychee**, the link checker: a single binary from its
  [releases](https://github.com/lycheeverse/lychee/releases). Give it the served pages,
  not the files in `dist/`: from the files it cannot tell `/cite` (the page, `cite.html`)
  from `cite/` (the folder holding `publications.html`), and resolves `#anchors` against
  the wrong page. In Git Bash:

  ```sh
  lychee --config lychee.toml $(cd dist && find . -name '*.html' ! -path './pagefind/*' \
    | sed 's|^\./||; s|index\.html$||; s|\.html$||; s|^|http://localhost:8888/|')
  ```

- **vnu**, the Nu HTML Checker: `vnu.windows.zip`, `vnu.linux.zip` or `vnu.osx.zip` from
  its [releases](https://github.com/validator/validator/releases) carries its own Java.
  (`vnu.jar` alone needs Java 11 or later.)

## What the checks do not settle

- Links to `github.com/slamballais/QDECR/tree/master/website/...` are 404 until the
  site's branch is merged into `master`.
- Publishers (Wiley, JAMA, OUP, MIT Press) and openalex.org answer an automated client
  with 403 however healthy the link; `lychee.toml` accepts 403 for that reason. A DOI that
  no longer resolves is still caught: doi.org itself answers 404.
- The Nu checker warns that `/design` has a second `<h1>`: the type specimen's, which
  shows what a page title looks like. It is left as it is.
