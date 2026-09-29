// A stylesheet for the design page's swatches, built from the tokens: one rule per colour
// token that paints an element marked data-fg="--name" in it, and one for data-bg. The
// page cannot set these colours in style attributes, which the CSP forbids, and a rule
// written by hand per token would fall behind the day a token is added.

import type { APIRoute } from 'astro'
import tokensCss from '../styles/tokens.css?raw'
import { themeTokens } from '../lib/palette'

export const GET: APIRoute = () => {
  const names = Object.keys(themeTokens(tokensCss).light)
  const css = names
    .map((name) => `[data-fg='${name}']{color:var(${name})}[data-bg='${name}']{background:var(${name})}`)
    .join('\n')
  return new Response(`${css}\n`, { headers: { 'Content-Type': 'text/css; charset=utf-8' } })
}
