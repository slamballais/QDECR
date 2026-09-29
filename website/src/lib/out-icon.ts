// The small arrow after a text link that leaves the site (plan 3.1: links out are marked
// with an icon, links within the site are not). Decorative: the link text already says
// where it goes. One copy of the markup, for OutIcon.astro and the Markdown plugin
// (external-links.ts) alike.

export const OUT_ICON =
  '<svg class="out-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false">' +
  '<path d="M5.5 3.5h7v7M12.5 3.5l-9 9" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"></path>' +
  '</svg>'
