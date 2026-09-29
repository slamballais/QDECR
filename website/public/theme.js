// Theme override. The stylesheet follows prefers-color-scheme on its own; this script only
// applies a choice the reader has made with the toggle, and it is loaded synchronously in
// <head> so that choice is in place before the first paint. With nothing stored it does no
// more than wire up the buttons. Once a choice is stored the page stops following the
// system until the toggle is pressed again.
//
// It is a plain file in public/ rather than a script Astro bundles, because a bundled
// script is deferred, and an inline one would need 'unsafe-inline' in the CSP.
;(() => {
  const KEY = 'theme'
  // The paper token of each scheme, for the theme-color metas. A copy of tokens.css, which
  // cannot be read this early; src/lib/palette.test.ts checks the two agree.
  const PAPER = { light: '#f7f5f0', dark: '#0d0f12' }
  const root = document.documentElement

  const stored = () => {
    try {
      const value = localStorage.getItem(KEY)
      return value === 'light' || value === 'dark' ? value : null
    } catch {
      return null
    }
  }

  const effective = () =>
    root.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light')

  // The button names the theme a press switches to, in words for screen readers and as a
  // tooltip, and shows it as a moon on paper and a sun in the dark (the CSS does that).
  const label = () => {
    const next = effective() === 'dark' ? 'light' : 'dark'
    for (const button of document.querySelectorAll('.theme-toggle')) {
      button.setAttribute('aria-label', `Switch to the ${next} theme`)
      button.title = `Switch to the ${next} theme`
    }
  }

  const apply = (theme) => {
    root.dataset.theme = theme
    // 'only light' also opts a chosen light theme out of Chrome's auto-dark on Android.
    root.style.colorScheme = theme === 'light' ? 'only light' : 'dark'
    // The media-scoped theme-color metas would otherwise follow the system, not the choice.
    for (const meta of document.querySelectorAll('meta[name="theme-color"]')) {
      meta.setAttribute('content', PAPER[theme])
    }
    label()
  }

  const choice = stored()
  if (choice) apply(choice)

  document.addEventListener('DOMContentLoaded', () => {
    label()
    for (const button of document.querySelectorAll('.theme-toggle')) {
      button.addEventListener('click', () => {
        const next = effective() === 'dark' ? 'light' : 'dark'
        try {
          localStorage.setItem(KEY, next)
        } catch {
          // Storage blocked: the switch still works for this page view.
        }
        apply(next)
      })
    }
    // While nothing is stored the page follows the system; keep the labels in step with it.
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
      if (!stored()) label()
    })
  })
})()
