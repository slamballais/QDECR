// The part of the home page viewer (src/components/Viewer.astro) that every visit loads:
// it shows "Explore in 3D" on the poster, and only when that is pressed imports the viewer
// itself, NiiVue and all, which Vite splits into a file of its own. The first page stays
// small, and nobody downloads a megabyte of surfaces they did not ask for.

import type { Hemisphere, Side } from '../lib/viewer'

for (const root of document.querySelectorAll<HTMLElement>('[data-viewer]')) {
  const open = root.querySelector<HTMLButtonElement>('.viewer-open')
  const firstCanvas = root.querySelector<HTMLCanvasElement>('.viewer-canvas')
  const controls = root.querySelector<HTMLElement>('.viewer-controls')
  const status = root.querySelector<HTMLElement>('.viewer-status')
  if (!open || !firstCanvas || !controls || !status) continue
  let canvas = firstCanvas

  // NiiVue needs WebGL 2. Without it the poster is the whole story, so no button. Asking
  // whether the browser knows WebGL 2 at all, rather than opening a context to find out,
  // keeps the GPU out of a visit that never opens the viewer; a GPU that then refuses is
  // caught below like any other failure.
  if (!('WebGL2RenderingContext' in window)) continue
  open.hidden = false

  const scale = { from: Number(root.dataset['scaleFrom']), to: Number(root.dataset['scaleTo']) }
  const hemiButtons = [...controls.querySelectorAll<HTMLButtonElement>('[data-hemi]')]
  const sideButtons = [...controls.querySelectorAll<HTMLButtonElement>('[data-side]')]
  const hemiNames: Record<Hemisphere, string> = { lh: 'left', rh: 'right' }
  // The caption's "left", which becomes "right" with the other hemisphere.
  const captionHemis = root.closest('figure')?.querySelectorAll<HTMLElement>('[data-viewer-hemi]') ?? []

  open.addEventListener('click', async (event) => {
    open.disabled = true
    open.setAttribute('aria-busy', 'true')
    open.textContent = 'Loading…'
    status.textContent = ''
    try {
      const { mount } = await import('./viewer-3d')
      canvas.hidden = false
      const viewer = await mount(canvas, scale)
      const live = canvas
      root.classList.add('is-live')
      open.remove()
      controls.hidden = false
      // The button is gone, so focus moves on. From the keyboard (a click with no mouse
      // behind it has no detail), into the viewer, where the arrow keys turn it. From the
      // mouse, to the first control: a focused canvas would take the next scroll of the
      // wheel over it as a zoom rather than scroll the page.
      if (event.detail === 0) live.focus()
      else hemiButtons[0]?.focus({ preventScroll: true })

      // One hemisphere at a time: a second press while one loads would be overtaken by
      // the first to arrive, and leave the buttons naming the other.
      let switching = false

      for (const button of hemiButtons) {
        button.addEventListener('click', async () => {
          const hemi = button.dataset['hemi'] as Hemisphere
          if (switching || button.getAttribute('aria-pressed') === 'true') return
          switching = true
          controls.setAttribute('aria-busy', 'true')
          const press = (pressed: HTMLButtonElement) => {
            for (const other of hemiButtons) other.setAttribute('aria-pressed', String(other === pressed))
          }
          const previous = hemiButtons.find((other) => other.getAttribute('aria-pressed') === 'true')
          press(button)
          status.textContent = `Loading the ${hemiNames[hemi]} hemisphere…`
          try {
            await viewer.show(hemi)
            live.setAttribute('aria-label', `The ${hemiNames[hemi]} hemisphere in 3D. Drag, or use the arrow keys, to turn it.`)
            for (const name of captionHemis) name.textContent = hemiNames[hemi]
            status.textContent = ''
          } catch (error) {
            console.error(error)
            if (previous) press(previous)
            status.textContent = `The ${hemiNames[hemi]} hemisphere could not be loaded. Try again.`
          } finally {
            switching = false
            controls.removeAttribute('aria-busy')
          }
        })
      }
      for (const button of sideButtons) {
        button.addEventListener('click', () => viewer.turnTo(button.dataset['side'] as Side))
      }
    } catch (error) {
      // A failed download or a GPU that refuses: the poster stays, and says why.
      console.error(error)
      // A fresh canvas for the next try (mount has let go of this one).
      const fresh = canvas.cloneNode() as HTMLCanvasElement
      fresh.hidden = true
      canvas.replaceWith(fresh)
      canvas = fresh
      open.disabled = false
      open.removeAttribute('aria-busy')
      open.textContent = 'Explore in 3D'
      status.textContent = 'The 3D view could not be loaded. Try again, or keep to the picture.'
    }
  })
}
