// The home page's 3D viewer: NiiVue drawing the poster's map, the age
// stack's −log10(p) on its significant clusters, on fsaverage6's inflated surface, in the
// poster's colours and from the poster's side, so the picture comes alive in place.
// src/scripts/viewer.ts imports this only when the reader asks for the viewer, so NiiVue
// and the surfaces cost the first page nothing.
//
// Everything comes from this origin: NiiVue is bundled, the surfaces are in public/viewer/.
// NiiVue draws its font and its lighting from images it carries as data: URLs, which is why
// the CSP allows img-src data: (netlify.toml).

import { Niivue, NVMeshLayerDefaults } from '@niivue/niivue'
import {
  FIRST_VIEW,
  GYRI,
  HEAT,
  SULCI,
  TURN_MS,
  azimuthFor,
  mirrorAzimuth,
  turnAzimuth,
  viewerFiles,
  type Hemisphere,
  type Scale,
  type Side,
} from '../lib/viewer'

export interface Viewer {
  /** Swaps in the other hemisphere, seen from the same side. */
  show(hemi: Hemisphere): Promise<void>
  /** Turns the camera to look at a side of the hemisphere shown. */
  turnTo(side: Side): void
}

// How far one press of an arrow key turns the camera, in degrees.
const KEY_STEP = 15

// How much larger than NiiVue's default fit the surface is drawn, so that it fills the
// frame as it does on the poster, with room left to turn it.
const ZOOM = 1.4

// An opaque layer read from a file. NiiVue's type asks for every field a loaded layer has;
// the rest are its defaults, as the loader would fill them in.
const layer = (url: string, colormap: string, cal_min: number, cal_max: number) => ({
  ...NVMeshLayerDefaults,
  values: [],
  url,
  colormap,
  cal_min,
  cal_max,
  opacity: 1,
})

const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)')

export async function mount(canvas: HTMLCanvasElement, scale: Scale): Promise<Viewer> {
  const nv = new Niivue({
    // The poster's black, which is also the stage's.
    backColor: [0, 0, 0, 1],
    show3Dcrosshair: false,
    isOrientCube: false,
    isColorbar: false,
    loadingText: '',
    // No dropping files onto the canvas to load them, which is NiiVue's default.
    dragAndDropEnabled: false,
    // The wheel zooms only once the canvas has focus, so scrolling the page past the
    // viewer scrolls the page.
    scrollRequiresFocus: true,
    // NiiVue's letter keys switch clip planes and the view mode, which would take the
    // viewer out of 3D; an empty key matches no key.
    clipPlaneHotKey: '',
    cycleClipPlaneHotKey: '',
    viewModeHotKey: '',
  })
  let hemi: Hemisphere = FIRST_VIEW.hemi

  async function load(next: Hemisphere) {
    const files = viewerFiles(next)
    for (const mesh of [...nv.meshes]) nv.removeMesh(mesh)
    await nv.loadMeshes([
      {
        url: files.mesh,
        rgba255: GYRI,
        layers: [
          // 1 in a sulcus: the dark grey from 0.5 up; 0 on a gyrus, below cal_min and so
          // transparent, which leaves the surface's own grey.
          layer(files.sulci, 'qdecr-sulci', 0.5, 1),
          // The map as the poster draws it: 0 off the clusters, and so transparent.
          layer(files.map, 'qdecr-heat', scale.from, scale.to),
        ],
      },
    ])
    // Diffuse light without NiiVue's default highlights, which dull the heat colours;
    // closer to Freeview's.
    nv.setMeshShader(nv.meshes[0].id, 'Matte')
    hemi = next
  }

  function look(azimuth: number, elevation = 0) {
    nv.setRenderAzimuthElevation(azimuth, elevation)
  }

  let turnFrame = 0
  function turnTo(side: Side) {
    cancelAnimationFrame(turnFrame)
    const from = nv.scene.renderAzimuth
    const to = azimuthFor(hemi, side)
    const fromElevation = nv.scene.renderElevation
    if (reducedMotion.matches) {
      look(to)
      return
    }
    const start = performance.now()
    const step = (now: number) => {
      const progress = Math.min(1, (now - start) / TURN_MS)
      look(turnAzimuth(from, to, progress), fromElevation * (1 - progress))
      if (progress < 1) turnFrame = requestAnimationFrame(step)
    }
    turnFrame = requestAnimationFrame(step)
  }

  // The arrow keys turn the surface while the canvas has focus. NiiVue's own arrow keys
  // work only while the mouse is over it, and by a degree a press.
  canvas.addEventListener('keydown', (event) => {
    const turns: Record<string, [number, number]> = {
      ArrowLeft: [-KEY_STEP, 0],
      ArrowRight: [KEY_STEP, 0],
      ArrowUp: [0, KEY_STEP],
      ArrowDown: [0, -KEY_STEP],
    }
    const turn = turns[event.key]
    if (!turn) return
    event.preventDefault()
    event.stopImmediatePropagation()
    cancelAnimationFrame(turnFrame)
    const elevation = Math.max(-90, Math.min(90, nv.scene.renderElevation + turn[1]))
    look((nv.scene.renderAzimuth + turn[0] + 360) % 360, elevation)
  }, { capture: true })

  try {
    await nv.attachToCanvas(canvas)
    nv.setSliceType(nv.sliceTypeRender)
    nv.addColormap('qdecr-heat', HEAT)
    nv.addColormap('qdecr-sulci', SULCI)
    await load(hemi)
  } catch (error) {
    // Let go of the canvas and of the listeners NiiVue added, so that the page can try
    // again on a fresh canvas.
    nv.cleanup()
    throw error
  }
  nv.scene.volScaleMultiplier = ZOOM
  look(azimuthFor(hemi, FIRST_VIEW.side))

  return {
    async show(next) {
      if (next === hemi) return
      cancelAnimationFrame(turnFrame)
      const { renderAzimuth, renderElevation } = nv.scene
      await load(next)
      // From the same angle, mirrored: a lateral view stays lateral.
      look(mirrorAzimuth(renderAzimuth), renderElevation)
    },
    turnTo,
  }
}
