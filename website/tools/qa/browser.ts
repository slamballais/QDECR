// Opens every page the build wrote in headless Chrome, 320 CSS pixels wide, twice: as a
// phone (320 by 800), and as a 1280 by 1024 desktop window at 400% zoom (320 by 256, four
// device pixels to each CSS pixel), which is why WCAG's reflow criterion names that width.
// Each in the light theme, the dark theme the system asks for, and the dark theme picked
// with the toggle (tokens.css writes the dark tokens twice, once for each). On each it:
//
// - checks nothing makes the page scroll sideways;
// - runs axe-core, the accessibility engine Lighthouse uses, which Lighthouse itself only
//   ever runs in the light theme, on a phone;
// - fails on anything the page logs as an error, a CSP violation included, and on any
//   request to another origin.
//
// A page with the 3D viewer then has it opened, as a reader would, and is checked again
// open: Lighthouse only ever sees the poster. A viewer that cannot open fails.
//
//   npm run build && npm run serve        (in one terminal)
//   node tools/qa/browser.ts [page...]    (in another; all pages when none are named)
//
// It drives Chrome over the DevTools protocol with Node's own WebSocket, so it needs no
// browser-automation package. Chrome is found in its usual place, or set CHROME.

import { spawn } from 'node:child_process'
import { existsSync } from 'node:fs'
import { mkdtemp, readFile, rm } from 'node:fs/promises'
import { createRequire } from 'node:module'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { ORIGIN, pagesToCheck } from './pages.ts'

const WIDTH = 320
const VIEWPORTS = [
  { name: 'phone', width: WIDTH, height: 800, deviceScaleFactor: 1, mobile: true },
  { name: '400% zoom', width: WIDTH, height: 256, deviceScaleFactor: 4, mobile: false },
]

const CHROMES = [
  process.env['CHROME'],
  'C:/Program Files/Google/Chrome/Application/chrome.exe',
  'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/usr/bin/google-chrome',
  '/usr/bin/chromium',
]

type Theme = 'light' | 'dark' | 'dark, by the toggle'
const THEMES: Theme[] = ['light', 'dark', 'dark, by the toggle']

/** A page in Chrome, spoken to over the DevTools protocol. */
class Tab {
  #socket: WebSocket
  #next = 0
  #waiting = new Map<number, { resolve: (value: any) => void; reject: (error: Error) => void }>()
  #listeners = new Map<string, ((params: any) => void)[]>()

  private constructor(socket: WebSocket) {
    this.#socket = socket
    socket.addEventListener('message', (event) => {
      const message = JSON.parse(String(event.data))
      if (message.id !== undefined) {
        const waiting = this.#waiting.get(message.id)
        this.#waiting.delete(message.id)
        if (message.error) waiting?.reject(new Error(message.error.message))
        else waiting?.resolve(message.result)
      } else {
        for (const listener of this.#listeners.get(message.method) ?? []) listener(message.params)
      }
    })
  }

  static async open(url: string): Promise<Tab> {
    const socket = new WebSocket(url)
    await new Promise((resolve, reject) => {
      socket.addEventListener('open', resolve, { once: true })
      socket.addEventListener('error', reject, { once: true })
    })
    return new Tab(socket)
  }

  on(method: string, listener: (params: any) => void) {
    this.#listeners.set(method, [...(this.#listeners.get(method) ?? []), listener])
  }

  send(method: string, params: object = {}): Promise<any> {
    const id = this.#next++
    this.#socket.send(JSON.stringify({ id, method, params }))
    return new Promise((resolve, reject) => this.#waiting.set(id, { resolve, reject }))
  }

  async navigate(url: string): Promise<void> {
    let loaded = () => {}
    const done = new Promise<void>((resolve) => (loaded = resolve))
    this.#listeners.set('Page.loadEventFired', [() => loaded()])
    await this.send('Page.navigate', { url })
    await done
  }

  async evaluate<T>(expression: string): Promise<T> {
    const { result, exceptionDetails } = await this.send('Runtime.evaluate', {
      expression,
      awaitPromise: true,
      returnByValue: true,
    })
    if (exceptionDetails) throw new Error(exceptionDetails.exception?.description ?? exceptionDetails.text)
    return result.value as T
  }

  close() {
    this.#socket.close()
  }
}

// What sticks out past the right edge of the viewport, and is not inside something that
// scrolls on its own (a code block, a wide table). Named by tag and class, the first few.
const OVERFLOW = `(() => {
  const width = document.documentElement.clientWidth
  if (document.documentElement.scrollWidth <= width) return []
  const scrolls = (el) => { for (let p = el.parentElement; p; p = p.parentElement) {
    if (/(auto|scroll|hidden|clip)/.test(getComputedStyle(p).overflowX)) return true } return false }
  return [...document.body.querySelectorAll('*')]
    .filter((el) => el.getBoundingClientRect().right > width + 0.5 && !scrolls(el))
    .slice(0, 5)
    .map((el) => el.tagName.toLowerCase() + (el.className && typeof el.className === 'string' ? '.' + el.className.trim().split(/\\s+/).join('.') : ''))
})()`

const chrome = CHROMES.find((path) => path && existsSync(path))
if (!chrome) throw new Error('No Chrome found; set CHROME to its path.')
const profile = await mkdtemp(join(tmpdir(), 'qdecr-qa-'))
const browser = spawn(chrome, [
  '--headless=new',
  '--remote-debugging-port=0',
  `--user-data-dir=${profile}`,
  '--no-first-run',
  'about:blank',
])

let failed = false
try {
  // Chrome writes the port it picked into the profile once it listens.
  let port = ''
  for (let tries = 0; !port && tries < 100; tries++) {
    await new Promise((resolve) => setTimeout(resolve, 100))
    port = (await readFile(join(profile, 'DevToolsActivePort'), 'utf8').catch(() => '')).split('\n')[0] ?? ''
  }
  const target = await (await fetch(`http://127.0.0.1:${port}/json/new?about:blank`, { method: 'PUT' })).json()
  const tab = await Tab.open(target.webSocketDebuggerUrl)
  await tab.send('Page.enable')
  // What the page logs as an error: its own console.error, an uncaught exception, and the
  // browser's own reports (a CSP violation, a failed request).
  let errors: string[] = []
  await tab.send('Runtime.enable')
  await tab.send('Log.enable')
  tab.on('Runtime.exceptionThrown', ({ exceptionDetails }) =>
    errors.push(`threw ${exceptionDetails.exception?.description ?? exceptionDetails.text}`),
  )
  tab.on('Runtime.consoleAPICalled', ({ type, args }) => {
    if (type === 'error') errors.push(`logged an error: ${args.map((arg: any) => arg.value ?? arg.description).join(' ')}`)
  })
  tab.on('Log.entryAdded', ({ entry }) => {
    if (entry.level === 'error') errors.push(`logged an error: ${entry.text}`)
  })
  // Everything the site loads comes from itself; data: URLs are the
  // viewer's own font and lighting (the CSP's img-src).
  await tab.send('Network.enable')
  tab.on('Network.requestWillBeSent', ({ request }) => {
    const url: string = request.url
    if (!url.startsWith(ORIGIN) && !url.startsWith('data:')) errors.push(`requested ${url} from another origin`)
  })
  const axe = await readFile(createRequire(import.meta.url).resolve('axe-core/axe.min.js'), 'utf8')

  /** What is wrong with the page as it stands, one line each. */
  async function check(dark: boolean): Promise<string[]> {
    const problems: string[] = []
    // Judged by the paper the page is actually painted on, so a dark block of tokens
    // that failed to apply shows up here.
    const isDark = await tab.evaluate<boolean>(
      `(() => { const [r, g, b] = getComputedStyle(document.body).backgroundColor.match(/\\d+/g).map(Number); return r + g + b < 384 })()`,
    )
    if (isDark !== dark) problems.push(`the page is not ${dark ? 'dark' : 'light'}`)
    const overflow = await tab.evaluate<string[]>(OVERFLOW)
    if (overflow.length) problems.push(`scrolls sideways at ${WIDTH}px: ${overflow.join(', ')}`)
    const violations = await tab.evaluate<{ id: string; help: string; nodes: number; target: string }[]>(
      `(async () => { ${axe}; const r = await axe.run(document, { resultTypes: ['violations'] });
        return r.violations.map((v) => ({ id: v.id, help: v.help, nodes: v.nodes.length, target: String(v.nodes[0]?.target ?? '') })) })()`,
    )
    for (const v of violations) problems.push(`axe ${v.id} (${v.nodes}×, first ${v.target}): ${v.help}`)
    for (const error of errors) problems.push(error)
    errors = []
    return problems
  }

  /**
   * Presses "Explore in 3D" and waits for the viewer: true once it is up, false on a page
   * without one. The button only shows where the browser has WebGL 2, so a page with a
   * viewer and no button means the check could not run, and fails.
   */
  async function openViewer(): Promise<boolean> {
    const found = await tab.evaluate<'none' | 'no button' | 'pressed'>(
      `(() => { if (!document.querySelector('[data-viewer]')) return 'none'
        const b = document.querySelector('.viewer-open'); if (!b || b.hidden) return 'no button'
        b.click(); return 'pressed' })()`,
    )
    if (found === 'none') return false
    if (found === 'no button') throw new Error('The 3D viewer shows no button: this Chrome has no WebGL 2')
    for (let tries = 0; tries < 300; tries++) {
      const state = await tab.evaluate<string>(
        `document.querySelector('[data-viewer].is-live') ? 'live' : document.querySelector('.viewer-status')?.textContent ?? ''`,
      )
      if (state === 'live') return true
      if (state) throw new Error(`The viewer did not open: ${state}`)
      await new Promise((resolve) => setTimeout(resolve, 100))
    }
    throw new Error('The viewer did not open within 30 seconds')
  }

  const pages = await pagesToCheck()
  for (const { name: viewportName, ...viewport } of VIEWPORTS) {
    await tab.send('Emulation.setDeviceMetricsOverride', viewport)
    for (const theme of THEMES) {
      const dark = theme !== 'light'
      await tab.send('Emulation.setEmulatedMedia', {
        features: [{ name: 'prefers-color-scheme', value: theme === 'dark' ? 'dark' : 'light' }],
      })
      // The toggle's choice is kept in localStorage, which public/theme.js reads before
      // the first paint; setting it there is what pressing the toggle does.
      const choice = await tab.send('Page.addScriptToEvaluateOnNewDocument', {
        source: theme === 'dark, by the toggle' ? `localStorage.setItem('theme', 'dark')` : `localStorage.removeItem('theme')`,
      })
      const report = (page: string, problems: string[]) => {
        if (problems.length) failed = true
        const where = `${viewportName}, ${theme}:`.padEnd(31)
        console.log(`${where} ${page}${problems.length ? '\n  ' + problems.join('\n  ') : ' ok'}`)
      }
      for (const page of pages) {
        errors = []
        await tab.navigate(ORIGIN + page)
        report(page, await check(dark))
        if (await openViewer()) report(`${page}, 3D viewer open`, await check(dark))
      }
      await tab.send('Page.removeScriptToEvaluateOnNewDocument', { identifier: choice.identifier })
    }
  }
  tab.close()
} finally {
  browser.kill()
  await rm(profile, { recursive: true, force: true }).catch(() => {})
}
if (failed) process.exitCode = 1
