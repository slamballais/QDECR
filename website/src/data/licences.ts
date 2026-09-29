// The licences of the code and fonts this site sends to the browser, printed on /licences.
// The build fails when its client bundle carries a package that is not listed here, or
// when a package listed as bundled has left it (src/lib/licences.ts), so adding a
// dependency to a page's scripts means adding its notice. The text of each licence is read
// from the package's own file at build time, so it is always the text of the version
// shipped; the two packages that ship no such file have theirs here.

import type { NoticeScope } from '../lib/licences.ts'

/** Where a notice's licence text comes from: the package's own file, or this file. */
export type LicenceSource =
  /** The licence file, relative to node_modules/<first package>/. */
  | { file: string }
  /** The licence itself, for a package that ships no file. */
  | { text: string }

export type Notice = NoticeScope &
  LicenceSource & {
    /** The name the page gives it. */
    name: string
    url: string
    /** The SPDX identifier. */
    licence: string
  }

export interface NoticeGroup {
  /** The heading on /licences: what on the site this code does. */
  title: string
  notices: Notice[]
}

const MIT = (copyright: string) => `MIT License

${copyright}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.`

// NiiVue's npm package has no licence file; this is the LICENSE of its repository.
const NIIVUE_LICENCE = `BSD 2-Clause License

Copyright (c) 2021, Niivue

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.`

export const noticeGroups: NoticeGroup[] = [
  {
    title: 'The 3D viewer',
    notices: [
      { name: 'NiiVue', packages: ['@niivue/niivue'], bundled: true, url: 'https://github.com/niivue/niivue', licence: 'BSD-2-Clause', text: NIIVUE_LICENCE },
      { name: 'gl-matrix', packages: ['gl-matrix'], bundled: true, url: 'https://github.com/toji/gl-matrix', licence: 'MIT', file: 'LICENSE.md' },
      { name: 'fflate', packages: ['fflate'], bundled: true, url: 'https://github.com/101arrowz/fflate', licence: 'MIT', file: 'LICENSE' },
      { name: 'NIFTI-Reader-JS', packages: ['nifti-reader-js'], bundled: true, url: 'https://github.com/rii-mango/NIFTI-Reader-JS', licence: 'MIT', file: 'LICENSE' },
      { name: 'zarrita', packages: ['zarrita'], bundled: true, url: 'https://github.com/manzt/zarrita.js', licence: 'MIT', file: 'LICENSE' },
      { name: '@zarrita/storage', packages: ['@zarrita/storage'], bundled: true, url: 'https://github.com/manzt/zarrita.js', licence: 'MIT', file: 'LICENSE' },
      { name: 'structured-clone', packages: ['@ungap/structured-clone'], bundled: true, url: 'https://github.com/ungap/structured-clone', licence: 'ISC', file: 'LICENSE' },
      { name: 'uuid', packages: ['@lukeed/uuid'], bundled: true, url: 'https://github.com/lukeed/uuid', licence: 'MIT', file: 'license' },
      { name: 'array-equal', packages: ['array-equal'], bundled: true, url: 'https://github.com/sindresorhus/array-equal', licence: 'MIT', file: 'LICENSE' },
    ],
  },
  {
    title: 'Search',
    notices: [
      // Copied into dist/pagefind/ by the Pagefind CLI after Astro has built, so never in
      // the bundle. The file also names vscode-ripgrep, which only the CLI uses.
      { name: 'Pagefind', packages: ['pagefind'], bundled: false, url: 'https://pagefind.app/', licence: 'MIT', file: 'LICENSE/LICENSE' },
    ],
  },
  {
    title: 'The bundler’s helpers',
    notices: [
      // Vite's LICENSE.md goes on to the licences of Vite's own dependencies, none of which
      // reach the browser, so only its first part is here.
      { name: 'Vite', packages: ['vite'], bundled: true, url: 'https://vite.dev/', licence: 'MIT', text: MIT('Copyright (c) 2019-present, VoidZero Inc. and Vite contributors') },
      { name: 'Rolldown', packages: ['rolldown'], bundled: true, url: 'https://rolldown.rs/', licence: 'MIT', file: 'LICENSE' },
    ],
  },
  {
    title: 'Fonts',
    // Newsreader and Plex Sans are served cut down to the weights the site uses
    // (tools/subset-fonts.mjs), which makes them Modified Versions under the OFL. That is
    // allowed with the licence attached, and keeping their names is too: these releases
    // (Google Fonts', by way of Fontsource) reserve no font name.
    notices: [
      { name: 'Newsreader', packages: ['@fontsource-variable/newsreader'], bundled: false, url: 'https://github.com/productiontype/Newsreader', licence: 'OFL-1.1', file: 'LICENSE' },
      { name: 'IBM Plex Sans', packages: ['@fontsource-variable/ibm-plex-sans'], bundled: false, url: 'https://github.com/IBM/plex', licence: 'OFL-1.1', file: 'LICENSE' },
      { name: 'IBM Plex Mono', packages: ['@fontsource/ibm-plex-mono'], bundled: false, url: 'https://github.com/IBM/plex', licence: 'OFL-1.1', file: 'LICENSE' },
    ],
  },
]

export const notices: Notice[] = noticeGroups.flatMap((group) => group.notices)
