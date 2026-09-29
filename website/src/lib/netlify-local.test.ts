import { test } from 'node:test'
import assert from 'node:assert/strict'
import { candidateFiles, headersFor, redirectFor, type HeaderRule } from './netlify-local.ts'

const redirects = [
  { from: 'https://www.qdecr.com/*', to: 'https://qdecr.com/:splat', status: 301, force: true },
  { from: '/about.html', to: '/about', status: 301, force: true },
  { from: '/ohbm2020', to: '/about#archive', status: 301, force: true },
  { from: '/old', to: '/new' },
]

test('a redirect matches its path exactly, and keeps its status', () => {
  assert.deepEqual(redirectFor('/about.html', redirects), { to: '/about', status: 301, force: true })
  assert.deepEqual(redirectFor('/ohbm2020', redirects), { to: '/about#archive', status: 301, force: true })
  assert.equal(redirectFor('/about', redirects), undefined)
  assert.equal(redirectFor('/about.html/x', redirects), undefined)
})

test('a redirect without a status is a 301 and without force is not forced, as on Netlify', () => {
  assert.deepEqual(redirectFor('/old', redirects), { to: '/new', status: 301, force: false })
})

test('rules for another host never match a local path', () => {
  assert.equal(redirectFor('/*', redirects), undefined)
  assert.equal(redirectFor('/', redirects), undefined)
})

const headers: HeaderRule[] = [
  { for: '/*', values: { 'X-Frame-Options': 'DENY', 'Cache-Control': 'public, max-age=0' } },
  { for: '/_astro/*', values: { 'Cache-Control': 'public, max-age=31536000, immutable' } },
  { for: '/robots.txt', values: { 'X-Robots': 'yes' } },
]

test('every matching header rule applies, a later one overriding an earlier one', () => {
  assert.deepEqual(headersFor('/', headers), { 'X-Frame-Options': 'DENY', 'Cache-Control': 'public, max-age=0' })
  assert.deepEqual(headersFor('/_astro/a.css', headers), {
    'X-Frame-Options': 'DENY',
    'Cache-Control': 'public, max-age=31536000, immutable',
  })
  assert.deepEqual(headersFor('/robots.txt', headers), {
    'X-Frame-Options': 'DENY',
    'Cache-Control': 'public, max-age=0',
    'X-Robots': 'yes',
  })
})

test('a splat matches anything below its folder, not a name that merely starts the same', () => {
  assert.equal(headersFor('/_astro', headers)['Cache-Control'], 'public, max-age=0')
  assert.equal(headersFor('/_astronaut.css', headers)['Cache-Control'], 'public, max-age=0')
})

test('a path is tried as a file, then with .html, then as a folder, as Netlify does', () => {
  assert.deepEqual(candidateFiles('/'), ['index.html'])
  assert.deepEqual(candidateFiles('/about'), ['about', 'about.html', 'about/index.html'])
  assert.deepEqual(candidateFiles('/tutorials/plotting'), [
    'tutorials/plotting',
    'tutorials/plotting.html',
    'tutorials/plotting/index.html',
  ])
  assert.deepEqual(candidateFiles('/cite/'), ['cite/index.html', 'cite.html'])
  assert.deepEqual(candidateFiles('/_astro/a.css'), ['_astro/a.css'])
})

test('a path is decoded, and one that climbs out of the site finds nothing', () => {
  assert.deepEqual(candidateFiles('/archive/my%20poster.pdf'), ['archive/my poster.pdf'])
  assert.deepEqual(candidateFiles('/../secret'), [])
  assert.deepEqual(candidateFiles('/a/%2e%2e/%2e%2e/secret'), [])
  assert.deepEqual(candidateFiles('/%zz'), [])
})

test('a backslash, which Windows reads as a separator, finds nothing either', () => {
  assert.deepEqual(candidateFiles('/..%5c..%5csecret'), [])
  assert.deepEqual(candidateFiles('/a%5cb.html'), [])
})
