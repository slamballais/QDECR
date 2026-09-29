// The copy button on code blocks. The buttons are in the markup but hidden, and only shown
// here, so a reader without JavaScript never meets a control that does nothing.
//
// It copies the code and leaves R's printed output (the #> lines) behind: what a reader
// wants to paste into a console is what they would type.

const announcer = document.createElement('span')
announcer.className = 'visually-hidden'
announcer.setAttribute('role', 'status')
document.body.append(announcer)

function codeOf(pre: HTMLPreElement): string {
  const lines = pre.querySelectorAll<HTMLElement>('.line')
  if (!lines.length) return pre.textContent ?? ''
  return [...lines]
    .filter((line) => !line.classList.contains('output'))
    .map((line) => line.textContent ?? '')
    .join('\n')
}

function say(button: HTMLButtonElement, message: string) {
  button.textContent = message
  announcer.textContent = message
  window.setTimeout(() => {
    button.textContent = 'Copy'
    announcer.textContent = ''
  }, 2000)
}

for (const button of document.querySelectorAll<HTMLButtonElement>('.code-copy')) {
  const pre = button.parentElement?.querySelector('pre')
  if (!pre) continue
  button.hidden = false
  button.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(codeOf(pre))
      say(button, 'Copied')
    } catch {
      // No clipboard access (an insecure origin, or a browser that asks and was refused):
      // select the code so a keyboard shortcut finishes the job.
      const range = document.createRange()
      range.selectNodeContents(pre)
      const selection = window.getSelection()
      selection?.removeAllRanges()
      selection?.addRange(range)
      say(button, 'Selected')
    }
  })
}
