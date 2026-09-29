// Numbers from the example run as the site's pages print them: the home page's table,
// and the tutorials' prose, which a test holds to run.json (example-prose.test.ts).

/** A whole number with thousands separated, as English prints it: 113,272. */
export const count = (n: number) => Math.round(n).toLocaleString('en-GB')

/** A number to a fixed number of decimals, with a true minus sign: −0.029. */
export const signed = (n: number, digits: number) => n.toFixed(digits).replace('-', '−')
