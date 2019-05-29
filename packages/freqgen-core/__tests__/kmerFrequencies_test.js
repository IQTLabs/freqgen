var kmerFrequencies = require('..').kmerFrequencies

// needed for testing warning output
global.console.error = jest.fn()
const chalk = require('chalk')
const warning = chalk.keyword('orange')

describe('k-mer frequencies', function () {
  test('Base case with empty k-mer count object', function () {
    expect(kmerFrequencies({})).toEqual({})
  })

  test('A single k-mer should have 100% frequency', function () {
    expect(kmerFrequencies({A: 1})).toEqual({A: 1})
  })

  test('A two k-mers should have 50% frequency each', function () {
    expect(kmerFrequencies({A: 1, T: 1})).toEqual({A: 0.5, T: 0.5})
  })

  test('An object with different k values should throw an error', () => {
    expect(() => { kmerFrequencies({A: 1, AA: 1}) }).toThrow()
  })

  test('An object with different k values should not throw an error if validation is off', () => {
    expect(kmerFrequencies({A: 1, AA: 1}, {validation: false})).toEqual({A: 0.5, AA: 0.5})
  })

  test('An object with different k values should not throw an error if validation is off but will warn the user if verbose', () => {
    kmerFrequencies({A: 1, AA: 1}, {validation: false, verbose: true})
    expect(global.console.error).toHaveBeenCalledWith(warning("Skipping validation... I hope you know what you're doing!"))
  })
})
