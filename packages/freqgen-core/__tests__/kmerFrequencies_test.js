var kmerFrequencies = require('..').kmerFrequencies

describe('k-mer frequencies', function () {
  test('Base case with empty k-mer count object', function () {
    expect(kmerFrequencies({})).toEqual(new Map())
  })

  test('A single k-mer should have 100% frequency', function () {
    expect(kmerFrequencies({A: 1})).toEqual(new Map(Object.entries({A: 1})))
  })

  test('A two k-mers should have 50% frequency each', function () {
    expect(kmerFrequencies({A: 1, T: 1})).toEqual(new Map(Object.entries({A: 0.5, T: 0.5})))
  })

  test('An object with different k values should throw an error', () => {
    expect(() => { kmerFrequencies({A: 1, AA: 1}) }).toThrow()
  })

  test('Objects and Maps should be treated equally by default', () => {
    expect(kmerFrequencies({A: 1, T: 1}))
      .toEqual(kmerFrequencies(new Map([['A', 1], ['T', 1]])))
  })

  test('An Object (instead of a Map) will throw an error if validation is off', () => {
    expect(() => kmerFrequencies({A: 1, AA: 1}, {validation: false})).toThrowError(TypeError)
  })

  test('An object with different k values should not throw an error if validation is off', () => {
    expect(kmerFrequencies(new Map([['A', 1], ['AA', 1]]), {validation: false})).toEqual(new Map(Object.entries({A: 0.5, AA: 0.5})))
  })
})
