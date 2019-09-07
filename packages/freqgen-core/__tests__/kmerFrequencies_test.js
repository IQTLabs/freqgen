var kmerFrequencies = require('..').kmerFrequencies

describe('k-mer frequencies', function() {
  test('Base case with empty k-mer count object', function() {
    expect(kmerFrequencies(new Map())).toEqual(new Map())
  })

  test('A single k-mer should have 100% frequency', function() {
    expect(kmerFrequencies(new Map([['A', 1]]))).toEqual(new Map([['A', 1]]))
  })

  test('A two k-mers should have 50% frequency each', function() {
    expect(kmerFrequencies(new Map([['A', 1], ['T', 1]]))).toEqual(
      new Map([['A', 0.5], ['T', 0.5]])
    )
  })

  test('An object with different k values should throw an error', () => {
    expect(() => {
      kmerFrequencies(new Map([['A', 1], ['AA', 1]]))
    }).toThrow()
  })

  test('An Object (instead of a Map) will throw an error', () => {
    expect(() => kmerFrequencies({ A: 1, AA: 1 })).toThrow()
  })

  test('An object with different k values should not throw an error if validation is off', () => {
    expect(
      kmerFrequencies(new Map([['A', 1], ['AA', 1]]), { validation: false })
    ).toEqual(new Map([['A', 0.5], ['AA', 0.5]]))
  })
})
