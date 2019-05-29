var kmerCounts = require('..').kmerCounts

describe('k-mer counting', function () {
  test('Base case with empty k-mer array', function () {
    expect(kmerCounts([])).toEqual({})
  })

  test('Count a single k-mer', function () {
    expect(kmerCounts(['A'])).toEqual({A: 1})
  })

  test('Count two different k-mers', function () {
    expect(kmerCounts(['A', 'T'])).toEqual({A: 1, T: 1})
  })

  test('Count a repeated k-mer', function () {
    expect(kmerCounts(['A', 'A'])).toEqual({A: 2})
  })

  test('Throw an error when k-mers aren\'t of equal length', function () {
    expect(() => kmerCounts(['A', 'AA'])).toThrow()
  })
})
