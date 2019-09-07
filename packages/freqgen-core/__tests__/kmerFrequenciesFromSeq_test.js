const kmerFrequenciesFromSeq = require('..').kmerFrequenciesFromSeq

describe('Test kmerFrequenciesFromSeq shortcut', () => {
  test('Basic test', () => {
    expect(kmerFrequenciesFromSeq('ATGC', [1])).toEqual(
      new Map([
        [1, new Map([['A', 0.25], ['T', 0.25], ['G', 0.25], ['C', 0.25]])],
      ])
    )
  })

  test('Multiple k-values', () => {
    expect(kmerFrequenciesFromSeq('ATGC', [1, 4])).toEqual(
      new Map([
        [1, new Map([['A', 0.25], ['T', 0.25], ['G', 0.25], ['C', 0.25]])],
        [4, new Map([['ATGC', 1]])],
      ])
    )
  })
})
