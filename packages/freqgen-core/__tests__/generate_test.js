const generate = require('..').generate

describe('Sequence generation', () => {
  describe('Throw if all target frequencies do not sum to 1', () => {
    test('One bad k-value', () => {
      expect(() => generate('M', { 1: { A: 1, T: 1 } })).toThrow()
    })

    test('Two k values, one of which is bad', () => {
      let targetFreqs = new Map([
        [1, new Map([['A', 1], ['T', 1]])],
        [2, new Map([['AA', 0.5], ['TT', 0.5]])],
      ])
      expect(() => generate('M', targetFreqs)).toThrow()
    })
  })

  // test('A simple 1 mer test', () => {
  //   let targetFreqs = new Map([
  //     [1, new Map([['A', 0.5], ['T', 0.5], ['G', 0], ['C', 0]])],
  //   ])
  //   expect(generate('FK', targetFreqs)).toEqual('TTTAAA')
  // })
})
