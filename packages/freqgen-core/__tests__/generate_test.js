const generate = require('..').generate

describe('Sequence generation', () => {
  describe('Throw if all target frequencies do not sum to 1', () => {
    test('One bad k-value', () => {
      expect(() => generate({1: {A: 1, T: 1}}, 'M')).toThrow()
    })

    test('Two k values, one of which is bad', () => {
      expect(() => generate({1: {A: 1, T: 1}, 2: {AA: 0.5, TT: 0.5}}, 'M')).toThrow()
    })
  })
})

