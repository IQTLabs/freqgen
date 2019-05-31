const distance = require('../distance')

// describe('Euclidean distance', () => {
//
// })

describe('Cosine similarity', () => {
  test('Base case with empty maps', () => {
    expect(distance.cosine(new Map(), new Map())).toEqual(NaN)
  })

  test('Base case with identical maps', () => {
    expect(distance.cosine(new Map([[1, 1]]), new Map([[1, 1]]))).toEqual(1)
  })

  test('Fully defined maps', () => {
    expect(
      distance.cosine(
        new Map([[0, 2], [1, 0], [2, 1]]),
        new Map([[0, 1], [1, 0], [2, 1]])
      )
    ).toBeCloseTo((3 * 10 ** 0.5) / 10)
  })

  test('Sparsely defined maps (in which there are no keys with value = 0)', () => {
    expect(
      distance.cosine(new Map([[0, 2], [2, 1]]), new Map([[0, 1], [2, 1]]))
    ).toBeCloseTo((3 * 10 ** 0.5) / 10)
  })

  test('One more sparely defined map', () => {
    expect(
      distance.cosine(
        new Map([[0, 2], [1, 4], [2, 1]]),
        new Map([[0, 1], [2, 1]])
      )
    ).toBeCloseTo(42 ** 0.5 / 14)
  })
})
