const distance = require('../distance')

describe('Euclidean distance', () => {
  test('Base case with empty maps', () => {
    expect(distance.euclidean(new Map(), new Map())).toEqual(0)
  })

  test('Fully defined maps', () => {
    map1 = new Map()
      .set(0, 1)
      .set(1, 2)
      .set(2, 0)
      .set(3, 5)
    map2 = new Map()
      .set(0, 3)
      .set(1, 0)
      .set(2, 4)
      .set(3, 6)
    expect(distance.euclidean(map1, map2)).toBeCloseTo(5)
  })

  test('Sparsely defined maps', () => {
    map1 = new Map()
      .set(0, 1)
      .set(1, 2)
      .set(3, 5)
    map2 = new Map()
      .set(0, 3)
      .set(2, 4)
      .set(3, 6)
    expect(distance.euclidean(map1, map2)).toBeCloseTo(5)
  })
})

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
