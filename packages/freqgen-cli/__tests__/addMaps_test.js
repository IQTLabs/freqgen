const addMaps = require('../addMaps')

describe('addMaps testing', () => {
  test('Basic test', () => {
    expect(
      addMaps(new Map([['A', 2], ['T', 4]]), new Map([['A', 1], ['C', 1]]))
    ).toEqual(new Map([['A', 3], ['C', 1], ['T', 4]]))
  })
})
