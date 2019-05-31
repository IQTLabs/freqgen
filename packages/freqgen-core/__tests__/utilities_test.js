const utilities = require('../utilities')

describe('validateKmerFrequencyMap testing', () => {
  test('A positive test case', () => {
    let freqMap = new Map()
    freqMap.set(1, new Map([['A', 0.5], ['T', 0.5]]))
    freqMap.set(2, new Map([['AT', 0.5], ['TA', 0.5]]))
    expect(() => utilities.validateKmerFrequencyMap(freqMap)).not.toThrow()
  })

  test('A positive test case but with codons as the key', () => {
    let freqMap = new Map()
    freqMap.set('codons', new Map([['AAA', 0.5], ['TTT', 0.5]]))
    expect(() => utilities.validateKmerFrequencyMap(freqMap)).not.toThrow()
  })

  test("A negative test case where values don't add up to one", () => {
    let freqMap = new Map()
    freqMap.set(1, new Map([['A', 0.5], ['T', 0.6]])) // sums to 1.1!!
    freqMap.set(2, new Map([['AT', 0.5], ['TA', 0.5]]))
    expect(() => utilities.validateKmerFrequencyMap(freqMap)).toThrow()
  })
})

test("A negative test case where k-mers don't have matching lengths", () => {
  let freqMap = new Map()
  freqMap.set(1, new Map([['A', 0.5], ['TT', 0.5]])) // TT is a 2-mer!
  freqMap.set(2, new Map([['AT', 0.5], ['TA', 0.5]]))
  expect(() => utilities.validateKmerFrequencyMap(freqMap)).toThrow()
})

describe('Dot product for map testing', () => {
  test('Both maps have all keys in common', () => {
    let map1 = new Map([[0, 1], [1, 2], [2, 0], [3, 5]])
    let map2 = new Map([[0, 3], [1, 0], [2, 4], [3, 6]])

    expect(utilities.mapDotProduct(map1, map2)).toEqual(33)
  })

  test("The maps don't all keys in common", () => {
    let map1 = new Map([[0, 1], [1, 2], [3, 5]])
    let map2 = new Map([[0, 3], [2, 4], [3, 6]])

    expect(utilities.mapDotProduct(map1, map2)).toEqual(33)
  })
})
