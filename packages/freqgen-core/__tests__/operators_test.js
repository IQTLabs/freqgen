const operators = require('../operators')

describe('Seed testing', () => {
  test('Ensure that the seed of one-codon AAs matches', () => {
    expect(new operators('FK', new Map(), 1, 2).seed()).toEqual([
      'TTTAAA',
      'TTTAAA',
    ])
  })
})

describe('Crossover testing', () => {
  test('One codon seqs should just get swapped', () => {
    expect(
      new operators(null, new Map()).crossover('ATG', 'GTA')
    ).toIncludeSameMembers(['GTA', 'ATG'])
  })

  test('Two codon seqs should have second codons swapped', () => {
    expect(
      new operators(null, new Map()).crossover('AAAAAA', 'TTTTTT')
    ).toIncludeSameMembers(['AAATTT', 'TTTAAA'])
  })

  test('Three codon seqs should be in list of possibilities', () => {
    expect([
      'AAAAAATTT',
      'AAATTTTTT',
      'TTTTTTAAA',
      'TTTAAAAAA',
    ]).toIncludeAllMembers(
      new operators(null, new Map()).crossover('AAAAAAAAA', 'TTTTTTTTT')
    )
  })
})

describe('Mutation testing', () => {
  test("Fail gracefully when sequence can't be mutated since there are no synonyms", () => {
    expect(new operators(null, new Map(), 1).mutate('ATG')).toEqual('ATG')
  })

  test('If a seq has one codon with one synonym, mutate should return that synonym', () => {
    expect(new operators(null, new Map(), 1).mutate('TGT')).toEqual('TGC')
  })
})
