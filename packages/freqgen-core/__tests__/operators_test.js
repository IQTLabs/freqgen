const operators = require('../operators')

describe('Seed testing', () => {
  test('Ensure that the seed of one-codon AAs matches', () => {
    expect(new operators('FK', null, 1, { populationSize: 2 }).seed()).toEqual([
      'TTTAAA',
      'TTTAAA',
    ])
  })
})

describe('Crossover testing', () => {
  test('One codon seqs should just get swapped', () => {
    expect(new operators().crossover('ATG', 'GTA')).toIncludeSameMembers([
      'GTA',
      'ATG',
    ])
  })

  test('Two codon seqs should have second codons swapped', () => {
    expect(new operators().crossover('AAAAAA', 'TTTTTT')).toIncludeSameMembers([
      'AAATTT',
      'TTTAAA',
    ])
  })

  test('Three codon seqs should be in list of possibilities', () => {
    expect([
      'AAAAAATTT',
      'AAATTTTTT',
      'TTTTTTAAA',
      'TTTAAAAAA',
    ]).toIncludeAllMembers(new operators().crossover('AAAAAAAAA', 'TTTTTTTTT'))
  })
})

describe('Mutation testing', () => {
  test("Fail gracefully when sequence can't be mutated since there are no synonyms", () => {
    expect(new operators(null, null, 1).mutate('ATG')).toEqual('ATG')
  })

  test('If a seq has one codon with one synonym, mutate should return that synonym', () => {
    expect(new operators(null, null, 1).mutate('TGT')).toEqual('TGC')
  })
})
