
var kmers = require('..').kmers

describe('k-mer calculation', function () {
  context('with overlap', function () {
    test('k=1', function () {
      let result = kmers('GATTACA', 1)
      expect(result).toEqual(['G', 'A', 'T', 'T', 'A', 'C', 'A'])
    })

    test('k=2', function () {
      let result = kmers('GATTACA', 2)
      expect(result).toEqual(['GA', 'AT', 'TT', 'TA', 'AC', 'CA'])
    })

    test('k=3', function () {
      let result = kmers('GATTACA', 3)
      expect(result).toEqual(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
    })

    test('k=4', function () {
      let result = kmers('GATTACA', 4)
      expect(result).toEqual(['GATT', 'ATTA', 'TTAC', 'TACA'])
    })
  })

  context('no overlap', function () {
    test('k=1', function () {
      let result = kmers('GATTACA', 1, overlap = false)
      result.should.eql(['G', 'A', 'T', 'T', 'A', 'C', 'A'])
    })

    test('k=2 with end cut off', function () {
      let result = kmers('GATTACA', 2, overlap = false)
      result.should.eql(['GA', 'TT', 'AC'])
    })

    test('k=2', function () {
      let result = kmers('GATTACAT', 2, overlap = false)
      result.should.eql(['GA', 'TT', 'AC', 'AT'])
    })

    test('k=3 with end cut off', function () {
      let result = kmers('GATTACA', 3, overlap = false)
      result.should.eql(['GAT', 'TAC'])
    })

    test('k=4', function () {
      let result = kmers('GATTACAT', 4, overlap = false)
      result.should.eql(['GATT', 'ACAT'])
    })
  })
})
