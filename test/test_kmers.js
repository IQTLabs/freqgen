var assert = require('assert')
var expect = require('chai').expect
var should = require('chai').should()

var kmers = require('../index').kmers

describe('k-mer calculation', function () {
  context('with overlap', function () {
    it('k=1', function () {
      let result = kmers('GATTACA', 1)
      result.should.eql(['G', 'A', 'T', 'T', 'A', 'C', 'A'])
    })

    it('k=2', function () {
      let result = kmers('GATTACA', 2)
      result.should.eql(['GA', 'AT', 'TT', 'TA', 'AC', 'CA'])
    })

    it('k=3', function () {
      let result = kmers('GATTACA', 3)
      result.should.eql(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
    })

    it('k=4', function () {
      let result = kmers('GATTACA', 4)
      result.should.eql(['GATT', 'ATTA', 'TTAC', 'TACA'])
    })
  })

  context('no overlap', function () {
    it('k=1', function () {
      let result = kmers('GATTACA', 1, overlap = false)
      result.should.eql(['G', 'A', 'T', 'T', 'A', 'C', 'A'])
    })

    it('k=2 with end cut off', function () {
      let result = kmers('GATTACA', 2, overlap = false)
      result.should.eql(['GA', 'TT', 'AC'])
    })

    it('k=2', function () {
      let result = kmers('GATTACAT', 2, overlap = false)
      result.should.eql(['GA', 'TT', 'AC', 'AT'])
    })

    it('k=3 with end cut off', function () {
      let result = kmers('GATTACA', 3, overlap = false)
      result.should.eql(['GAT', 'TAC'])
    })

    it('k=4', function () {
      let result = kmers('GATTACAT', 4, overlap = false)
      result.should.eql(['GATT', 'ACAT'])
    })
  })
})
