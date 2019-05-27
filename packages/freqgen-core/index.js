var freqgen = module.exports

var _ = require('lodash')

freqgen.kmers = function (seq, k, overlap = true) {
  if (k == null) {
    throw new Error('k value is required')
  }

  let numKmers = overlap ? seq.length - k + 1 : Math.floor(seq.length / k)
  let result = new Array(numKmers)
  let spacing = overlap ? 1 : k
  for (let i = 0; i < seq.length - k + 1; i += spacing) {
    result[Math.floor(i / spacing)] = seq.substring(i, i + k)
  }
  return result
}

freqgen.kmerCounts = function (kmers) {
  let counts = {}

  for (let i = 0; i < kmers.length; i++) {
    let kmer = kmers[i]
    counts[kmer] = (counts[kmer] || 0) + 1 // counts[kmer] || 0 is 0 or the k-mer count, whichever is greater
  }
  return counts
}

freqgen.kmerFrequencies = function (counts) {
  let totalKmers = _.sum(_.values(counts))
  return _.mapValues(counts, x => x / totalKmers)
}
