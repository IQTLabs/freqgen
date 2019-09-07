const utilities = require('./utilities')

module.exports.kmers = function(seq, k, { overlap = true } = {}) {
  if (k === undefined || k < 1) {
    throw new Error('k value >= 0 is required')
  }

  if (seq.length % k != 0 && !overlap) {
    throw new Error(
      `Sequence is of length ${seq.length}, which is not unabiguously divisible into overlapping ${k}-mers.`
    )
  }
  let numKmers = overlap ? seq.length - k + 1 : Math.floor(seq.length / k)
  let result = new Array(numKmers)
  let spacing = overlap ? 1 : k
  for (let i = 0; i < seq.length - k + 1; i += spacing) {
    result[Math.floor(i / spacing)] = seq.substring(i, i + k)
  }
  return result
}

module.exports.kmerCounts = function(kmers) {
  let counts = new Map()

  let k = kmers[0] ? kmers[0].length : 0

  for (let i = 0; i < kmers.length; i++) {
    let kmer = kmers[i]
    if (kmer.length !== k) {
      throw new Error(
        `Not all k-mers are of length ${k}. At index ${i}, got ${kmer}, which is of length ${kmer.length}.`
      )
    }
    counts.set(kmer, (counts.get(kmer) || 0) + 1) // counts.get(kmer) || 0 is 0 or the k-mer count, whichever is greater
  }
  return counts
}

module.exports.kmerFrequencies = function(counts, { validation = true } = {}) {
  /* Optionally check that all of the k-mers in counts Object are of the same
     length. Because this adds overhead, we can skip it if we generate the
     counts Object from module.exports.kmerCounts since it already did the checking. */

  if (validation) {
    if (counts.size === 0) {
      return new Map()
    }

    utilities.validateKmerCountMap(counts)
  }

  let totalKmers = utilities.sumMapValues(counts)

  counts.forEach((v, k, m) => m.set(k, v / totalKmers))

  return counts
}

// a helper function that utilizes the performance boost of skipping validation in kmerFrequencies
module.exports.kmerFrequenciesFromSeq = function(seq, k) {
  let result = new Map()
  for (let _k of k) {
    if (_k == 'codons') {
      _k = 3
      overlap = { overlap: false }
    } else {
      overlap = { overlap: true }
    }
    result.set(
      _k,
      module.exports.kmerFrequencies(
        module.exports.kmerCounts(module.exports.kmers(seq, _k, overlap)),
        {
          validation: false,
        }
      )
    )
  }
  return result
}
