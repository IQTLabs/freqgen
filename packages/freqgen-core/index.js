var freqgen = module.exports
freqgen.generate = require('./generate')

// configure warnings a la parcel/packages/core/test-utils/src/utils.js#L15 @ 30624f7
const chalk = require('chalk')
const warning = chalk.keyword('orange')
console.warn = (...args) => {
  console.error(warning(...args))
}

const validateKmerMap = require('./validateKmerMap')

freqgen.kmers = function (seq, k, { overlap = true } = {}) {
  if (k === undefined || k < 1) {
    throw new Error('k value >= 0 is required')
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
  let counts = new Map()

  let k = kmers[0] ? kmers[0].length : 0

  for (let i = 0; i < kmers.length; i++) {
    let kmer = kmers[i]
    if (kmer.length !== k) {
      throw new Error(`Not all k-mers are of length ${k}. At index ${i}, got ${kmer}, which is of length ${kmer.length}.`)
    }
    counts.set(kmer, (counts.get(kmer) || 0) + 1) // counts.get(kmer) || 0 is 0 or the k-mer count, whichever is greater
  }
  return counts
}

freqgen.kmerFrequencies = function (counts, { validation = true } = {}) {
  /* Optionally check that all of the k-mers in counts Object are of the same
     length. Because this adds overhead, we can skip it if we generate the
     counts Object from freqgen.kmerCounts since it already did the checking. */

  if (validation) {
    if (!(counts instanceof Map)) {
      counts = new Map(Object.entries(counts))
    }

    if (counts.size === 0) {
      return new Map()
    }

    validateKmerMap(counts)
  }

  // sum up all of the values in the Map
  let totalKmers = 0
  for (let kmerCount of counts.values()) {
    totalKmers += kmerCount
  }

  counts.forEach((v, k, m) => m.set(k, v / totalKmers))

  return counts
}

// a helper function that utilizes the performance boost of skipping validation in kmerFrequencies
freqgen.kmerFrequenciesFromSeq = (seq, k) => freqgen.kmerFrequencies(freqgen.kmerCounts(freqgen.kmers(seq, k)), { validation: false })
