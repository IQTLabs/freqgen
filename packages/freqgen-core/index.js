var freqgen = module.exports

const _ = require('lodash')
// configure warnings a la parcel/packages/core/test-utils/src/utils.js#L15 @ 30624f7
const chalk = require('chalk')
const warning = chalk.keyword('orange')
console.warn = (...args) => {
  console.error(warning(...args))
}

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
  let counts = {}

  let k = kmers[0] ? kmers[0].length : 0

  for (let i = 0; i < kmers.length; i++) {
    let kmer = kmers[i]
    if (kmer.length != k) {
      throw new Error(`Not all k-mers are of length ${k}. At index ${i}, got ${kmer}, which is of length ${kmer.length}.`)
    }
    counts[kmer] = (counts[kmer] || 0) + 1 // counts[kmer] || 0 is 0 or the k-mer count, whichever is greater
  }
  return counts
}

freqgen.kmerFrequencies = function (counts, {validation = true, verbose = false} = {}) {
  /* Optionally check that all of the k-mers in counts Object are of the same
     length. Because this adds overhead, we can skip it if we generate the
     counts Object from freqgen.kmerCounts since it already did the checking. */
  if (validation) {
    let kmers = Object.keys(counts)
    let kValues = kmers.map(key => key.length)
    k = kValues[0]
    for (let i = 0; i < kValues.length; i++) {
      let kmer = kmers[i]
      if (kmer.length != k) {
        throw new Error(`Not all k-mers are of length ${k}. Got ${kmer}, which is of length ${kmer.length}.`)
      }
    }
  } else if (verbose) {
    console.warn("Skipping validation... I hope you know what you're doing!")
  }

  let totalKmers = _.sum(Object.values(counts))
  return _.mapValues(counts, x => x / totalKmers)
}
