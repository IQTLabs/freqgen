var utilities = module.exports

/* This is an internal method use to ensure that a k-mer Map only contains
   k-mers of the same size as keys. For example, {A: 1, T: 1} would be fine, but
   {A: 1, TT: 1} would not. */
utilities.validateKmerCountMap = function(kmerCountMap) {
  let kmers = kmerCountMap.keys()
  let k = kmers.next().value.length
  for (let kmer of kmers) {
    if (kmer.length !== k) {
      throw new Error(
        `Not all k-mers are of length ${k}. Got ${kmer}, which is of length ${kmer.length}.`
      )
    }
  }
}

utilities.validateKmerFrequencyMap = function(kmerFrequencyMap) {
  kmerFrequencyMap.forEach((v, k) => {
    utilities.validateKmerCountMap(v) // first check to make sure that all the k-mers are of the same length
    // for float imprecision, make it fixed and compare to a string
    if (utilities.sumMapValues(v).toFixed(10) !== '1.0000000000') {
      throw new Error(`Values for k=${k} do not sum to 1.`)
    }
  })
}

// just a little utility to get the total of a map's values
utilities.sumMapValues = function(map) {
  let sum = 0
  for (let value of map.values()) {
    sum += value
  }
  return sum
}
