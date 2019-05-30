/* This is an internal method use to ensure that a k-mer Map only contains
/* k-mers of the same size as keys. For example, {A: 1, T: 1} would be fine, but {A: 1, TT: 1} would not. */
module.exports = function(kmerMap) {
  let kmers = kmerMap.keys()
  let k = kmers.next().value.length
  for (let kmer of kmers) {
    if (kmer.length !== k) {
      throw new Error(
        `Not all k-mers are of length ${k}. Got ${kmer}, which is of length ${
          kmer.length
        }.`
      )
    }
  }
}
