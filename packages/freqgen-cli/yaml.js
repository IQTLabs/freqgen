// So, why not use one of the many existing YAML writers? None of them support
// Maps, which Freqgen uses extensively. Sinice Maps can have arbitrary objects
// as keys, it's not possible to serialize them into YAML easily. However, for
// our purposes, our keys will only be numbers (k values, to be precise) and
// strings (only the string "codons"). Hence, we'll iterate over the keys in
// sorted order and output it to a depth of 1, since that's all we'll require.

module.exports = function(mapToSerialize) {
  let result = ''
  for (let entry of mapToSerialize.entries()) {
    result += `${entry[0]}:\n` // add the key with no indentation
    for (let kmer of [...entry[1].entries()].sort()) {
      result += `  ${kmer[0]}: ${kmer[1]}\n`
    }
  }
  return result.trimEnd()
}
