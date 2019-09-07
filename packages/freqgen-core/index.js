var freqgen = module.exports

// add all of the kmer methods to the top level object
Object.assign(freqgen, require('./kmers'))

freqgen.generate = require('./generate')
