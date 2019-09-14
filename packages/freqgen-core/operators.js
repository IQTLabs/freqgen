const yaml = require('js-yaml')
const fs = require('fs')
const kmers = require('./kmers')
const distance = require('./distance')
const path = require('path')
const random = require('lodash.random')
const memoize = require('fast-memoize')

const codonsForAminoAcid = yaml.load(
  fs.readFileSync(path.resolve(__dirname, './data/codons_for_aa.yaml'), 'utf8')
)

const synonmyousCodons = yaml.load(
  fs.readFileSync(
    path.resolve(__dirname, './data/synonymous_codons.yaml'),
    'utf8'
  )
)

const codonsWithoutSynonyms = yaml.load(
  fs.readFileSync(
    path.resolve(__dirname, './data/codons_without_synonyms.yaml'),
    'utf8'
  )
)

class Operators {
  constructor(
    targetAminoAcidSeq,
    targetFreqs,
    geneticCode,
    cache,
    { populationSize = 100 }
  ) {
    this.targetAminoAcidSeq = targetAminoAcidSeq
    this.targetFreqs = targetFreqs
    this.populationSize = populationSize

    this.targetFreqsFlat = new Map(
      [...targetFreqs.values()].map(x => Array.from(x.entries())).flat()
    )

    // store some useful mappings
    this.codonsForAminoAcid = codonsForAminoAcid[geneticCode]
    this.codonsWithoutSynonyms = codonsWithoutSynonyms[geneticCode]
    this.synonmyousCodons = synonmyousCodons[geneticCode]
    this.k = Array.from(targetFreqs.keys())

    this.seed = () => {
      let population = []
      for (let index = 0; index < populationSize; index++) {
        let dnaSeq = ''
        for (let letter of this.targetAminoAcidSeq) {
          dnaSeq += this.codonsForAminoAcid[letter][
            Math.floor(Math.random() * this.codonsForAminoAcid[letter].length)
          ]
        }
        population.push(dnaSeq)
      }
      return population
    }

    let fitness = seq => {
      // first, convert it to a flat map (e.g. {A => 0.5, T => 0.5, AT => 1.0})
      let freqs = kmers.kmerFrequenciesFromSeq(seq, this.k)
      freqs = new Map(
        [...freqs.values()].map(x => Array.from(x.entries())).flat()
      )

      // then compare with flat map we precalculated
      return distance.cosine(freqs, this.targetFreqsFlat)
    }

    // in the real world, use the fast-memoize package
    if (cache) {
      this.fitness = memoize(fitness)
    }
    // otherwise, use our own (probably slow) implementation that tracks stats
    // else if (cache && verbose) {
    //   let cacheMap = new Map()
    //   this.stats = { hits: 0, misses: 0 }
    //   this.fitness = seq => {
    //     if (cacheMap.has(seq)) {
    //       this.stats.hits++
    //     } else {
    //       this.stats.misses++
    //       cacheMap.set(seq, fitness(seq))
    //     }
    //     return cacheMap.get(seq)
    //   }
    // }
    // or, turn off caching altogether
    else {
      this.fitness = fitness
    }

    this.crossover = (parent_1, parent_2) => {
      let idx = random(1, parent_1.length / 3 - 1) * 3

      let child_1 = parent_1.slice(0, idx) + parent_2.slice(idx)
      let child_2 = parent_2.slice(0, idx) + parent_1.slice(idx)

      return [child_1, child_2]
    }

    this.mutate = sequence => {
      let codons = kmers.kmers(sequence, 3, { overlap: false })

      // if none of the codons have synonyms, fail fast
      if (
        codons.every(codon => this.codonsWithoutSynonyms.indexOf(codon) > -1)
      ) {
        return sequence
      }

      let codon = ''

      let idx = -1 // a dummy value
      do {
        idx = Math.floor(Math.random() * codons.length)
        codon = codons[idx]
      } while (this.codonsWithoutSynonyms.indexOf(codon) > -1)

      let choices = this.synonmyousCodons[codon].filter(key => key != codon)
      codons[idx] = choices[Math.floor(Math.random() * choices.length)]
      return codons.join('')
    }
  }
}

module.exports = Operators
