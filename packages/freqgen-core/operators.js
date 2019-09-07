const yaml = require('js-yaml')
const fs = require('fs')
const kmers = require('./kmers')
const distance = require('./distance')
const path = require('path')
const random = require('lodash.random')

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
  constructor(targetAminoAcidSeq, targetFreqs, geneticCode, populationSize) {
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
      let dnaSeq = ''

      for (let letter of this.targetAminoAcidSeq) {
        dnaSeq += this.codonsForAminoAcid[letter][0]
      }

      return Array.from({ length: this.populationSize }, () => dnaSeq)
    }

    this.fitness = seq => {
      // first, convert it to a flat map (e.g. {A => 0.5, T => 0.5, AT => 1.0})
      let freqs = kmers.kmerFrequenciesFromSeq(seq, this.k)
      freqs = new Map(
        [...freqs.values()].map(x => Array.from(x.entries())).flat()
      )

      // then compare with flat map we precalculated
      return distance.cosine(freqs, this.targetFreqsFlat)
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
