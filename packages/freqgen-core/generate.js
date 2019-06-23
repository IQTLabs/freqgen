const yaml = require('js-yaml')
const fs = require('fs')
const path = require('path')
const GenAlgo = require('GenAlgo')

// local packages
const kmers = require('./kmers')
const Operators = require('./operators')
const distance = require('./distance')
const utilities = require('./utilities')

// configure warnings a la parcel/packages/core/test-utils/src/utils.js#L15 @ 30624f7
const chalk = require('chalk')
const warning = chalk.keyword('orange')
console.warn = (...args) => {
  console.error(warning(...args))
}

const DNA_BASES = new Set(['A', 'T', 'G', 'C'])

module.exports = function generate(
  targetAminoAcidSeq,
  targetFreqs,
  {
    populationSize = 100,
    mutationProb = 0.3,
    crossoverProb = 0.8,
    maxGensSinceImprovement = 50,
    improvementRelThreshold = 0.0,
    geneticCode = 11,
    verbose = false,
    fitnessFunction = null,
  } = {}
) {
  // Do some quick sequence validation
  if (
    targetAminoAcidSeq.split('').every(aminoAcid => DNA_BASES.has(aminoAcid))
  ) {
    console.warn(
      'This appears to be a DNA sequence, not an amino acid sequence. Ensure that you are passing an amino acid sequence.'
    )
  }

  for (let k of targetFreqs.keys()) {
    if (!(typeof k == 'number' || k == 'codons')) {
      throw new Error(
        'targetFreqs must be a map from k values to frequency maps.'
      )
    }
  }

  utilities.validateKmerFrequencyMap(targetFreqs)

  let target = new Operators(
    targetAminoAcidSeq,
    targetFreqs,
    geneticCode,
    populationSize
  )

  const algo = new GenAlgo.GenAlgo({ iterationNumber: 1000000 })

  algo.setCrossoverFunction(target.crossover)
  algo.setFitnessEvaluator(target.fitness)
  algo.setMutationFunction(target.mutate)
  algo.setSeed(target.seed())

  // Will be called at each iteration
  const iterationCallback = ({
    bestIndividual,
    elapsedTime,
    iterationNumber,
  }) => {
    // if (iterationNumber % 10 == 0) {
    //   console.log('Iteration ' + iterationNumber)
    //   console.log('Best fitness : ' + bestIndividual.fitness)
    //   console.log('Elapsed time : ' + elapsedTime)
    //   console.log('Gens since improvement : ' + algo.gensSinceImprovement)
    // }

    if (bestIndividual.fitness > algo.bestFitness) {
      algo.bestFitness = bestIndividual.fitness
      algo.gensSinceImprovement = 0
    } else {
      algo.gensSinceImprovement += 1
    }
    if (algo.gensSinceImprovement > maxGensSinceImprovement) {
      return false
    }
    return true
  }
  algo.setIterationCallback(iterationCallback)

  algo.gensSinceImprovement = 0
  algo.bestFitness = 0

  algo
    .start()
    .then(res =>
      console.log(res[0], kmers.kmerFrequenciesFromSeq(res[0].entity, [1, 2]))
    )

  return
}
