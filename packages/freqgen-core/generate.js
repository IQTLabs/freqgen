const GenAlgo = require('GenAlgo')

// local packages
const Operators = require('./operators')
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
    mutationProbability = 0.3,
    crossoverProbability = 0.8,
    maxGensSinceImprovement = 50,
    maxGensTotal = 10000,
    geneticCode = 11,
    emitter = null,
    operators = null,
    cache = true,
    populationSize = 100,
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

  operators = new Operators(
    targetAminoAcidSeq,
    targetFreqs,
    geneticCode,
    cache,
    { populationSize }
  )
  const algo = new GenAlgo.GenAlgo({
    iterationNumber: maxGensTotal,
    crossoverProbability,
    mutationProbability,
  })

  algo.setCrossoverFunction(operators.crossover)
  algo.setSelectPairFunction(GenAlgo.tournament3Pair)
  algo.setFitnessEvaluator(operators.fitness)
  algo.setMutationFunction(operators.mutate)
  algo.setSeed(operators.seed())

  // Will be called at each iteration
  const iterationCallback = ({
    bestIndividual,
    elapsedTime,
    iterationNumber,
  }) => {
    // send an optional emitter information during optimization
    if (emitter !== null) {
      emitter.emit('generation', {
        iterationNumber,
        bestIndividualFitness: bestIndividual.fitness,
        elapsedTime,
        gensSinceImprovement: algo.gensSinceImprovement,
      })
    }
    // stop early if the fitness has stopped improving
    if (bestIndividual.fitness > algo.bestFitness) {
      algo.bestFitness = bestIndividual.fitness
      algo.gensSinceImprovement = 0
    } else {
      algo.gensSinceImprovement += 1
    }
    if (algo.gensSinceImprovement > maxGensSinceImprovement) {
      if (emitter != null) {
        emitter.emit('complete')
        // emitter.emit('stats', operators.stats)
      }
      return false
    }
    return true
  }
  algo.setIterationCallback(iterationCallback)

  algo.gensSinceImprovement = 0
  algo.bestFitness = 0

  return algo.start()
}
