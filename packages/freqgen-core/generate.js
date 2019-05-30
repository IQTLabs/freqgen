// const Genetic = require('genetic-js')
// const yaml = require('js-yaml')
// const fs = require('fs')

// configure warnings a la parcel/packages/core/test-utils/src/utils.js#L15 @ 30624f7
const chalk = require('chalk')
const warning = chalk.keyword('orange')
console.warn = (...args) => {
  console.error(warning(...args))
}

const utilities = require('./utilities')

const DNA_BASES = new Set(['A', 'T', 'G', 'C'])

module.exports = function generate(
  targetFreqs,
  targetAminoAcidSeq,
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

  // although we want a Map, don't break if given an Object
  if (!(targetFreqs instanceof Map)) {
    targetFreqs = new Map(Object.entries(targetFreqs))
  }

  utilities.validateKmerFrequencyMap(targetFreqs)
}

// codons_for_aa = yaml.load(fs.readFileSync('data/codons_for_aa.yaml', 'utf8'))
// console.log(codons_for_aa)

// var genetic = Genetic.create()
// genetic.optimize = Genetic.Optimize.Maximize
// genetic.select1 = Genetic.Select1.Tournament2
// genetic.select2 = Genetic.Select2.FittestRandom
//
// genetic.seed = function () {
//   var a = []
//   // create coefficients for polynomial with values between (-0.5, 0.5)
//
//   var i
//   for (i = 0; i < this.userData.foods.length; ++i) {
//     a.push(Math.random() < 0.5)
//   }
//
//   return a
// }
//
// genetic.mutate = function (entity) {
//   var idx = Math.floor(Math.random() * entity.length)
//   entity[idx] = !entity[idx]
//
//   return entity
// }
//
// genetic.crossover = function (mother, father) {
//   var idx = Math.floor(Math.random() * mother.length)
//
//   var m1 = mother.slice(0, idx),
//     m2 = mother.slice(idx),
//     f1 = father.slice(0, idx),
//     f2 = father.slice(idx)
//
//   var son = m1.concat(f2)
//   var daughter = f1.concat(m2)
//
//   return [son, daughter]
// }
//
// genetic.fitness = function (entity) {
//   let putative = {
//     calories: 0,
//     carbs: 0,
//     fat: 0,
//     protein: 0
//   }
//
//   var sumSqErr = 0
//
//   // calculate the total macros of the solution
//   for (let i = 0; i < entity.length; i++) {
//     if (entity[i]) {
//       for (let macro of ['calories', 'protein', 'fat', 'carbs']) {
//         putative[macro] += this.userData.foods[i][macro]
//       }
//     }
//   }
//
//   for (let macro of ['calories', 'protein', 'fat', 'carbs']) {
//     sumSqErr += (this.userData.target[macro] - putative[macro]) ** 2
//   }
//
//   return 1 / Math.sqrt(sumSqErr)
// }
//
// genetic.notification = function (pop, gen, stats, isFinished) {
//   if (isFinished) {
//     // console.log(pop, gen, stats);
//     console.log(solutionToFoods(pop[0], this.userData.foods))
//   }
// }
//
// genetic.generation = function (pop, gen, stats) {
//   if (stats.maximum > this.best_fitness) {
//     this.best_fitness = stats.maximum
//     this.gens_since_improvement = 0
//   } else {
//     this.gens_since_improvement += 1
//   }
//   if (this.gens_since_improvement > 50) {
//     return false
//   }
// }
//
// genetic.best_fitness = 0
// genetic.gens_since_improvement = 0
//
// // console.log(genetic);
// var config = {
//   'iterations': 500,
//   'size': foods.length > 500 ? 500 : foods.length * 10,
//   'crossover': 0.8,
//   'mutation': 0.2,
//   'skip': 10
// }
//
// var userData = {
//   // genomeSize: 10,
//   foods: foods,
//   target: target
// }
//
// genetic.evolve(config, userData)
