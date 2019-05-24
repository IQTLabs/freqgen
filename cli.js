const fs = require('fs')
const yaml = require('js-yaml')
const freqgen = require('./index.js')
const ora = require('ora')
const Fasta = require('biojs-io-fasta')
const _ = require('lodash')

var program = require('commander')

function commaSeparatedIntList (value, dummyPrevious) {
  return [...new Set(value.split(',').map(x => parseInt(x)))] // deduplicate the list
}

program
  .command('featurize [files...]')
  .description('Featurize one or more FASTA files')
  .option('-k, --k-mers <int>', 'Comma separated list of k values to featurize. Example: -k 1,2,3', commaSeparatedIntList)
  .option('-c, --codons', 'Whether to featurize codons')
  .option('-o, --output <file>', 'The output YAML file')
  .action(function (files, options) {
    const spinner = ora(`Parsing ${files.length} files...`).start()

    if (options.kMers == null) {
      spinner.fail('No k-mers or codons specified to featurize. Provide at least one k value after -k or use the -c flag to featurize codons.')
      return
    }

    let seqs = files.map(file => Fasta.parse(fs.readFileSync(file, 'utf8')).map(obj => obj.seq))
    seqs = _.flatten(seqs)

    let totalKmerCounts = {}

    for (let k of options.kMers) {
      spinner.text = `Counting ${k}-mers...`
      for (let seq of seqs) {
        let counts = freqgen.kmerCounts(freqgen.kmers(seq, k))
        totalKmerCounts[k] = _.mergeWith({}, totalKmerCounts[k],
          counts, function (objValue, srcValue) {
            return _.isNumber(objValue) ? objValue + srcValue : srcValue
          })
      }
    }

    if (!(options.codons == null)) {
      spinner.text = `Counting codons...`
      for (let seq of seqs) {
        let counts = freqgen.kmerCounts(freqgen.kmers(seq, 3, overlap = false))
        totalKmerCounts['codons'] = _.mergeWith({}, totalKmerCounts['codons'],
          counts, function (objValue, srcValue) {
            return _.isNumber(objValue) ? objValue + srcValue : srcValue
          })
      }
    }

    let kmerFrequencies = _.mapValues(totalKmerCounts, kmerCounts => freqgen.kmerFrequencies(kmerCounts))

    // console.log(files, options.kMers || 0, !(options.codons == null))

    spinner.succeed('Done!')
    if (options.output == null) {
      console.log(yaml.safeDump(kmerFrequencies))
    } else {
      fs.writeFileSync(options.output, yaml.safeDump(kmerFrequencies))
    }
  })

program.parse(process.argv)

// if (program.kMers !== undefined) console.log(program.kMers)
//
// if (program.codons) {
//   console.log('featurizing codons!')
// }

// console.log(yaml.safeDump(f reqgen.kmers()))
