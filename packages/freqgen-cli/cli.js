const fs = require('fs')
const yaml = require('./yaml')
const freqgen = require('@freqgen/core')
const ora = require('ora')
const Fasta = require('biojs-io-fasta')

var program = require('commander')

function commaSeparatedIntList(value) {
  return [...new Set(value.split(',').map(x => parseInt(x)))] // deduplicate the list
}

program
  .command('featurize [files...]')
  .description('Featurize one or more FASTA files')
  .option(
    '-k, --k-mers <int>',
    'Comma separated list of k values to featurize. Example: -k 1,2,3',
    commaSeparatedIntList
  )
  .option('-c, --codons', 'Whether to featurize codons')
  .option('-o, --output <file>', 'The output YAML file')
  .action(function(files, options) {
    // start up a pretty spinner (and correctly use file vs files!)
    const spinner = ora(
      `Parsing ${files.length} file${files.length > 1 ? 's' : ''}...`
    ).start()

    if (options.kMers == null) {
      spinner.fail(
        'No k-mers or codons specified to featurize. Provide at least one k value after -k or use the -c flag to featurize codons.'
      )
      return
    }

    // parse the FASTA files into a flat list. Ex: ["ATGC...", "GTCAA...", ...]
    let seqs = files
      .map(file =>
        Fasta.parse(fs.readFileSync(file, 'utf8')).map(obj => obj.seq)
      )
      .flat()

    let totalKmerCounts = new Map() // will map k values to maps with k-mers and counts. Ex: {1: {"A": n, "T": n}...}
    for (let k of options.kMers) {
      spinner.text = `Counting ${k}-mers...`
      for (let seq of seqs) {
        let counts = freqgen.kmerCounts(freqgen.kmers(seq, k, true))

        if (totalKmerCounts.get(k)) {
          totalKmerCounts.set(
            k,
            freqgen.utilities.addMaps(totalKmerCounts.get(k), counts)
          )
        } else {
          totalKmerCounts.set(k, counts)
        }
      }
    }

    // if codon featurization is requested, count the codons of every seq
    if (!(options.codons == null)) {
      spinner.text = 'Counting codons...'
      for (let seq of seqs) {
        let counts = freqgen.kmerCounts(freqgen.kmers(seq, 3, false))

        if (totalKmerCounts.get('codons')) {
          totalKmerCounts.set(
            'codons',
            freqgen.utilities.addMaps(totalKmerCounts.get('codons'), counts)
          )
        } else {
          totalKmerCounts.set('codons', counts)
        }
      }
    }

    let kmerFrequencies = new Map()

    for (let entry of totalKmerCounts.entries()) {
      kmerFrequencies.set(entry[0], freqgen.kmerFrequencies(entry[1]))
    }

    spinner.succeed(
      `Done featurizing ${files.length} file${files.length > 1 ? 's' : ''}! ${
        options.output == null
          ? ''
          : 'Output written to ' + options.output + '.'
      }`
    )

    // either write to a file or print it out
    if (options.output == null) {
      console.log(yaml(kmerFrequencies))
    } else {
      fs.writeFileSync(options.output, yaml(kmerFrequencies))
    }
  })

program.parse(process.argv)
