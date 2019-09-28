<template>
  <b-steps>
    <b-step-item label="Determine Frequencies" icon="chart-bar">
      <CalculateFreqs />
    </b-step-item>
    <b-step-item label="Protein Target" icon="target">
      <div class="box">
        <h2 class="title has-text-centered is-4">Paste FASTA</h2>
        <PasteOrUpload filetype-name="FASTA" acceptable-filetypes=".fna,.fa,.fasta"></PasteOrUpload>
      </div>
    </b-step-item>
    <b-step-item label="Generate DNA" icon="dna">
      <div class="box">
        <h1 class="title">DNA goes here.</h1>
      </div>
    </b-step-item>
  </b-steps>
</template>

<script>
import CalculateFreqs from './CalculateFreqs.vue'
import PasteOrUpload from './PasteOrUpload.vue'
export default {
  name: 'Steps',
  components: {
    CalculateFreqs,
    PasteOrUpload,
  },
}
</script>

<style>
.b-steps .steps .step-items .step-item {
  flex-basis: 1em !important;
}
</style>
