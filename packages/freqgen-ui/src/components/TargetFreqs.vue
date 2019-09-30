<template>
  <div>
    <b-card-text>
      <b-container>
        <p>
          To get started, we need to know what the target frequencies are. You
          can either use a FASTA file containing DNA sequences as to calculate
          frequencies from, use a Freqgen-generated YAML file with precalculated
          frequencies, or manually enter frequencies.
        </p>
      </b-container>
      <b-container>
        <b-card no-body>
          <b-tabs pills card vertical>
            <b-tab title="FASTA" active>
              <TargetFreqsFasta />
            </b-tab>
            <b-tab title="YAML">
              <TargetFreqsYaml />
            </b-tab>
          </b-tabs>
        </b-card>
      </b-container>
    </b-card-text>
  </div>
</template>

<script>
import TargetFreqsFasta from './TargetFreqsFasta.vue'
import TargetFreqsYaml from './TargetFreqsYaml.vue'
export default {
  name: 'TargetFreqs',
  data() {
    return {
      text: '',
    }
  },
  components: { TargetFreqsFasta, TargetFreqsYaml },
}
</script>

<style></style>
