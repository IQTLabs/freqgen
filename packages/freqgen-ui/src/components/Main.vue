<template>
  <div class="mb-5">
    <b-card no-body>
      <b-tabs card justified>
        <b-tab title="1. Target Frequencies" active>
          <TargetFreqs></TargetFreqs>
        </b-tab>
        <b-tab title="2. Target Amino Acid Sequence">
          <TargetAminoAcid></TargetAminoAcid>
        </b-tab>
        <b-tab title="3. Generate DNA Sequence" disabled>
          <b-card-text>Tab Contents 2</b-card-text>
        </b-tab>
      </b-tabs>
    </b-card>
  </div>
</template>

<script>
import TargetFreqs from './TargetFreqs.vue'
import TargetAminoAcid from './TargetAminoAcid.vue'
export default { name: 'Main', components: { TargetFreqs, TargetAminoAcid } }
</script>
