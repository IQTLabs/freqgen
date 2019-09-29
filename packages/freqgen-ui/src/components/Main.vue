<template>
  <b-container>
    <b-card no-body>
      <b-tabs card justified>
        <b-tab title="1. Target Frequencies" active>
          <TargetFreqs></TargetFreqs>
        </b-tab>
        <b-tab title="2. Target Amino Acid Sequence" disabled>
          <b-card-text>Tab Contents 2</b-card-text>
        </b-tab>
        <b-tab title="3. Generate DNA Sequence" disabled>
          <b-card-text>Tab Contents 2</b-card-text>
        </b-tab>
      </b-tabs>
    </b-card>
  </b-container>
</template>

<script>
import TargetFreqs from './TargetFreqs.vue'
export default { name: 'Main', components: { TargetFreqs } }
</script>

<style>
</style>