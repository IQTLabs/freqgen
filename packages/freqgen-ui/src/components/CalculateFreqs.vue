<template>
  <b-tabs position="is-centered" class="box" type="is-toggle" :animated="false">
    <b-tab-item label="FASTA" icon="file-document">
      <PasteOrUpload />
    </b-tab-item>
    <b-tab-item label="YAML" icon="file-chart">
      <PasteOrUpload></PasteOrUpload>
    </b-tab-item>
  </b-tabs>
</template>

<script>
import PasteOrUpload from './PasteOrUpload.vue'
export default {
  name: 'CalculateFreqs',
  components: {
    PasteOrUpload,
  },
}
</script>
