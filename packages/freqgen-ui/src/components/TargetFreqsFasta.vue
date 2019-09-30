<template>
  <b-card-text>
    <b-form-textarea
      id="textarea"
      v-model="text"
      placeholder=">Enter FASTA..."
      rows="3"
      max-rows="6"
      :disabled="Boolean(file.length)"
    ></b-form-textarea>
    <p class="text-center my-3">or</p>
    <!-- <b-container> -->
    <b-row>
      <b-col cols="10">
        <b-form-file
          v-model="file"
          :state="Boolean(file.length) ? true : null"
          placeholder="Choose a file or drop it here..."
          drop-placeholder="Drop file here..."
          accept=".fna, .fasta, .fa"
          multiple
        ></b-form-file>
      </b-col>
      <b-col cols="2">
        <b-button
          variant="outline-secondary"
          block
          :disabled="!Boolean(file.length)"
          @click="file = []"
        >Clear file{{ file.length > 1 ? 's' : '' }}</b-button>
      </b-col>
    </b-row>
    <!-- </b-container> -->

    <p class="text-center my-3">
      then, select
      <i>k</i>-mer targets
    </p>

    <div class="text-center">
      <b-form-checkbox-group v-model="selected" :options="options" switches name="buttons-1"></b-form-checkbox-group>
    </div>
    <div class="text-center my-3">
      <b-button variant="outline-secondary" :disabled="!Boolean(text)">Next →</b-button>
    </div>
  </b-card-text>
</template>

<script>
export default {
  name: 'TargetFreqsFasta',
  data() {
    return {
      file: [],
      options: [1, 2, 3, 'codons', 4, 5, 6, 7, 8, 9],
      selected: [],
    }
  },
}
</script>

<style>
</style>