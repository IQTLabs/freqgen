<template>
  <b-tab title="1. Target Frequencies" active>
    <b-card-text>
      <b-container>
        <p>To get started, we need to know what the target frequencies are. You can either use a FASTA file containing DNA sequences as to calculate frequencies from, use a Freqgen-generated YAML file with precalculated frequencies, or manually enter frequencies.</p>
      </b-container>
      <b-container>
        <b-card no-body>
          <b-tabs pills card vertical>
            <b-tab title="FASTA" active>
              <b-card-text>
                <b-form-textarea
                  id="textarea"
                  v-model="text"
                  placeholder="Enter something..."
                  rows="3"
                  max-rows="6"
                ></b-form-textarea>
                <p class="text-center my-3">or</p>

                <b-form-file
                  v-model="file"
                  :state="Boolean(file) ? true : null"
                  placeholder="Choose a file or drop it here..."
                  drop-placeholder="Drop file here..."
                ></b-form-file>
              </b-card-text>
            </b-tab>
            <b-tab title="YAML">
              <b-card-text>Tab Contents 2</b-card-text>
            </b-tab>
          </b-tabs>
        </b-card>
      </b-container>
    </b-card-text>
  </b-tab>
</template>

<script>
export default {
  name: 'TargetFreqs',
  data() {
    return {
      file: null,
    }
  },
}
</script>

<style>
</style>