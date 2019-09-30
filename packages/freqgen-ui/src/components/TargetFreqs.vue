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
                  <b-form-checkbox-group
                    v-model="selected"
                    :options="options"
                    switches
                    name="buttons-1"
                  ></b-form-checkbox-group>
                </div>
                <div class="text-center my-3">
                  <b-button variant="outline-secondary" :disabled="!Boolean(text)">Next →</b-button>
                </div>
              </b-card-text>
            </b-tab>
            <b-tab title="YAML">
              <b-card-text>
                <b-form-textarea
                  id="textarea"
                  v-model="text"
                  placeholder="Paste YAML here..."
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
                      accept=".yml, .yaml"
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
              </b-card-text>
            </b-tab>
          </b-tabs>
        </b-card>
      </b-container>
    </b-card-text>
  </div>
</template>

<script>
export default {
  name: 'TargetFreqs',
  data() {
    return {
      file: [],
      options: [1, 2, 3, 'codons', 4, 5, 6, 7, 8, 9],
      selected: [],
      text: '',
    }
  },
}
</script>

<style></style>
