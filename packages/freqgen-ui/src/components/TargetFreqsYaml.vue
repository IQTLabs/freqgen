<template>
  <b-card-text>
    <b-form-textarea
      id="textarea"
      v-model="text"
      placeholder="Paste YAML here..."
      rows="3"
      max-rows="6"
      :disabled="Boolean(file)"
      :state="textareaState"
      @input="loadText"
    ></b-form-textarea>
    <p class="text-center my-3">or</p>
    <b-row>
      <b-col cols="10">
        <b-form-file
          v-model="file"
          :state="fileState"
          placeholder="Choose a file or drop it here..."
          drop-placeholder="Drop file here..."
          accept=".yml, .yaml"
          @change="loadFile"
        ></b-form-file>
      </b-col>
      <b-col cols="2">
        <b-button
          variant="outline-secondary"
          block
          :disabled="!Boolean(file)"
          @click="file = null; fileState = null"
        >Clear file</b-button>
      </b-col>
    </b-row>
  </b-card-text>
</template>

<script>
const yaml = require('js-yaml')
export default {
  name: 'TargetFreqsYaml',
  data() {
    return {
      file: null,
      yamlObj: undefined,
      text: '',
      fileState: null,
    }
  },
  methods: {
    loadText() {
      try {
        this.yamlObj = yaml.safeLoad(this.text)
      } catch (error) {
        this.yamlObj = undefined
      }
    },
    loadFile(ev) {
      let reader = new FileReader()
      reader.readAsText(ev.target.files[0])
      reader.onload = e => {
        // if things worked, update the object
        try {
          this.yamlObj = yaml.safeLoad(e.target.result)
          console.log(this.yamlObj)

          this.fileState = true
        } catch (error) {
          // otherwise, create an error message and reset things to try again
          this.yamlObj = undefined
          this.fileState = false
          this.$bvModal.msgBoxOk('Unable to parse YAML file.').then(() => {
            this.file = null
            this.fileState = null
          })
        }
      }
    },
  },
  computed: {
    // decide whether to make the textarea outline red, green, or blank
    textareaState() {
      if (typeof this.yamlObj == 'object' && this.text.length > 0) {
        return true
      } else if (typeof this.yamlObj != 'object' && this.text.length != 0) {
        return false
      } else {
        return null
      }
    },
  },
}
</script>

<style>
</style>