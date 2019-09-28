<template>
  <div>
    <b-field position="is-centered">
      <b-input
        type="textarea"
        :disabled="textareaDisabled"
        :placeholder="'Paste ' + filetypeName + ' here'"
      ></b-input>
    </b-field>
    <h4 class="title is-4 has-text-centered">or</h4>
    <div class="has-text-centered">
      <b-field position="is-centered">
        <b-upload v-model="dropFiles" multiple drag-drop expanded :accept="acceptableFiletypes">
          <section class="section">
            <div class="content has-text-centered">
              <p>
                <b-icon icon="upload" size="is-large"></b-icon>
              </p>
              <p>Drop your {{filetypeName}} files here or click to upload</p>
            </div>
          </section>
        </b-upload>
      </b-field>
      <div class="tags is-centered">
        <span v-for="(file, index) in dropFiles" :key="index" class="tag is-primary">
          {{file.name}}
          <button class="delete is-small" type="button" @click="deleteDropFile(index)"></button>
        </span>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'PasteOrUpload',
  props: ['filetypeName', 'acceptableFiletypes'],
  data() {
    return {
      dropFiles: [],
    }
  },
  methods: {
    deleteDropFile(index) {
      this.dropFiles.splice(index, 1)
    },
  },
  computed: {
    textareaDisabled() {
      return this.dropFiles.length > 0
    },
  },
}
</script>

<style>
</style>