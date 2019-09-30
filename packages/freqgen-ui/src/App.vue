<template>
  <div id="app">
    <Navbar></Navbar>
    <b-container>
      <b-row>
        <b-col>
          <b-jumbotron
            header="DNA optimization made easy."
            bg-variant="transparent"
            lead="Generate DNA sequences with target k-mer and codon usage frequencies."
            class="mb-0"
          ></b-jumbotron>
        </b-col>
      </b-row>
      <b-row>
        <b-col>
          <Main></Main>
        </b-col>
      </b-row>
    </b-container>
  </div>
</template>

<script>
import Main from './components/Main.vue'
import Navbar from './components/Navbar.vue'
export default {
  name: 'app',
  components: {
    Main,
    Navbar,
  },
}
</script>
