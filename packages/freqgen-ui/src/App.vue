<template>
  <div id="app">
    <link
      rel="stylesheet"
      href="https://cdn.materialdesignicons.com/2.5.94/css/materialdesignicons.min.css"
    />
    <NavBar />
    <div class="container">
      <Steps />
    </div>
  </div>
</template>

<script>
import NavBar from './components/NavBar.vue'
import Steps from './components/Steps.vue'

export default {
  name: 'app',
  components: {
    NavBar,
    Steps,
  },
}
</script>
<style >
body {
  text-align: center;
}
</style>