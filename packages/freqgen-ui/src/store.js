import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
  state: { aminoAcidSeq: '', targetFreqs: new Map() },
  mutations: {
    setAminoAcidSeq(state, seq) {
      state.aminoAcidSeq = seq
    },
    setTargetFreqs(state, targetFreqs) {
      state.targetFreqs = targetFreqs
    },
  },
  actions: {},
})
