module.exports = function(map1, map2) {
  let result = new Map()
  for (let [k, v] of map2.entries()) {
    map1.get(k) ? result.set(k, map1.get(k) + v) : result.set(k, v)
  }
  for (let [k, v] of map1) {
    if (!result.has(k)) {
      result.set(k, v)
    }
  }
  return result
}
