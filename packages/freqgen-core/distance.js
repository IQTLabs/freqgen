module.exports.euclidean = function(map1, map2) {
  let seenKeys = new Set()

  let sum = 0

  if (map2 === undefined) {
    map2 = new Map()
  }

  // first pass through map1
  for (let entry of map1.entries()) {
    seenKeys.add(entry[0])
    sum +=
      (entry[1] - (map2.get(entry[0]) == null ? 0 : map2.get(entry[0]))) ** 2
  }

  // then pass through
  for (let entry of map2.entries()) {
    if (!seenKeys.has(entry[0])) {
      sum += entry[1] ** 2
    }
  }

  return Math.sqrt(sum)
}

module.exports.cosine = function(map1, map2) {
  let dotProduct = 0
  // first pass through map1
  for (let entry of map1.entries()) {
    dotProduct +=
      entry[1] * (map2.get(entry[0]) == null ? 0 : map2.get(entry[0]))
  }

  return dotProduct / (this.euclidean(map1) * this.euclidean(map2))
}
