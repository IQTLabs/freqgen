var distance = module.exports

const utilities = require('./utilities')
distance.euclidean = function(map1, map2) {
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

distance.cosine = function(map1, map2) {
  return (
    utilities.mapDotProduct(map1, map2) /
    (distance.euclidean(map1) * distance.euclidean(map2))
  )
}
