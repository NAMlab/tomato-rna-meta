paths = read.csv("output/shortest_paths.csv")
paths$length = paths$length - 1 # remove the starting node
print(t.test(paths$length ~ paths$set))