paths = read.csv("output/shortest_paths.csv")
paths$length = paths$length - 1 # remove the starting node
print(t.test(paths$length ~ paths$set))

print("Proportion of direct interactions in random set")
print(sum(paths[paths$set == "random",]$length == 1) / sum(paths$set == "random"))

print("Proportion of direct interactions in heat-core set")
print(sum(paths[paths$set == "heat-core",]$length == 1) / sum(paths$set == "heat-core"))