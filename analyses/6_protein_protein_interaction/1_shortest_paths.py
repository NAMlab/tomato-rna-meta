# In this script we will calculate the shortest path between each gene in the heat response core set as well as a random
# baseline of 200 genes from ITAG4.1 to any of the validated genes. 
import igraph as ig
import csv

# Load the graph
g = ig.Graph.Read_Ncol("intermediary/string_edges.lfs.txt", directed=False)
g.simplify() # there are apparently some self-loops or duplicates in the graph
g.ecount()

# Load our relevant gene sets
heat_core = []
heat_core_names = {}
with open("../1_core_response/input/hs_core_genes_internal_names.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the first row (header)
    for row in reader:
        heat_core.append(row[0])
        heat_core_names[row[0]] = row[1]


upregulated_proteins = []
downregulated_proteins = []
with open("../1_core_response/input/core_genes.csv", "r") as file:
    reader = csv.reader(file)
    for row in reader:
        if row[4] == "upregulated":
            upregulated_proteins.append(row[0])
        else:
            downregulated_proteins.append(row[0])


validated_proteins = []
with open("input/validated_heat_proteins.txt", "r") as file:
    for line in file:
        validated_proteins.append(line.strip())

random_proteins = []
with open("input/random_proteins.txt", "r") as file:
    for line in file:
        random_proteins.append(line.strip())

# Remove any entries that are not in the graph
validated_proteins = [protein for protein in validated_proteins if protein in g.vs["name"]]

# Get the direct neighbors of the proteins in validated_proteins
direct_neighbors = set()
for protein in validated_proteins:
    neighbors = g.neighbors(protein)
    direct_neighbors.update(neighbors)

# Calculate the proportion of entries in heat_core and in the entire proteome that are in the direct neighborhood
direct_neighbor_names = [g.vs[index]["name"] for index in direct_neighbors]
proportion = len(set(heat_core).intersection(direct_neighbor_names)) / len(heat_core)
print("Proportion of entries in heat_core in the direct neighborhood of validated_proteins:", proportion)
print("Proportion of entries in heat_core in the direct neighborhood of validated_proteins:", len(direct_neighbors) / 34688.0)

# Now remove any entries that are not in the graph from the other sets
heat_core = [gene for gene in heat_core if gene in g.vs["name"]]
random_proteins = [gene for gene in random_proteins if gene in g.vs["name"]]

# This is for plotting a subnetwork later, we'll add the intermediate genes as well
subnetwork = heat_core
# Calculate and save the shortest path of each core and random gene to any of the validated genes
with open("output/shortest_paths.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow(["set", "source", "target", "paths", "length"])
    for source in heat_core:
        print(source)
        current_shortest_path = None
        current_target = None
        current_paths = None
        for target in validated_proteins:
            print(" " + target)
            paths = g.get_shortest_paths(source, target)
            if paths:
                if not current_shortest_path or len(paths[0]) < current_shortest_path:
                    current_shortest_path = len(paths[0])
                    current_target = target
                    current_paths = [g.vs[x]["name"] for x in paths]
        subnetwork = subnetwork + [g for path in current_paths for g in path] # flatten the list
        writer.writerow(["heat-core", source, current_target, current_paths, current_shortest_path or -1])
    for source in random_proteins:
        print(source)
        current_shortest_path = None
        current_target = None
        current_paths = None
        for target in validated_proteins:
            print(" " + target)
            paths = g.get_shortest_paths(source, target)
            if paths:
                if not current_shortest_path or len(paths[0]) < current_shortest_path:
                    current_shortest_path = len(paths[0])
                    current_target = target
                    current_paths = [g.vs[x]["name"] for x in paths]
        writer.writerow(["random", source, current_target, current_paths, current_shortest_path or -1])

subnetwork = list(set(subnetwork)) # remove duplicates

s = g.subgraph(subnetwork)
visual_style = {}
layout = s.layout("kk")
visual_style["layout"] = layout
color_dict = {"none": "#00000049", "upregulated": "#D55E0099", "downregulated": "#0072B299", "validated": "#009E7399", "random": "#F0E442aa"}
vertex_colors = []
for vertex in s.vs:
    if vertex["name"] in heat_core:
        if heat_core_names[vertex["name"]] in upregulated_proteins:
            vertex_colors.append(color_dict["upregulated"])
        elif heat_core_names[vertex["name"]] in downregulated_proteins:
            vertex_colors.append(color_dict["downregulated"])
    elif vertex["name"] in validated_proteins:
        vertex_colors.append(color_dict["validated"])
    elif vertex["name"] in random_proteins:
        vertex_colors.append(color_dict["random"])
    else:
        vertex_colors.append(color_dict["none"])

visual_style["vertex_size"] = [7 if v["name"] in heat_core + validated_proteins else 5 for v in s.vs]
visual_style["vertex_frame_width"] = 0
visual_style["vertex_color"] = vertex_colors # [color_dict["heat-core"] if v["name"] in heat_core else color_dict["validated"] if v["name"] in validated_proteins else color_dict["random"] if v["name"] in random_proteins else color_dict["none"] for v in s.vs]
visual_style["edge_width"] = 0.5
visual_style["edge_color"] = "#00000049" # [color_dict["validated"] if s.vs["name"][e.source] in validated_proteins or s.vs["name"][e.target] in validated_proteins else "#00000011" for e in s.es]
visual_style["vertex_label"] = [ heat_core_names[n] if n in heat_core_names else n for n in s.vs["name"] ]
visual_style["vertex_label_color"] = visual_style["vertex_color"]
visual_style["vertex_label_size"] = 8
visual_style["vertex_label_dist"] = 1.5
visual_style["bbox"] = (600, 400)

ig.plot(s, "output/plots/6.1_shortest_path_subgraph.pdf", **visual_style)
