import igraph as ig
import csv
import gzip

# Load the graph
g = ig.Graph.Read_Ncol("output/hs_core_gene_family_annotations.csv", directed=False)

interpro_info = {}
with gzip.open("output/family_info.csv.gz", "rt") as file:
    reader = csv.reader(file)
    for row in reader:
        interpro_info[row[0]] = row[1:]

heat_core_names = {}
with open("../analyses/1_core_response/input/hs_core_genes_internal_names.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the first row (header)
    for row in reader:
        heat_core_names[row[0]] = row[1]


upregulated_genes = []
downregulated_genes = []
with open("../analyses/1_core_response/input/core_genes.csv", "r") as file:
    reader = csv.reader(file)
    for row in reader:
        if row[4] == "upregulated":
            upregulated_genes.append(row[0])
        else:
            downregulated_genes.append(row[0])

hs_gene_info = {}
with open("../analyses/1_core_response/input/core_genes.csv", "r") as file:
    reader = csv.reader(file)
    for row in reader:
        hs_gene_info[row[0]] = row[1:]


visual_style = {}

color_dict = {"none": "#00000049", "upregulated": "#D55E0099", "downregulated": "#0072B299", "superfamily": "#80008099" }
vertex_colors = []
for vertex in g.vs:
    if vertex["name"] in heat_core_names:
        if heat_core_names[vertex["name"]] in upregulated_genes:
            vertex_colors.append(color_dict["upregulated"])
        elif heat_core_names[vertex["name"]] in downregulated_genes:
            vertex_colors.append(color_dict["downregulated"])
    elif vertex["name"] in interpro_info and interpro_info[vertex["name"]][2] == "HOMOLOGOUS_SUPERFAMILY":
        vertex_colors.append(color_dict["superfamily"])
    else:
        vertex_colors.append(color_dict["none"])

visual_style["vertex_label"] = [ heat_core_names[n] if n in heat_core_names else n for n in g.vs["name"] ]
visual_style["vertex_label"] = [ interpro_info[n][1] if n in interpro_info else n for n in visual_style["vertex_label"] ]
# Introduce line break after 26 characters
visual_style["vertex_label"] = [label[:26] + '\n' + label[26:] if len(label) > 26 else label for label in visual_style["vertex_label"]]
visual_style["vertex_label_size"] = 8
visual_style["vertex_label_dist"] = 1.5
visual_style["vertex_size"] = 7 #[7 if v["name"] in heat_core + validated_proteins else 5 for v in s.vs]
visual_style["vertex_frame_width"] = 0
visual_style["vertex_color"] = vertex_colors # [color_dict["heat-core"] if v["name"] in heat_core else color_dict["validated"] if v["name"] in validated_proteins else color_dict["random"] if v["name"] in random_proteins else color_dict["none"] for v in s.vs]
visual_style["vertex_label_color"] = visual_style["vertex_color"]
visual_style["vertex_shape"] = ["square" if n in interpro_info else "circle" for n in g.vs["name"]]
visual_style["edge_width"] = 0.5
visual_style["edge_color"] = "#00000049" # [color_dict["validated"] if s.vs["name"][e.source] in validated_proteins or s.vs["name"][e.target] in validated_proteins else "#00000011" for e in s.es]
visual_style["margin"] = 60
visual_style["bbox"] = (800, 600)

ig.plot(g, "output/plots/3_hs_core_interpro_families.pdf", **visual_style)

# Print the vertices sorted by number of edges
[x["name"] for x in sorted(g.vs, key=lambda v: v.degree(), reverse=True)]

len(g.connected_components())