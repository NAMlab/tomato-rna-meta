This part of the analysis looks into evolutionary history and co-evolution among the stress core genes.
We first find the Hierarchical Orthologous Groups (HOGs) in OMA (https://omabrowser.org/oma/home/) for each of the core gene sets and the random baseline (sequences in the output of the analysis 1):

1. Querying OMA with the protein sequence and getting the genes it finds for it in tomato
2. Finding out which HOGs these genes belong to
3. Determining what the deepest level of that HOG is, i.e. how much of the phylogenetic tree it spans.

The hypothesis is that the genes we found in our core set are "ancient", i.e. they should exist across all kinds of species, at least across a wider range than a randomly drawn sample of the tomato genome.
We plot the results of this analysis in the plot 4.3.

Next we look into the co-evolution (and co-expression) of our heat stress core genes using the HogProf algorithm (https://doi.org/10.1371/journal.pcbi.1007553). The results are visualized in plot 4.6.