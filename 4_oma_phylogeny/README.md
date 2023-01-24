The hypothesis is that the genes we found in our core set are "ancient", i.e. they should exist across all kinds of species, at least across a wider range than a randomly drawn sample of the tomato genome.
This does not necessarily mean that the orthologues have the same function in all these species -- we would hope that is the case but so far we have no method of testing this second hypothesis.
For the first one, we are using the OMA Database. For each core response gene + a random sample of 500 genes from the tomato genome, we are

1. Querying OMA with the protein sequence and getting the genes it finds for it in tomato
2. Finding out which HOGs these genes belong to
3. Determining what the deepest level of that HOG is, i.e. how much of the phylogenetic tree it spans.

`clades_hierarchy.txt` was created manually by inspecting the species phylogeny provided by OMA at https://omabrowser.org/oma/current/ .
