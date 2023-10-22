# This is a super-slimmed down, minimal version of the HogProf library which I thought was easier
# to re-create than get the full HogProf software running (which is still in early Alpha).
# There are many code snippets taken and simplified from HogProf, mostly the lshbuilder.py, pyhamutils.py, 
# file_utils.py, and hashutils.py files.
import pyham
import numpy as np
import pandas as pd
from datasketch import WeightedMinHashGenerator

STREPTOPHYTA_ONLY = False

# The oma-hogs.orthoXML is too big for this repo, you can download the latest versions of both files from https://omabrowser.org/oma/current/
pyham_analysis = pyham.Ham("input/speciestree.nwk", "input/oma-hogs.orthoXML", tree_format="newick", use_internal_name=True, species_resolve_mode="OMA")

full_tree = pyham_analysis.taxonomy.tree
taxa_index = {}
# Search term here limits analysis to specific clade, e.g. Eukaryota or Streptophyta
tree_traverser = enumerate(full_tree.search_nodes(name="Streptophyta")[0].traverse()) if STREPTOPHYTA_ONLY else enumerate(full_tree.traverse())
for i, n in tree_traverser:
  # For some reason newick chars are replaced by "_" automatically in the tree
  # profile but not here, so we have to do it manually to avoid a key error later.
  name = n.name.replace("(", "_").replace(")", "_").replace(",", "_").replace(":", "_")
  taxa_index[name] = i-1

wmg = WeightedMinHashGenerator(3*len(taxa_index), sample_size = 256 , seed=1)

def hog2hash(hog_id):
  # The hog_id needs to be the pure numerical int ID, not "HOG:XYZ"!
  # you can get this ID through the OMA RESt interface HOG/read
  hog = pyham_analysis.get_hog_by_id(hog_id) 
   
  # This treemap is now what is actually used for phylogenetic profiling, i.e. it's the actual profile.
  tp = pyham_analysis.create_tree_profile(hog=hog).treemap

  # This is going to be a matrix of 0s and 1s marking whether presence, loss, and/or duplication has
  # happened at each node in the tree for this hog.
  # In HogProf there is also a weighted version but tee-weights are all 1 by default so we can 
  # stay with the binary matrix. It is initialized to all 0s and then later 1s are filled in where the
  # corresponding events happened. It has a size of 3*the amount of taxa, the first third marks in which
  # the HOG is present, the second one where losses occured and the third one duplucations. It is the
  # vector which will then be converted to a minHash.
  hog_matrix_binary = np.zeros((1, 3*len(taxa_index)))

  losses = [ taxa_index[n.name]  for n in tp.traverse() if n.lost and n.name in taxa_index ]
  dupl = [ taxa_index[n.name]  for n in tp.traverse() if n.dupl and n.name in taxa_index ]
  presence = [ taxa_index[n.name]  for n in tp.traverse() if n.nbr_genes > 0 and n.name in taxa_index ]
  indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
  for i,event in enumerate(indices):
      if len(indices[event])>0:
          hogindex = np.asarray(indices[event])+i*len(taxa_index)
          hog_matrix_binary[:,hogindex] = 1

  return wmg.minhash(list(hog_matrix_binary.flatten()))

# TEST CASES
# hash1 = hog2hash(178595)
# hash2 = hog2hash(666324)
# print(hash1.jaccard(hash2))

# Read HS core hogs and perform all v all jaccard distances
df = pd.read_csv("output/combined_deepest_levels.csv")
df = df[(df.set == "heat") & (df.hog_id.notnull())]
df.hog_id = df.hog_id.astype(int)

# Save minhashes once for each HOG
df["minhash"] = df.apply(lambda r: hog2hash(r["hog_id"]), axis=1) 

jaccard_matrix = np.zeros((len(df), len(df)))
for i, hog1 in enumerate(df.hog_id):
  for j, hog2 in enumerate(df.hog_id):
    if i < j:
      jaccard_matrix[i,j] = df.minhash.values[i].jaccard(df.minhash.values[j])

jaccard_matrix = jaccard_matrix + jaccard_matrix.T
np.fill_diagonal(jaccard_matrix, 1)

df_out = pd.DataFrame(jaccard_matrix, columns = df.target)
df_out.insert(0, "target", list(df.target))

outfile_name = "output/heat_hogprof_jaccard_matrix.csv" if STREPTOPHYTA_ONLY else "output/heat_hogprof_jaccard_matrix_all_taxa.csv"
df_out.to_csv(outfile_name, index=False)
