# This is a super-slimmed down, minimal version of the HogProf library which I thought was easier
# to re-create than get the full HogProf software running (which is still in early Alpha).
# There are many code snippets taken and simplified from HogProf, mostly the lshbuilder.py, pyhamutils.py, 
# file_utils.py, and hashutils.py files.
import pyham
import numpy as np
from datasketch import WeightedMinHashGenerator

pyham_analysis = pyham.Ham("speciestree.nwk", "oma-hogs.orthoXML", tree_format="newick", use_internal_name=True, species_resolve_mode="OMA")

full_tree = pyham_analysis.taxonomy.tree
taxa_index = {}
for i, n in enumerate(full_tree.traverse()):
    taxa_index[n.name] = i-1

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

  losses = [ taxa_index[n.name]  for n in tp.traverse() if n.lost  ]
  dupl = [ taxa_index[n.name]  for n in tp.traverse() if n.dupl  ]
  presence = [ taxa_index[n.name]  for n in tp.traverse() if n.nbr_genes > 0 ]
  indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )
  for i,event in enumerate(indices):
      if len(indices[event])>0:
          hogindex = np.asarray(indices[event])+i*len(taxa_index)
          hog_matrix_binary[:,hogindex] = 1

  return wmg.minhash(list(hog_matrix_binary.flatten()))

hash1 = hog2hash(178595)
hash2 = hog2hash(666324) #@TODO there is a key error here because "Capsaspora owczarzaki (strain ATCC 30864)" is transformed into "Capsaspora owczarzaki _strain ATCC 30864_" somehow but not in the taxa_index... Not sure why but it's an easy workaround (just transform it when building the taxa_index).
print(hash1.jaccard(hash2))

#@TODO get all HOGs and properly compare all vs all and save it out.

