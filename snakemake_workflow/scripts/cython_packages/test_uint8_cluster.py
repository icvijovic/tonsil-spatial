import numpy as np
from MST_SINGLE_UINT8 import *

from scipy.cluster.hierarchy import cut_tree
import scipy.spatial.distance as ssd

D = np.arange(100, dtype = np.uint8).reshape((10,10))
D = D**2
for i in range(10):
	for j in range(i):
		D[i,i] = 0
		D[j,i] = D[i,j]

print(D)

C = ssd.squareform(D)
print(C.dtype)
CLONE_DISTANCE_CUTOFF = 2
n = D.shape[0]

def clone_ids(condensed_distance_matrix, n, cutoff = CLONE_DISTANCE_CUTOFF):
        if len(condensed_distance_matrix) == 0:
                return [0]
        else:
            #print("current memory usage:", getCurrentMemoryUsage())
            Z = mst_single_linkage(C,n)
            clones = cut_tree(Z, n_clusters=None, height=cutoff)
            clones = np.asarray([x[0] for x in clones])
            return clones

print(clone_ids(C,n))
