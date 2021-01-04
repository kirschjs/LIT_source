from bridge import *
import gc
import itertools

tree = scipy.spatial.cKDTree(mynumbers)  #build k-dimensional trie
pairs = tree.query_pairs(r)  #find all pairs closer than radius: r
neighbors = {}  #dictionary of neighbors
for i, j in pairs:  #iterate over all pairs
    neighbors.get(i, set()).add(j)
    neighbors.get(j, set()).add(i)
keep = []
discard = set(
)  #a list would work, but I use a set for fast member testing with `in`
for node in mynumbers:
    if node not in discard:  #if node already in discard set: skip
        keep.append(node)  #add node to keep list
        discard.update(neighbors.get(
            node, set()))  #add node's neighbors to discard set

fig = plt.figure(figsize=[12, 6])
#fig.subplots_adjust(hspace=0.4, wspace=0.4)
ax = plt.subplot(111)

ax.plot(
    wgri,
    'o',
    alpha=1.,
    color='black',
    label=r'full basis set',
    marker='.',
    markersize=14)

ax.plot(
    wgri_strat,
    'o',
    alpha=1.,
    color='orange',
    label=r'reduced basis set',
    marker='.',
    markersize=6)

plt.xlabel(r'$\gamma(\rho_1)\;\;\; [fm^{-2}]$')
plt.ylabel(r'$\gamma(\rho_2)\;\;\; [fm^{-2}]$')
plt.legend(loc='best')

fig.savefig('WidthXY.pdf')