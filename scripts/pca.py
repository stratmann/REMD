from __future__ import print_function
import pyemma
print(pyemma.__version__)

import os
import sys
import argparse
import numpy as np

import pyemma.coordinates as coor

import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from itertools import combinations
from copy import deepcopy

#################################
#               PCA             #
#################################
t = md.load(traj, top=topfile)
print(t)
topology = t.topology
print(topology)

pca = PCA(n_components=8)
#Select backbone (C1 is for the nmethy residue)
serial_backbone = topology.select("backbone")
#In this example we have a customized amino acid
#We have to selected them
serial_nmethyl = topology.select("(name N or name CA or name C or name O)\
 and (resid 0 or resid 7 or resid 8 or resid 15)")
backbone = np.concatenate((serial_backbone, serial_nmethyl), axis=0)
backbone.sort()
for i in backbone:
    print (topology.atom(i))



atom_pairs = list(combinations(backbone, 2))



pca = PCA(n_components=8)
atom_pairs = list(combinations(range(t.n_atoms), 2))
pairwise_distances = md.geometry.compute_distances(t, atom_pairs)
print(pairwise_distances.shape)
reduced_distances = pca.fit_transform(pairwise_distances)
pca.components_
pca.explained_variance_ratio_
tmp = deepcopy(pca.components_)



dist_space = coor.cluster_regspace(pairwise_distances, dmin=55)
centers_space = np.sort(dist_space.clustercenters, axis=0)
Sspace = coor.assign_to_centers(pairwise_distances, centers_space)
test = dist_space.dtrajs


plt.figure()
#plt.scatter(reduced_distances[:, 0], reduced_distances[:,1], marker='x', c=t.time)
#plt.scatter(reduced_distances[:, 0], reduced_distances[:,1], marker='x', c=[1]*len(reduced_distances[:,1]))
plt.scatter(reduced_distances[:, 0], reduced_distances[:,2], marker='o', c=test[0]*5, alpha=0.3)
plt.xlabel('PC1')
plt.ylabel('PC3')
plt.title('Pairwise distance PCA: cyclic peptide')
#cbar = plt.colorbar()
#cbar.set_label('Time [ps]')
plt.show()


####################################



feat = coor.featurizer(topfile)
feat.add_backbone_torsions(selstr= None, deg=True, cossin=True) # in degrees
#feat.add_chi1_torsions(deg=True)

#List of all the angles
#print(feat.describe())
#Number of dregree of freedom
#print(feat.dimension())

inp = coor.source(traj, feat)
sincos = inp.get_output()
#############
cl_space2 = coor.cluster_regspace(sincos, dmin=5)
centers_space2 = np.sort(cl_space2.clustercenters, axis=0)

#We now discretize the trajectory to either set of cluster centers
Sspace = coor.assign_to_centers(sincos, centers_space2)

pca2 = PCA(n_components=8)
reduced_distances = pca2.fit_transform(sincos[0])
pca2.components_
pca2.explained_variance_ratio_
tmp = deepcopy(pca2.components_)
