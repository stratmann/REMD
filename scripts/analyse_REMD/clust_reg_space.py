# coding: utf8
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


def Boxplot(clust, clustmean):
    """
    The function plots the boxplot corresponding of the clusters which
    represent more 10% of total counts.
    """
    boxplt=[]
    boxplt5=[]
    for i in range(len(clustmean)):
        if clustmean[i] >= 0.1*len(clust):
            boxplt.append(dt[:,i])
        if clustmean[i] < 0.1*len(clust) and clustmean[i] >= 0.05*len(clust):
            boxplt5.append(dt[:,i])
    plt.boxplot(boxplt)
    plt.show()
    plt.boxplot(boxplt5)
    plt.show()

def PieChart(clust, threshold):
    chart = deepcopy(clust)
    chart.sort()
    labels = []
    count = []
    tmp = 0
    effct = 0
    explode=[0]
    for i in range(int(np.max(chart))):
        tmp = len(np.where(chart ==i+1)[0])
        if tmp >= threshold*len(chart):
            effct += tmp
            count.append(tmp)
            labels.append("clust"+str(i+1))
    count.append(int(len(clust)-effct))
    labels.append("other")
    return count, labels


def PlotPieChart(clust, threshold):
    size = len(clust)
    the_grid = GridSpec(2, 3)
    #row, col
    fig = plt.figure(tight_layout=True)
    #Global pie chart
    plt.subplot(the_grid[0, 1], aspect=1)
    count, labels = PieChart(clust, threshold)
    plt.pie(count, labels=labels, autopct='%1.1f%%', shadow=True)
    plt.title('All frames')
    #Sub pie chart
    for i in range(3):
        plt.subplot(the_grid[1, i], aspect=1)
        count1, labels1 = PieChart(clust[(i*size)/3:((i+1)*size)/3], threshold)
        plt.pie(count1, labels=labels1, autopct='%.0f%%', shadow=True)
        plt.title('Frames '+str(i*size/3)+' to '+str((i+1)*size/3))
    plt.show()


def ExportResults(clustmean, clustsd):
    """
    Export boostrap result to a log file
    Column 1: Index
    Column 2: Mean
    Column 3: Standard deviation
    """
    with open("./boot.log", "w") as results:
        results.write("Cluster_index\tMean\tSD\n")
        for i in range(len(clustmean)):
            results.write("{0}\t{1:.2f}\t{2:.2f}\n".format(i+1,clustmean[i],clustsd[i]))

def logfile(logfile, line):
    with open(logfile, "w+") as filin:
        filin.write(line)

def distance(array1, array2):
    dist = 0
    #print(array1.shape)
    #for i in range(array1.shape):
    return np.sum((array1 - array2)**2)


def cleanFolder(PathOut):
	"""
	clean erase the older files and folder
	Arguments:
		_PathOut: folder
	"""
	if os.path.exists(PathOut):
		remove = "rm -rf "+PathOut
		os.system(remove)
	os.makedirs(PathOut)


#################################
#			Main				#
#################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
    parser.add_argument('-g', action="store", dest="g", type=str,
                    help="reference.pdb")
    parser.add_argument('-f', action="store", dest="f", type=str,
                    help="trajectory file (xtc)")
    parser.add_argument('-dmin', action="store", dest="dmin", type=float, default=2.0,
                    help="minimum distance between cluster centers.")
    parser.add_argument('-log', action="store", dest="log", type=str,\
     default = "clust.log", help="log file's name: default clust.log")
    parser.add_argument('-o', action="store", dest="o", type=str, default = "./",\
    help="output path filename ")
    arg = parser.parse_args()

    #topfile = 'peptide_example.gro'
    topfile = arg.g
    #traj = 'md_example.xtc'
    traj = arg.f
    cleanFolder(arg.o)



    feat = coor.featurizer(topfile)
    feat.add_backbone_torsions(selstr= None, deg=True, cossin=True) # in degrees

    #List of all the angles
    #print(feat.describe())
    #Number of dregree of freedom
    #print(feat.dimension())

    inp = coor.source(traj, feat)
    sincos = inp.get_output()[0]

    #############
    #Use a regular space clustering. Cluster centers are at least in distance of
    #dmin to each other according to the given metric.Then Voronoi discretization
    #with the computed centers is used to partition the data
    cl_space = coor.cluster_regspace(sincos, dmin=2, max_centers = 100000)
    clustCenters = cl_space.clustercenters #angle for each centroid
    #We now discretize the trajectory to either set of cluster centers
    #assign structure's cluster number
    Sspace = coor.assign_to_centers(sincos, clustCenters)
    #assign for each cluster their frames number
    indexClusters = cl_space.index_clusters
    clustCentersFrameNo = -1*np.ones(len(clustCenters),dtype='int32')

    #Find the centroid (euclideen distance)
    for ind_clust in range(len(indexClusters)):
        #Frames which compose the cluster
        frameNumber = indexClusters[ind_clust][:,1]
        cosinusSinus = sincos[frameNumber,:]
        print(str(ind_clust))
        for j in range(len(cosinusSinus)):
            #If the frame feature is close to the centroid, then save it
            if (distance(clustCenters[ind_clust], cosinusSinus[j]) < 1e-7):
                line = 'Cluster ',str(ind_clust), 'Centroid ', str(frameNumber[j])
                if (clustCentersFrameNo[ind_clust] == -1):
                    clustCentersFrameNo[ind_clust] = frameNumber[j]
                    print(line)
                    #logfile(logfile, line)
                else:
            #Inform if there is an another potential centroid                print('Found an another potential centroid ! You could be more strigean')
                    print('Error ',line)
                    #logfile(logfile, line)


    #np.savetxt(arg.o+'_angles4_clustCentersFrameNo.txt', clustCentersFrameNo, fmt='%d')
    t = md.load(traj, top=topfile)
    #check results:
    #print(clustCentersFrameNo)
    with open(arg.o+'angles4_clustCentersFrameNo.txt', "w") as filout:
        #print('cluster no., centroid, nb_elemnt, pourcentage')
        filout.write('cluster no., centroid, nb_elemnt, pourcentage\n')
        for i in range(len(clustCenters)):
            pourcentage = np.shape(indexClusters[i])[0] *100.0/len(Sspace[0])
            if pourcentage > 5:
                #print(i, clustCentersFrameNo[i], np.shape(indexClusters[i])[0], np.round(pourcentage,2))
                filout.write(str(i)+','+str(clustCentersFrameNo[i])+','+str(np.shape(indexClusters[i])[0])+','+str( np.round(pourcentage,2) )+'\n' )
                t[clustCentersFrameNo[i]].save_pdb(arg.o+"centroid_clust"+str(i)+".pdb")
