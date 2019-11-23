# -*- coding: utf-8 -*-

import pyemma
#pyemma.__version__
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os.path
import sys
import subprocess
from subprocess import Popen, PIPE, STDOUT
#import commands
import argparse
import mdtraj as md


gmx = "/commun/gromacs/512/bin/gmx"
#gmx = "/usr/local/gromacs/bin/gmx_mpi"

def parse_xvg(filexvg):
    """
    Parse xvg file and extract the number of hydrogen bonds
    and return the nulber of hydrogen for each frame
    """
    hbond = []
    with open(filexvg, "r") as filin:
        for line in filin:
            if line[0] == "@" or line[0] == "#" or line[0] == "&":
                continue
            hbond.append(line.split()[1])
    hbond = np.array(hbond)
    return hbond.astype(np.float)



#################################
#			Functions			#
#################################


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

def save_figure(name, fig_dir):
	do_save = True
	if do_save:
		plt.savefig(fig_dir + name, bbox_inches='tight')


def centroid(traj, frameClust, nb, output_folder):
	"""
	find the centroid for each cluster and save them
	"""
	centroid = 0
	MeanRMSD = 10
	#Extract conformation the most close to the cluster
	topology = traj.topology
	#indice of each carbon alpha
	CA = traj.topology.select('name CA')
	#Transpose the numpy matrix to extract frame indice inside each cluster.
	frameClust = frameClust.T[1]
	#détermine which frame is the centroid
	print("determine the centroid for the cluster "+str(nb))
	for IndStruct in frameClust:
		rmsds = md.rmsd(target = traj[frameClust], reference = traj[IndStruct], atom_indices=CA)
		if np.mean(rmsds) < MeanRMSD:
			centroid = IndStruct
			MeanRMSD = np.mean(rmsds)
	traj[centroid].save_pdb(output_folder+"clust"+str(nb)+".pdb")



def compute_effectif_cluster(ind_clust, time0, dt, PathOut, replica,struct):
	"""
	Generate a PDB file for the first structure inside each clusters. Save also the
	caracteristic of cluster inside the file:  effectif_clust
	Arguments:
		_ind_clust: array which contains the list of each frame inside the cluster
		_time0: The time corresponding of the first frame
		_dt: Time step
		_PathOut: Path where md_ex.xtc is localizated
		_replica: Path where is the fist trajectory
	"""
	traj = md.load(replica, top=struct)
	effectif_clust = []
	tot_effectif = 0
	for i in range(len(ind_clust)):
		effectif_clust.append(len(ind_clust[i]))
		tot_effectif += len(ind_clust[i])
	effectif_clust = np.array(effectif_clust)
	filout = open(PathOut+'/effectif_clust', 'w')
	filout.write("% effectif\tcluster\tframe (ps)\trange_RMSD (A)\t range_Rgyra(A)\n")
	for i in range(len(ind_clust)):
		pourcent = (100*effectif_clust/float(tot_effectif))[i]
		filout.write(str(round(pourcent,2))+"\t"+str(i+1)+"\t"+str(time0+ind_clust[i][0,1]*dt)+"\n")
		centroid(traj, ind_clust[i], i+1, PathOut)


def gyrateXVG(path_file):
    """
    Read xvg
    """
    time = []
    gyrate = []
    with open(path_file, 'r') as fichier1:
        for line in fichier1:
            if line[0] != "#" and line[0] != "@":
                columns = line.split()
                time.append(float(columns[0]))
                gyrate.append(float(columns[1]))
    return time, gyrate

def save_figure(name, fig_dir):
	do_save = True
	if do_save:
		plt.savefig(fig_dir + name, bbox_inches='tight')

def exportAvergaeStructure(aver_strucut, traj, PathOut):
    with open(PathOut+"average.pdb", "w") as filout:
        n_atoms = traj.n_atoms
        for i in range(n_atoms):
            atom = topology.atom(i)
            #convert nm to angstrom
            x = aver_coord[i][0]*10
            y = aver_coord[i][1]*10
            z = aver_coord[i][2]*10
            #atom.index
            #atom.name
            #atom.n_bonds
            #atom.residue.name
            resid = str(atom.residue)[3:]
            resn = str(atom.residue)[:3]
            filout.write("ATOM  {:5d} {:^4s} {:3s} A{:4d} \
    {:8.3f}{:8.3f}{:8.3f}{:6.2f}      {:>2s}  \n".format(\
                         atom.index+1, atom.name, resn, int(resid),\
                         x, y,z, 1.00, atom.name[0] ) )

#################################
#			Main				#
#################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description="This script generate free energy graph based on the radus of gyration and hydrogen bonds.")
    parser.add_argument('-f', action="store",dest="f", type=str,    help='name  of the trajectory analyzed.\nExemple: -f traj.xtc')
    parser.add_argument('-g', action="store", dest="g", type=str,
    help='structure file: .gro, .pdb\nExemple: -g struct.gro')
    parser.add_argument('-s', action="store", dest="s", type=str,
    help='tpr file: .tpr')
    parser.add_argument('-o', action="store", dest="o", type=str,
    help='folder where the output are saved', default="./free_energy_map/")
    parser.add_argument('-kmeans', action="store", dest="kmeans", type=int,
    help='Number of clusters', default=10)
    parser.add_argument('-stru4RMSD', action="store", dest="stru4RMSD", type=str,
    help='Number of clusters', default=None)
    args = parser.parse_args()


    print(args)
    replica = args.f
    struct = args.g
    tpr = args.s
    PathOut = args.o

    cleanFolder(PathOut)

    traj = md.load(replica, top=struct)
    topology = traj.topology
    n_frames = traj.n_frames
        #center coordinates to the origin
    traj.center_coordinates()
    print("Simulation has {0} frames and length is {1} ps".format(n_frames, traj.time[-1]))
    print("Remove the first 10%")
    traj[int(n_frames/10.0):].save_xtc(PathOut+"ex_md.xtc")
    replica = PathOut+"ex_md.xtc"



    replica=PathOut+"ex_md.xtc"
    traj = md.load(replica, top=struct)
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms

    #backbone's atoms indice
    bb = topology.select('name N or name CA or name C')
    #Fit all structure's backbone to the first frame
    traj.superpose(reference = traj[0], frame = 0, atom_indices = bb)
    topology = traj.topology

        #Compute radius of gyration
    gyrateArray = md.compute_rg(traj)
        #Compute average structure if there is no reference structure for the RMSD
    if args.stru4RMSD is None:
        cmd = gmx+" covar -f "+replica+" -s "+tpr+" -av "+PathOut+"average.gro"
        #Structure moyenne obtenu sur la traj ex_md
        Popen("echo \"4 0\" | "+cmd, shell=True).wait()
        ref_struct = md.load(PathOut+"average.gro")
        rms = md.rmsd(traj, ref_struct, atom_indices=bb)*10
    else:
        ref_struct = md.load(args.stru4RMSD)
        topology2 = ref_struct.topology
        bb2 = topology2.select('name N or name CA or name C')
        print("Superpose atoms:\n{0}\nwith\n{0}\n".format(bb,bb2))
        rms = md.rmsd(traj, ref_struct, atom_indices=bb)*10
        print("There is {0:.2f} %of conformation with a RMSD below 2.0 A".format(len(np.where(rms<2.0)[0])/len(rms)*100.0))
    #Convert in A
    gyrateArray = gyrateArray*10
    peptide_name = "peptide"
    #Passe d'un tableau python à un tableau numpy, obligatoire

    xmin=round(np.min(gyrateArray),2)-0.01
    xmax=round(np.max(gyrateArray),2)+0.01
    #plt.xlim([xmin,xmax])	#Borne
    #on fixe les borne pour avoir les mêmes échelles entre les graphes
    ymin=round(np.min(rms),2)
    ymax=round(np.max(rms),2)+0.05
    plt.xlim([xmin,xmax])


    #Trace la carte d'énergie libre, abscisse : rayon de gyration, ordonnee : RMSD
    plt.xlim([xmin,xmax])	#Borne
    plt.ylim([ymin,ymax])
    mplt.plot_free_energy(gyrateArray, rms)
    #plt.plot([Refgyrate],[ymin], '+')
    plt.ylabel('RMSD (A)')
    plt.xlabel('Radius of gyration (A)')
    save_figure('free'+peptide_name+'.pdf',PathOut+"/")	#Par défaut, image au format pdf

    #Mise en place du k-means
    n_clusters = args.kmeans
    Y = np.vstack((gyrateArray, rms))
    X = np.transpose(Y)
    clustering = coor.cluster_kmeans(X,k=n_clusters, max_iter=100)
    dtrajs = clustering.dtrajs

    cc_x = clustering.clustercenters[:,0]
    cc_y = clustering.clustercenters[:,1]
    ind_clust = clustering.index_clusters



    plt.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=5, color='black')
    for i in range(len(cc_x)):
        plt.text(cc_x[i], cc_y[i], str(i+1), color='grey', fontsize=12)

    save_figure('free'+peptide_name+'_clusters.pdf',PathOut)
    compute_effectif_cluster(ind_clust, int(traj.time[0]), int(traj.timestep), PathOut, PathOut+"ex_md.xtc", struct)
