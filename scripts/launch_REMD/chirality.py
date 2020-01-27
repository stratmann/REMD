#!/usr/bin/env python3
# coding: utf8

import numpy as np
import mdtraj as md
import argparse
import shutil
import os
import copy



def atoms_residue(topology, i):
	resid = []
	resid.append(topology.select('name N and backbone and resid '+str(i))[0])
	resid.append(topology.select('name CA and backbone and resid '+str(i))[0])
	resid.append(topology.select('name HA and resid '+str(i))[0])
	resid.append(topology.select('name C and backbone and resid '+str(i))[0])
	resid.append(topology.select('name CB and resid '+str(i))[0])
	return resid



def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def MatRotX(thetaX):
	Rx = np.array((
					(1, 0, 0),
					(0,  np.cos(thetaX), -np.sin(thetaX)),
					(0,  np.sin(thetaX),  np.cos(thetaX))
					))
	return Rx

def MatRotY(thetaY):
	Ry = np.array((
					(np.cos(thetaY), 0, np.sin(thetaY)),
					(0,  1, 0),
					(-np.sin(thetaY),  0,  np.cos(thetaY))
					))
	return Ry

def MatRotZ(thetaZ):
	Rz = np.array((
					(np.cos(thetaZ), -np.sin(thetaZ), 0),
					(np.sin(thetaZ),  np.cos(thetaZ), 0),
					(0,  0,  1)
					))
	return Rz


def rotate_matrix(mat, angle, radian = True):
	"""
	The function center a matrix to zero, then applies a generic rotation
	Arguments:
		mat: coordinates
		angle: angle for the rotation (the right-hand rule)
		axis: which axis the rotation is made (1=x;2=y;3=z)
	Retrun:
		mat
	"""
	if radian is True:
		#print("use radian angle")
		thetaX = angle[0]
		thetaY = angle[1]
		thetaZ = angle[2]
	else:
		#print("use degree angle")
		thetaX = np.radians(angle[0])
		thetaY = np.radians(angle[1])
		thetaZ = np.radians(angle[2])
	newMat=[]
	#Rotate on X axis
	Rx = MatRotX(thetaX)
	#Rotate on Y axis
	Ry = MatRotY(thetaY)
	#Rotate on Z axis
	Rz = MatRotZ(thetaZ)
	#print('apply the rotation matrix to structure: r*v')
	#print( Rx.dot((np.transpose(mat)) ))
	#center to origin and apply translation
	R = np.dot(Rx, Ry)
	#print(R)
	R = np.dot(R, Rz)
	#print(R)
	return np.transpose( R.dot(np.transpose(mat)))



def minimize_angle(vector, vect_start, vect2, start=0, end=360, step=2):
    """
    Find angle which maximize coordinate Z
    Argument:
        _vector_HA: vector to rotate
        _vect_start: vector use to made rotation
        _vect2: vector use to increase vect_start
    """
    tmp = rotate_matrix(vector, vect_start + vect2*start, False)
    Rangle = start
    for angle in np.arange(start,end,step):
        tmp2 = rotate_matrix(vector, vect_start + vect2*angle, False)
        #If we increase Z dimension (it means we are closer Z axis)
        if tmp2[2] > tmp[2]:
            tmp = tmp2
            Rangle = angle
            #print("Angle {0} ° coordinate {1}".format(vect_start + vect2*angle, tmp))
        #print("Angle {0} ° coordinate {1}".format(angle, tmp2))
    #print("Angle {0} ° new coordinates {1}".format(vect_start + vect2*Rangle, tmp))
    return Rangle


def Align2Z(vector):
    X = minimize_angle(vector, np.array([0,0,0]), np.array([1,0,0]))
    #optimize X value
    X = minimize_angle(vector, np.array([0,0,0]), np.array([1,0,0]), X-2, X+2, 0.005)
    Y = minimize_angle(vector, np.array([X,0,0]), np.array([0,1,0]))
    Y = minimize_angle(vector, np.array([X,0,0]), np.array([0,1,0]), Y-2, Y+2, 0.005)
    Z = minimize_angle(vector, np.array([X,Y,0]), np.array([0,0,1]))
    Z = minimize_angle(vector, np.array([X,Y,0]), np.array([0,0,1]), Z-2, Z+2, 0.005)
    return np.array([X,Y,Z])

def compute_dist(v1, v2):
    tmp = np.square(v1[0]-v2[0]) + np.square(v1[1]-v2[1]) + np.square(v1[2]-v2[2])
    return np.sqrt(tmp)

def minimize_dist(coord, vector, axis, start = 0, end = 360, step = 10):
    Rangle = start
    rotated_coord = rotate_matrix(coord, vector, False)
    dist = compute_dist(coord[2], coord[4])
    for angle in np.arange(start, end, step):
        rotated_coord = rotate_matrix(coord, vector + axis*angle, False)
        #compute distance between HA ref and CB rotate
        tmp_dist = compute_dist(coord[2], rotated_coord[4])
        if tmp_dist < dist:
            Rangle = angle
            dist = tmp_dist
            #print("Dist: {0} angle {1}°".format(tmp_dist, Rangle))
    return Rangle

def swap_CB_2_HA(coord, ind_bb):
#order NA CA HA C CB
    ref_coord = copy.deepcopy(coord)[ind_bb]
    Rx = minimize_dist(ref_coord, np.array([0,0,0]), np.array([1,0,0]), 0,360,10)
    Rx = minimize_dist(ref_coord, np.array([0,0,0]), np.array([1,0,0]), Rx-10,Rx+10,0.5)
    Ry = minimize_dist(ref_coord, np.array([Rx,0,0]), np.array([0,1,0]), 0,360,10)
    Ry = minimize_dist(ref_coord, np.array([Rx,0,0]), np.array([0,1,0]), Ry-10,Ry+10,0.5)
    Rz = minimize_dist(ref_coord, np.array([Rx,Ry,0]), np.array([0,0,1]), 0,360,10)
    Rz = minimize_dist(ref_coord, np.array([Rx,Ry,Rz]), np.array([0,0,1]), Rz-10,Rz+10,0.5)
    return np.array([Rx,Ry,Rz])



def align_CO_2_X(mat):
    """
    Rotate alpha carbon on Z axis to fit CO to X axis
    Arguments:
        mat:cordinate Matrix
    Return
        newMat: cordinate matrix
    """
    Rangle = 0
    newMat = rotate_matrix(mat, np.array([0,0,0]) + np.array([0,0,1])*0, False)
    for angle in np.arange(0,360,5):
        tmpMat = rotate_matrix(mat, np.array([0,0,0]) + np.array([0,0,1])*angle, False)
        if tmpMat[3][0] > newMat[3][0]:
            Rangle = angle
            newMat = tmpMat
    #print("Angle {0} ° new coordinates {1}".format(np.array([0,0,1])*Rangle, tmpMat))
    return newMat, Rangle



def export_carbon(amino_acid, filenmae, coord):
    with open(filenmae, "w") as filout:
        filout.write("ATOM      1  N   {:3s} A   1  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n".format(amino_acid[:3], coord[0,0], coord[0,1], coord[0,2]))
        filout.write("ATOM      2  CA  {:3s} A   1  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n".format(amino_acid[:3], coord[1,0], coord[1,1], coord[1,2]))
        filout.write("ATOM      3  HA  {:3s} A   1  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n".format(amino_acid[:3], coord[2,0], coord[2,1], coord[2,2]))
        filout.write("ATOM      4  C   {:3s} A   1  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n".format(amino_acid[:3], coord[3,0], coord[3,1], coord[3,2]))
        if amino_acid[:3] != "GLY":
            filout.write("ATOM      5  CB  {:3s} A   1  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n".format(amino_acid[:3], coord[4,0], coord[4,1], coord[4,2]))





def export_residue(topology, i, coord, path = "./"):
    """
    Save residue's coordinates into a pdb file
    Arguments:
        _topology: mdtraj trajectory
        _i: residue number in according mdtraj (0 is the first residue)
        _coord: Matrix coordinate
    """
    ind_CA = topology.select('name CA and backbone and resid '+str(i))[0]
    atoms = topology.select('resid '+str(i))
    atom = str(topology.atom(atoms[0]))
    aa = atom.split("-")[0]
    with open(path+aa+".pdb", "w") as filout:
        for elmt in atoms:
            atom = str(topology.atom(elmt))
            aa = atom.split("-")[0]
            atom_type = atom.split("-")[1]
            filout.write("ATOM  {:5d} {:^4s} {:3s} A{:4d}  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00           \n".format(elmt+1, atom_type, aa[:3], int(aa[3]),\
        coord[elmt][0], coord[elmt][1], coord[elmt][2]))



def center_residue(topology, i, coord):
    ind_CA = topology.select('name CA and backbone and resid '+str(i))[0]
    atoms = topology.select('resid '+str(i))
    atom = str(topology.atom(atoms[0]))
    aa = atom.split("-")[0]
    cood_CA = coord[ind_CA]
    return coord - cood_CA



def change_chirality(ref_struct, aa):
    """
    Swap chirality for one amino acid
    Arguments:
        _ref_struct: topology (mdtraj)
        _aa: amino acid number
    Return coordinates matrix
    """
    topology = ref_struct.topology
    print("change chirality residue: {0}".format(topology.residue(aa)))
    #atoms's indice in the amino acid
    atoms_in_aa = topology.select('resid '+str(aa))
    #center amino acid's carbon alpha to the origin
    ind_CA = topology.select('name CA and resid '+str(aa))[0]
    #center2origin = center_residue(topology, aa, ref_struct.xyz[0]*10)
    center2origin = ref_struct.xyz[0] - ref_struct.xyz[0][ind_CA]
    #backbone's atoms indice
    ind_bb = atoms_residue(topology, aa) #order NA CA HA C CB
    #sidechain's atoms indice
    ind_sidechain = topology.select('sidechain and resid '+str(aa))
    ind_HA = topology.select('name HA and resid '+str(aa))
    #find angle to use to swap position
    angles = swap_CB_2_HA(center2origin, ind_bb)
    #Apply rotation matrix and translation to CB and HA
    new_sidechain = rotate_matrix(center2origin, angles, False) + ref_struct.xyz[0][ind_CA]
    new_HA = rotate_matrix(center2origin, -angles, False) + ref_struct.xyz[0][ind_CA]
    #save new coordinates
    new_coor = ref_struct.xyz[0]
    new_coor[ind_sidechain] = new_sidechain[ind_sidechain]
    new_coor[ind_HA] = new_HA[ind_HA]
    return new_coor


def export_pdb(ref_struct, coord, filout):
    """
    Arguments:
        _ref_struct: mdtraj.Trajectory
        _coord: cordinates matrix
        _filout: filename for the export
    """
    topology = ref_struct.topology
    with open(filout, "w") as output:
        for i in range(topology.n_atoms):
            tmp = str(topology.atom(i))
            atom_type = tmp.split("-")[-1]
            amino_acid = tmp.split("-")[0]
            output.write("ATOM  {:5d} {:^4s} {:3s} A{:4d}  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00           \n".format(i+1, atom_type, amino_acid[:3], int(amino_acid[3:]),\
        coord[i][0], coord[i][1], coord[i][2]))
        output.write("TER")


def automotize_swap(pdb, list_aa, output):
    cpt = 0
    for aa in list_aa:
        if cpt == 0:
            new_coord = change_chirality(md.load(pdb), int(aa))
            export_pdb(md.load(pdb), new_coord*10, output_pdb)
        if cpt > 0:
            new_coord = change_chirality(md.load(output_pdb), int(aa))
            export_pdb(md.load(output_pdb), new_coord*10, output_pdb)
        cpt += 1



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modify chirlity of amino acid')
    parser.add_argument('-g', action="store", dest="g", type=str, default = None, help="pdb file (default None)")
    parser.add_argument('-o', action="store", dest="o", type=str, \
    default = None, help="output filename ")
    parser.add_argument('-aa', action="store", dest="aa", type=str, \
    help="list of residues to change (first is 0). Example -aa \"0 2 4\"")
    arg = parser.parse_args()
    print(arg)
    pdb = arg.g

    pdb = arg.g
    #pdb = "mini1.pdb"

    if arg.o is None:
        output_pdb = "tmp-"+pdb
    else:
        output = arg.o

    list_aa = arg.aa.split()
    automotize_swap(pdb, list_aa, output_pdb)
