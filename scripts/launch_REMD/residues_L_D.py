# coding: utf8

import numpy as np
import mdtraj as md
import argparse


def XYZ(coord):
	"""
	Convert coordinate to string and check the good number of caracters
	Argument:
		_coord: float
	"""
	while len(str(coord).split(".")[1]) < 3:
		coord = str(coord)+"0"
	while len(str(coord).split(".")[1]) > 3:
		coord = str(coord)[:-1]
	while len(str(coord).split(".")[0]) < 3:
		coord = " "+str(coord)
	return str(coord)



def atoms_residue(topology, i):
	resid = []
	resid.append(topology.select('name N and backbone and resid '+str(i))[0])
	resid.append(topology.select('name CA and backbone and resid '+str(i))[0])
	resid.append(topology.select('name HA and resid '+str(i))[0])
	resid.append(topology.select('name C and backbone and resid '+str(i))[0])
	resid.append(topology.select('name CB and resid '+str(i))[0])
	return resid


def coordinates(resid_atm, topol):
	"""
	Arguments:
		_resid_atm: dictionary countain backbone's atom number
		_topology: trajectory file
	Return: Matrix with atomistic position
	"""
	#print(topol.xyz[0][resid_atm])
	return topol.xyz[0][resid_atm]



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

def align_CO_2_X(mat):
    """
    Rotate alpha carbon on Z axis to fit CO to X axis
    Arguments:
        mat:cordinate Matrix
    Return
        tmpMat: cordinate matrix
    """
    Rangle = 0
    tmpMat = rotate_matrix(mat, np.array([0,0,0]) + np.array([0,0,1])*0, False)
    for angle in np.arange(0,360,5):
        tmp = rotate_matrix(mat, np.array([0,0,0]) + np.array([0,0,1])*angle, False)
        if tmp[3][0] > tmpMat[3][0]:
            Rangle = angle
            tmpMat = tmp
    #print("Angle {0} ° new coordinates {1}".format(np.array([0,0,1])*Rangle, tmpMat))
    return tmpMat, angle

def export_carbon(amino_acid, coord):
    with open(amino_acid+".pdb", "w") as filout:
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





def export_residue(topology, i, angles, coord, Rz):
    ind_CA = topology.select('name CA and backbone and resid '+str(i))[0]
    atoms = topology.select('resid '+str(i))
    atom = str(topology.atom(atoms[0]))
    aa = atom.split("-")[0]
    cood_CA = coord[ind_CA]
    center_coord = (coord - cood_CA)
    #rotate coordinates
    tmp = rotate_matrix(center_coord, angles, False)
    #rotate Z axis
    rotate_coord = rotate_matrix(tmp, np.array(
                        [0,0,0]) + np.array([0,0,1])*Rz,\
                        False)
    with open(aa+".pdb", "w") as filout:
        for elmt in atoms:
            atom = str(topology.atom(elmt))
            aa = atom.split("-")[0]
            atom_type = atom.split("-")[1]
            filout.write("ATOM  {:5d} {:^4s} {:3s} A{:4d}  \
{:8.3f}{:8.3f}{:8.3f}  1.00  0.00           \n".format(elmt+1, atom_type, aa[:3], int(aa[3]),\
        rotate_coord[elmt][0], rotate_coord[elmt][1], rotate_coord[elmt][2]))

def center_residue(topology, i, coord):
    ind_CA = topology.select('name CA and backbone and resid '+str(i))[0]
    atoms = topology.select('resid '+str(i))
    atom = str(topology.atom(atoms[0]))
    aa = atom.split("-")[0]
    cood_CA = coord[ind_CA]
    return coord - cood_CA



parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
parser.add_argument('-g', action="store", dest="g", type=str, default = None, help="pdb file (default None)")
arg = parser.parse_args()
print(arg)
pdb = arg.g

ref_struct = md.load(pdb)
topology = ref_struct.topology
#Take backbone's atoms index
resid_atm = atoms_residue(topology)
for frame in range(ref_struct.n_frames):
    for ind in resid_atm.keys():
        #for frame in in range(ref_struct.n_frames):
        #convert to angstrom unit useless. Avoid only little number
        center2origin = center_residue(topology, ind, ref_struct.xyz[frame]*10)
        #residue's atoms (N, CA, HA, C, CB)
        ind_bb = atoms_residue(topology, ind)
        #method not optimal but it is working
        #Align HA to Z axis
        angles = Align2Z(center2origin[ind_bb[2]])
        tmp = rotate_matrix(center2origin, angles, False)
        #then compute angle between CO-Ca-N
        amino_acid_rotated, Rz = align_CO_2_X(tmp)
        print(amino_acid_rotated)
        #Compare position of N, if it is on the left then amino acid is L form
        #in other hand if it on right it is D form
        if str(topology.residue(ind))[:3] == "GLY":
            print("{0} is L form".format(topology.residue(ind)))
            #export_carbon(str(topology.residue(ind)), center2origin)
            coninue
        elif str(topology.residue(ind))[:3] != "PRO":
            if amino_acid_rotated[0][1] > 0:
                print("{0} is L form".format(topology.residue(ind)))
            else:
                print("{0} is D form".format(topology.residue(ind)))
        elif str(topology.residue(ind))[:3] == "PRO":
            if amino_acid_rotated[4][1] > 0:
                print("{0} is D form".format(topology.residue(ind)))
            else:
                print("{0} is L form".format(topology.residue(ind)))
        export_residue(topology, ind, angles, ref_struct.xyz[frame]*10, Rz)
