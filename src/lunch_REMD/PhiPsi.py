#!/usr/bin/env python
# coding: utf-8
import numpy as np
import math
import mdtraj

#########################################################
#       Compute Phi Psi angle foe each residue          #
#   Work as well for cyclic peptide and regular peptide #
#########################################################

def normalVect(p1,p2,p3):
    """
    Return vector normal to the plane
    """
    # These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1
    # the cross product is a vector normal to the plane
    return np.cross(v1, v2)

def ComputeAngle(vect1, vect2):
    """
    compute the angle between two vector
    """
    a1 = vect1[0]
    b1 = vect1[1]
    c1 = vect1[2]
    a2 = vect2[0]
    b2 = vect2[1]
    c2 = vect2[2]
    d = ( a1 * a2 + b1 * b2 + c1 * c2 )
    e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    return math.degrees(math.acos(round(d,4)))


def ComputeDihedral(traj, IJKL):
    """
    Compute dihedral angle by using 4 points
    """
    angle = np.array([])
    for frame in xrange(len(traj)):
        points = traj.xyz[frame][IJKL]
        #points = np.array([[1,0,0], [0,0,0], [0,1,0], [0,0,1]])
        #print points
        p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        p4 = points[3]
        vect1 = normalVect(p1,p2,p3)
        vect2 = normalVect(p2,p3,p4)
		#print vect1
		#print vect2
        angle = np.append(angle, ComputeAngle(vect1, vect2))
    return angle

def info(traj, IJKL):
    topology = traj.topology
    print "Atoms used to computed the dihedral angle:"
    for atom in IJKL:
        print "{0}\t{1}".format(atom,topology.atom(atom))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Use 4 points (ijkl) to compute\
the dihedral angle between two plans (IJK, JKL).\n/!\ Be careful number begin\
to 0. The scrit can read trajectory files to product histogram and even made a\
comparaison frequency if a second trajectory is provided.\n\
Example: python pathCV_S_Z.py -g input.pdb')
    parser.add_argument('-g', action="store", dest="g", type=str,
help="Topolgy file (PDB).\nExample: reference.pdb")
    parser.add_argument('-IJKL', action="store", dest="id", type=str,
help="Atom serial to compute the dihedral angle.\nExample: -IJKL \"1 2 3 4\"")
    parser.add_argument('-f1', action="store", dest="f1", type=str,
help="trajectory file")

    arg = parser.parse_args()
	#print "Parameters used (∩｀-´)⊃━☆ﾟ\n{0}.".format(arg)
    pdb = arg.g
    IJKL = np.fromstring(arg.id, dtype=int, sep=' ')
	#print "Wrong: N-    CA-     C-     O (1     60     72     73)\n\
#Correct: CA-     N-     C-     O (60      1     72     73)"
    if arg.f1 is None:
        traj = mdtraj.load(pdb)
        topology = traj.topology
    else:
        traj = mdtraj.load(arg.f1, top=pdb)
        topology = traj.topology


print "Phi angle"
phi = []
psi = []
AtomN = topology.select('name N and backbone')
AtomC = topology.select('name C and backbone')
AtomCA = topology.select('name CA and backbone')
#Phi atoms are Ci-1-N-Ca-C
#Psi N-Ca-Ci+1-Ni+1
for i in  xrange(topology.n_residues-1):
    IJKL = [AtomC[i],AtomN[i+1], AtomCA[i+1], AtomC[i+1]]
    ComputeDihedral(traj, IJKL)


mdtraj.compute_phi(traj)
val1 = ComputeDihedral(traj, IJKL)

mdtraj.compute_phi(traj)


    info(traj, IJKL)
    print "Compute dihedral angle:"
    val1 = ComputeDihedral(traj, IJKL)
    print "{0}\t{1}".format(np.mean(val1),np.std(val1))

    if val2 is not None:
        info(traj2, IJKL)
        print "{0}\t{1}".format(np.mean(val2),np.std(val2))
