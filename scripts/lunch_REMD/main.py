#!/usr/bin/env python
# coding: utf8

import argparse
from subprocess import Popen, PIPE
import sys
import re
from CyclicPeptide import *
from MakeTopology import *

gmx="/usr/local/gromacs/bin/gmx_mpi"
#path for the force field
leap="/amber18/dat/leap/cmd/oldff/leaprc.ff96"
acpype="/home/REMD/src/acpype.py"

#gmx="/commun/gromacs/512/bin/gmx"
#path for the force field
#leap="/home/jaysen/amber14/dat/leap/cmd/oldff/leaprc.ff96"
#acpype="/bigdata/jaysen/Rapport_MAUD/Scripts/acpype/acpype.py"
#scwrl="/bigdata/jaysen/these/github/Docker/REMD/src/scwrl3_lin/./scwrl3"


parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
parser.add_argument('-g', action="store", dest="g", type=str, help="pdb file")
parser.add_argument('-p', action="store", dest="p", type=str, default = None, help="topology file")
parser.add_argument('-cyclic', action="store", dest="cyclic",default=False, type=bool, help="flag for cyclic peptide")
parser.add_argument('-temperature', action="store", dest="temp",default=None, type=str, help="\"300 313 329 347 367 391 418 450\"")
parser.add_argument('-time', action="store", dest="time",default=200000, type=int, help="Time simulation (ps)")
parser.add_argument('-log', action="store", dest="log", type=str,\
 default = "clust.log", help="log file's name: default clust.log")
parser.add_argument('-o', action="store", dest="o", type=str, default = "./",\
help="output path filename ")
parser.add_argument('-seq', action="store", dest="s", type=str, help="sequence file", default= None)

arg = parser.parse_args()
print arg
pdb = arg.g
sequence = arg.s
topology = arg.p

if arg.s is not None:
    print "Structure creation..."
    parseSeq(arg.s)
    for elmt in glob.glob("./seqLibrary/*.txt"):
        newfolder=elmt.split("/")[-1][:-4]
        MakePeptideGreatAgain(scwrl, elmt, arg.cyclic)
        pdb = "ref.pdb"
if pdb[-3:] != "pdb":
    print "please provide a PDB file"
    sys.exit(0)
if topology is None:
    pdb, topology = makeTopology(pdb, arg.cyclic, gmx, leap, acpype, "amber96")


if arg.temp is None:
    print "No temperature given...\n"
    print "Final files:\nTopology: {0}\nStructure: {1}\n".format(topology, pdb)
    sys.exit(0)
refTemp = arg.temp.split()
#####Realisation de la minimisation
print "Minimization"
Popen("rm -rf ./mini1 ./mini2 ./REMD *#", shell=True).wait()
Popen("mkdir mini1 mini2 REMD", shell=True).wait()
Popen("cp *.itp *.top *.gro "+pdb+" ./mini1", shell=True).wait()
Popen("cp /home/REMD/src/minim_good.mdp ./mini1", shell=True).wait()
Popen("cp /home/REMD/src/Equil.mdp ./mini2", shell=True).wait()
Popen("cp /home/REMD/src/md_good.mdp ./REMD", shell=True).wait()
Popen(gmx+" editconf -f ./mini1/"+pdb+" -o ./mini1/center_mini_systeme.gro  -bt cubic -d 1.0 -c yes", shell=True).wait()
Popen(gmx+" grompp -f ./mini1/minim_good.mdp -p ./mini1/"+topology+" -c ./mini1/center_mini_systeme.gro -o ./mini1/mini1.tpr", shell=True).wait()
Popen(gmx+" mdrun -deffnm ./mini1/mini1", shell=True).wait()

Popen("cp ./mini1/mini1.gro ./mini2/mini1.gro", shell=True).wait()
Popen("cp ./mini1/"+topology+" ./mini2/"+topology, shell=True).wait()
for i in xrange(len(refTemp)):
    Popen("sed \"s|ref_t = 300 ; .*|ref_t = "+str(refTemp[i])+" ;|\" ./mini2/Equil.mdp>./mini2/Equil_"+str(i)+".mdp", shell=True).wait()
    Popen(gmx+" grompp -f ./mini2/Equil_"+str(i)+".mdp -c ./mini2/mini1.gro -p ./mini2/"+topology+" -o ./mini2/Equil_"+str(i)+".tpr", shell=True).wait()
    Popen(gmx+" mdrun -v -deffnm ./mini2/Equil_"+str(i), shell=True).wait()


Popen("cp ./mini2/Equil_* ./REMD", shell=True).wait()
Popen("cp ./mini2/"+topology+" ./REMD", shell=True).wait()

print "Time simulation given by the user {0}".format(arg.time)
Popen("sed -i \"s|nsteps = 100000000 .*|nsteps = "+str(arg.time*500)+" ;|\" ./REMD/md_good.mdp", shell=True).wait()
for i in xrange(len(refTemp)):
    Popen("sed \"s|ref_t = 300 .*|ref_t = "+str(refTemp[i])+" ;|\" ./REMD/md_good.mdp>./REMD/md_good_"+str(i)+".mdp", shell=True).wait()
    Popen(gmx+" grompp -f md_good_"+str(i)+".mdp -c ./REMD/Equil_"+str(i)+".gro -p ./REMD/"+topology+" -o ./REMD/md_good1_"+str(i)+".tpr -t ./REMD/Equil_"+str(i)+".cpt", shell=True).wait()
    Popen("rm ./REMD/*#", shell=True).wait()

print "utilisation de {0} replica".format(len(refTemp))
