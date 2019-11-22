#!/usr/bin/env python3
# coding: utf8

import argparse
import subprocess
from subprocess import Popen, PIPE
import sys
import re
from CyclicPeptide import *
from MakeTopology import *
import random

gmx="/usr/local/gromacs/bin/gmx_mpi"
#path for the force field
leap="/amber18/dat/leap/cmd/oldff/leaprc."
acpype="/home/REMD/src/acpype/acpype.py"
scwrl="/home/REMD/src/scwrl3_lin/./scwrl3"

#gmx="/commun/gromacs/512/bin/gmx"
#path for the force field
#leap="/home/jaysen/amber14/dat/leap/cmd/oldff/leaprc.ff96"
#acpype="/bigdata/jaysen/Rapport_MAUD/Scripts/acpype/acpype.py"
#scwrl="/bigdata/jaysen/these/github/Docker/REMD/src/scwrl3_lin/./scwrl3"

leapFF={"amber96":"ff96", "amber99sb":"ff99SB", "amber99sb-ildn":"ff99SBildn", "amber03":"ff03"}


def minimization(path="./"):
    os.chdir(path)
    print("Minimization")
    Popen("rm -rf mini1 mini2 REMD *#", shell=True).wait()
    Popen("mkdir mini1 mini2 REMD", shell=True).wait()
    Popen("cp *.itp *.top *.gro ./mini1/", shell=True).wait()
    Popen("cp /home/REMD/src/minim_good.mdp "+path+"mini1", shell=True).wait()
    Popen("cp /home/REMD/src/Equil.mdp "+path+"mini2", shell=True).wait()
    Popen("cp /home/REMD/src/md_good.mdp "+path+"REMD", shell=True).wait()
    Popen(gmx+" editconf -f "+path+"mini1/peptide.gro -o "+path+"mini1/center_mini_systeme.gro  -bt cubic -d 1.0 -c yes", shell=True).wait()
    Popen(gmx+" grompp -f "+path+"mini1/minim_good.mdp -p "+path+"mini1/*.top -c "+path+"mini1/center_mini_systeme.gro -o "+path+"mini1/mini1.tpr", shell=True).wait()
    Popen(gmx+" mdrun -deffnm "+path+"mini1/mini1", shell=True).wait()
    Popen("cp "+path+"mini1/mini1.gro "+path+"mini2/mini1.gro", shell=True).wait()
    Popen("cp "+path+"mini1/*.top "+path+"mini2/", shell=True).wait()
    Popen("cp "+path+"mini1/*.itp "+path+"mini2/", shell=True).wait()



def equilibration(path="./"):
    #os.chdir(path)
    print("Equilibration")
    for i in range(len(refTemp)):
        Popen("sed \"s|ref_t = 300 ; .*|ref_t = "+str(refTemp[i])+" ;|\" "+path+"mini2/Equil.mdp>"+path+"mini2/Equil_"+str(i)+".mdp", shell=True).wait()
        Popen(gmx+" grompp -f "+path+"mini2/Equil_"+str(i)+".mdp -c "+path+"mini2/mini1.gro -p "+path+"mini2/*.top -o "+path+"mini2/Equil_"+str(i)+".tpr", shell=True).wait()
        Popen(gmx+" mdrun -v -deffnm "+path+"mini2/Equil_"+str(i), shell=True).wait()
    Popen("cp "+path+"mini2/Equil_* "+path+"REMD", shell=True).wait()
    Popen("cp "+path+"mini2/*.top "+path+"REMD", shell=True).wait()


def REMD(refTemp, path="./"):
    print("Time simulation given by the user {0}".format(arg.time))
    Popen("sed -i \"s|nsteps = 100000000 .*|nsteps = "+str(arg.time*500)+" ;|\" "+path+"REMD/md_good.mdp", shell=True).wait()
    for i in range(len(refTemp)):
        Popen("sed \"s|ref_t = 300 .*|ref_t = "+str(refTemp[i])+" ;|\" "+path+"REMD/md_good.mdp>"+path+"REMD/md_good_"+str(i)+".mdp", shell=True).wait()
        Popen(gmx+" grompp -f "+path+"REMD/md_good_"+str(i)+".mdp -c "+path+"REMD/Equil_"+str(i)+".gro -p "+path+"REMD/*.top -o "+path+"REMD/md_good1_"+str(i)+".tpr -t "+path+"REMD/Equil_"+str(i)+".cpt", shell=True).wait()
        Popen("rm "+path+"REMD/*#", shell=True).wait()
    print("utilisation de {0} replica".format(len(refTemp)))
    #command
    cmd = "mpirun -np "+str(len(refTemp))+" --allow-run-as-root "+gmx+" mdrun -s "+path+"REMD/md_good1_ -deffnm "+path+"REMD/md_good1_ -replex 500 -multi "+str(len(refTemp))
    Popen(cmd, shell=True).wait()

def lunch_REMD(path, refTEMP):
    minimization(path)
    equilibration(path)
    REMD(refTemp,path)
    Popen("rm "+path+"*#", shell=True).wait()

def lunch_analyze(subfold):
    os.makedirs(subfold+"analyze/")
    Popen("cp "+subfold+"REMD/md_good1_0* "+subfold+"analyze/", shell=True).wait()
    print("python /home/REMD/scripts/analyze_REMD/free_energy_map.py -f "+subfold+"analyse/md_good1_0.xtc -s "+subfold+"analyze/md_good1_0.tpr -g "+subfold+"analyze/md_good1_0.gro -o "+subfold+"analyze/")
    os.system("python /home/REMD/scripts/analyse_REMD/free_energy_map.py -f "+subfold+"analyze/md_good1_0.xtc -s "+subfold+"analyze/md_good1_0.tpr -g "+subfold+"analyze/md_good1_0.gro -o "+subfold+"analyze/free_energy_map/")
    Popen("cp "+subfold+"analyze/free_energy_map/ex_md.xtc "+subfold+"analyze/", shell=True).wait()
    os.system("python /home/REMD/scripts/analyse_REMD/clust_reg_space.py -f "+subfold+"analyze/ex_md.xtc -g "+subfold+"analyze/md_good1_0.gro -o "+subfold+"analyze/reg_space/")



parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
parser.add_argument('-g', action="store", dest="g", type=str, default = None, help="pdb file (default None)")
parser.add_argument('-p', action="store", dest="p", type=str, default = None, help="topology file (default None)")
parser.add_argument('-cyclic', action="store", dest="cyclic",default=False, type=bool, help="flag for cyclic peptide (default False)")
parser.add_argument('-temperature', action="store", dest="temp",default=None, type=str, help="\"300 313 329 347 367 391 418 450\"")
parser.add_argument('-time', action="store", dest="time",default=1500000, type=int, help="Time simulation in ps (default 1,500,000 ps)")
parser.add_argument('-log', action="store", dest="log", type=str,\
 default = "clust.log", help="log file's name: default clust.log")
parser.add_argument('-o', action="store", dest="o", type=str, \
default = None, help="output path filename ")
parser.add_argument('-seq', action="store", dest="s", type=str, help="sequence file", default= None)
parser.add_argument('-ff', action="store", dest="ff", type=str, \
default = "amber96", help="Force field use for the MD")


arg = parser.parse_args()
print(arg)
pdb = arg.g
sequence = arg.s
topology = arg.p
if arg.o is None:
    outputs = os.getcwd()+"/"+str(random.randint(10000,99999))+"/"
#sequence = "seq-example.txt"
else:
    outputs = arg.o
if arg.ff not in leapFF.keys():
    #if the force field is not implemented, use amber96
    leap += leapFF["amber96"]
else:
    leap += leapFF[arg.ff]


if arg.s is not None:
    print("Structure creation...")
    parseSeq(arg.s, outputs)
    for elmt in glob.glob(outputs+"seqLibrary/*.txt"):
        newfolder=elmt.split("/")[-1][:-4]
        MakePeptideGreatAgain(scwrl, elmt, arg.cyclic, outputs)
        pdb = "ref"+newfolder+".pdb"
        #new foldders for each peptides
        subfold = outputs+newfolder+"/"
        os.makedirs(subfold)
        Popen("mv "+outputs+pdb+" "+subfold, shell=True).wait()
        #Check if there is some D- amino acid
        if newfolder.isupper() is not True:
            ind = [] #position of the D- amino acid
            #convert L- amino acid to D- amino acid
            #Except proline (they are not convert to D form)
            for i in range(len(newfolder)):
                if newfolder[i].islower():
                    ind.append(i+1) #first residue is 1 in amber
            cmd = gmx+" editconf -f "+subfold+pdb+" -resnr 1 -o "+subfold+pdb
            Popen(cmd, shell=True).wait()
            ChangeChirality(pdb, ind, leap, outputs, subfold)
            #Replace amino acid to D-proline
            for i in range(len(newfolder)):
                if newfolder[i].islower() and newfolder[i].islower() == "p":
                    print("D-Pro in position {0}".format(i+1))
                    PutThatDpro(subfold+pdb, "/home/REMD/src/scwrl3_lin/patron_D-pro.pdb", int(i+1))
        if pdb[-3:] != "pdb":
                print("No PDB found")
                sys.exit(0)
        if topology is None:
                makeTopology(pdb, arg.cyclic, gmx, leap, acpype, "amber96", subfold)
        if arg.temp is None:
            print("No temperature given...\n")
            print("Final files:\nTopology: {0}\nStructure: {1}\n".format(topology, pdb))
            continue
        elif len(arg.temp) <= 1:
            print("Need more than 1 temperature to do REMD\n")
            continue
        refTemp = arg.temp.split()
        lunch_REMD(subfold, refTemp)
        lunch_analyze(subfold)
else:
    if arg.g is not None:
        pdb = arg.g
    if topology is None:
        subfold = outputs
        os.makedirs(subfold)
        Popen("cp "+pdb+" "+subfold, shell=True).wait()
        makeTopology(pdb, arg.cyclic, gmx, leap, acpype, "amber96", subfold)
    if arg.temp is None:
        print("No temperature given...\n")
        print("Final files:\nTopology: {0}\nStructure: {1}\n".format(topology, pdb))
        sys.exit()
    elif len(arg.temp) <= 1:
        print("Need more than 1 temperature to do REMD\n")
        sys.exit()
    else:
        refTemp = arg.temp.split()
        lunch_REMD(subfold, refTemp)
        lunch_analyze(subfold)
