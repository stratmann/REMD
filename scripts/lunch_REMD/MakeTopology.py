#!/usr/bin/env python
# coding: utf8

import argparse
from subprocess import Popen, PIPE
import sys
import re
import time
import os

RgFlagCat = re.compile('^\[ [a-z]+ \]')
RgFlagDefault = re.compile('\[ defaults \]')
RgFlagAtomType = re.compile('\[ atomtypes \]')


def ReadTopFile(Filename) :
    """
    Read a gromacs top file and add or delete informations.
    """
    Flag = 0
    newFile = ''
    with open(Filename, 'r') as filin :
        for li in filin :
            if RgFlagCat.search(li) :
                if RgFlagDefault.search(li) or RgFlagAtomType.search(li) :
                    Flag = 1
                else :
                    Flag = 0
            if Flag == 0 :
                newFile += li
            elif Flag == 1 :
                newFile += ';'
                newFile += li
            else :
                print("Error in Flag")
    return newFile


def WriteNewFile(Filename, newFile, forceField=96) :
    """
    Write a new file readable for gromacs.
    """
    with open(Filename.replace(".top", "_corrected.top"),'w') as filout :
        filout.write(';force field\n#include "amber')
        filout.write(str(forceField))
        filout.write('.ff/forcefield.itp"\n\n')
        filout.write(newFile)


def Regex4topol(line, topfile):
    """
    Modify order atoms to correct improper dihedral
    """
    AtomID = {}
    atomName = line.split(";")[1]
    atomNumber = line.split(";")[0]
    for i in range(len(atomName.split())):
        AtomID[atomName.split()[i]] = atomNumber.split()[i]
    if line.split(";")[1].split()[0][0] == "N":
        #Modifiy initial order of N-    CA-     C-     O to CA-     N-     C-     O
        newline = '{0:>6} {1:>6} {2:>6} {3:>6}'.format(AtomID["CA-"], AtomID["N-"], AtomID["C-"],AtomID["O"])
        newline = newline+line[27:60]+"    CA-     N-     C-     O "
    elif line.split(";")[1].split()[0][0] == "H":
        #Modifiy initial order of H-     N-     C-    CA to C-    CA     N-     H
        newline = '{0:>6} {1:>6} {2:>6} {3:>6}'.format(AtomID["C-"], AtomID["CA"], AtomID["N-"],AtomID["H-"])
        newline =newline+line[27:60]+"     C-    CA      N-     H "
    else:
        print("Regex not found !\nPlease send a email with your PDB file")
        #print line
        sys.exit(0)
    #newline = newline+line[27:]
    print(newline)
    cmd = "sed -i \"s/"+line[:-1]+"/"+newline[:-1]+"/\" "+topfile
    Popen(cmd, shell=True).wait()


def CorrectDihedral(topfile1):
    """
    Amber define a wrong improper dihedral (the fisrt one)
    The function modify topology file to correct it
    Argument:
        _topfile: GROMACS topology generated by acpype (from amber topology)
    """
    flag = False
    cpt = 0
    with open(topfile1, "r") as filin1:
        for line1 in filin1:
            if line1[:-1] == "[ dihedrals ] ; impropers":
                flag = True
                continue
            elif flag is True and line1[0] == ";":
                continue
            elif flag is True and line1[0] != ";":
                tmp = line1
                cpt += 1
                Regex4topol(line1, topfile1)
                print(line1)
                if cpt > 1:
                    flag = False




def CorrectFirstRes(groFile):
    with open(groFile, "r") as file:
        for line in file:
            if len(line.split()) != 7:
                continue
            else:
                if line.split()[1] == "PRO":
                    print("first amino acid is a proline, no need to correct topology file")
                    return False
                else:
                    print("Improper dihedral correction")
                    return True

def ParseCA(pdb):
    cpt = 0
    with open(pdb, "r") as filin:
        for line in filin:
            if len(line) < 16:
                continue
            if line[12:16].split()[0] == "CA":
                cpt += 1
    return str(cpt)


def ChangeChirality(pdb, residues, tleap, outputs="./", subfold="./"):
    """
    Create a tleap script to modify the chirality ofselected residues
    /!\ Proline chirality is not modify
    Arguments:
        _pdb: pdb file name you want to modify
        _residues: residues you want to change (integer list)
        _outputs: the current path
        _subfold: path where pdb is
    Return:
        _new pdb files
    """
    os.chdir(subfold)
    print("generate script for amber")
    with open("chirality.leap", "w") as filin:
        filin.write("source "+tleap)
        filin.write("\nset default PBradii bondi")
        filin.write("\nclearpdbresmap")
        filin.write("\nmol = loadpdb "+pdb)
        for i in residues:
            filin.write("\nselect mol."+str(i)+".CA")
        filin.write("\nflip mol")
        filin.write("\nsavepdb mol "+pdb)
        filin.write("\nquit")
        print("generate amber's topology")
    Popen("tleap -f chirality.leap", shell=True).wait()
    os.chdir(outputs)


def makeTopology(structure, peptide, gmx, tleap, acpype, forcefield, output = "./"):
    """
    Generate topopoly files for gromacs from intial structure
    Arguments:
        _structure: PDB file
        _peptide: Flag to specify if the peptide is cyclic or not
        _gmx: path for gromacs
        _tleap: path for leap
    """
    os.chdir(output)
    if peptide == False:
        print("No cyclic peptide")
        cmd = gmx+" pdb2gmx -p peptide.top -ignh yes -ff "+forcefield+" -water none\
     -o peptide.gro -f "+structure
        Popen(cmd, shell=True).wait()
    else:
        print("Cyclic peptide")
        #First step we remove hhydrogen
        Popen("grep -v 'H' "+structure+" > pept-H.pdb", shell=True).wait()
        #Atoms's number beging at 1
        Popen(gmx+" editconf -f pept-H.pdb -o pept-good.pdb -resnr 1", shell=True).wait()
        residue = ParseCA("pept-good.pdb")
        #Create leap script to generate amber topology files
        #gromacs can not make correct topology file for cyclic peptide)
        print("generate script for amber")
        with open("script.leap", "w") as filin:
            filin.write("source "+tleap)
            filin.write("\nset default PBradii bondi")
            filin.write("\nclearpdbresmap")
            filin.write("\nmol = loadpdb pept-good.pdb")
            print("\nbond mol.1.N mol."+residue+".C")
            filin.write("\nbond mol.1.N mol."+residue.strip(" ")+".C")
            filin.write("\nsaveamberparm mol pept_amber.prmtop pept_amber.inpcrd")
            filin.write("\nsavepdb mol pept_amber.pdb")
            filin.write("\nquit")
        print("generate amber's topology")
        Popen("tleap -f script.leap", shell=True).wait()
        #delete unacessary files
        Popen("rm pept-H.pdb pept-good.pdb pept_amber.pdb", shell=True).wait()
        Popen(acpype+" -x pept_amber.inpcrd -p pept_amber.prmtop", shell=True).wait()
        #If the first residu is a proline, no need to correct topology file
        NewFile = ReadTopFile("pept_amber_GMX.top")
        WriteNewFile("pept_amber_GMX.top", NewFile)
        Popen("rm pept_amber.inpcrd pept_amber.prmtop pept_amber_GMX.top", shell=True).wait()
        if CorrectFirstRes("pept_amber_GMX.gro"):
            CorrectDihedral("pept_amber_GMX_corrected.top")
        Popen("mv pept_amber_GMX.gro peptide.gro", shell=True).wait()
        Popen("mv pept_amber_GMX_corrected.top peptide.top", shell=True).wait()




if __name__ == '__main__':
    #path for gromacs
    #gmx="/usr/local/gromacs/bin/gmx_mpi"
    #path for the force field
    #leap="/amber18/dat/leap/cmd/oldff/leaprc.ff96"
    #acpype="/home/MM_PBSA/src/acpype/acpype.py"
    #path for gromacs
    gmx="/commun/gromacs/512/bin/gmx"
    #path for the force field
    leap="/home/jaysen/amber14/dat/leap/cmd/oldff/leaprc.ff96"
    acpype="/bigdata/jaysen/Rapport_MAUD/Scripts/acpype/acpype.py"
    parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
    parser.add_argument('-g', action="store", dest="g", type=str, help="pdb file", default= None)
    parser.add_argument('-cyclic', action="store", dest="cyclic",default=False, type=bool, help="flag for cyclic peptide")
    parser.add_argument('-o', action="store", dest="o", type=str, default = "./",\
    help="output path filename ")
    parser.add_argument('-ff', action="store", dest="ff",default="amber96", type=str, help="force field use")
    arg = parser.parse_args()
    pdb = arg.g

    if pdb[-3:] != "pdb":
        print("please provide a PDB file")
        sys.exit(0)
    makeTopology(pdb, arg.cyclic, gmx, leap, acpype, arg.ff)
