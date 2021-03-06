#!/usr/bin/env python
# coding: utf8

import argparse
from subprocess import Popen, PIPE
import sys
import re
import time
import os
import Bio.PDB
import numpy as np


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



def numb(resid):
    resid = str(resid)
    while(len(resid) < 4):
        resid = " "+resid
    return resid


def ChangePDB(filin, new_aa, Pro2Swap):
    """
    Replace amino acid coordinate by the one in new_aa
    """
    newPDB = ""
    with open(filin, "r") as pdb:
        for line in pdb:
            if line.split()[0] == "CONECT":
                continue
            if line.split()[0] != "ATOM":
                newPDB += line
            elif int(line[22:26]) == Pro2Swap and line[12:16].replace(" ", "") == "N":
                #open the new_aa
                with open(new_aa, "r") as newproline:
                    for line2 in newproline:
                        if line2.split()[0] == "ATOM":
                            newPDB += line2[:22]+numb(Pro2Swap)+line2[26:]
            elif int(line[22:26]) == Pro2Swap:
                continue
            else:
                newPDB += line
    #Export the newPDB
    with open(filin, "w") as filout:
    #with open("test.pdb", "w") as filout:
        filout.write(newPDB)
    #print(newPDB)





def PutThatDpro(pdb, patron, Pro2Swap):
    """
    Swapp amino acid to D-proline
    Arguments:
        _pdb: structure you want to modify
        _patron: D-proline model use
        _Pro2Swap:  residue you want to change
    Return:
        _pdb file: D-proline's new coordinate (that
    Source code:
    #http://combichem.blogspot.com/2013/08/aligning-pdb-structures-with-biopython.html
    """
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    peptide = pdb_parser.get_structure("reference", pdb)
    Dpro = pdb_parser.get_structure("Dpro", patron)
    ref_model = peptide[0]
    patron_model = Dpro[0]
    # Make a list of the atoms (in the structures) you wish to align.
    # In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    sample_atoms = []
    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
            print(ref_res.get_id()[1])
            # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in np.array([Pro2Swap]):
                # Append backbone atoms to list
                ref_atoms.append(ref_res['N'])
                ref_atoms.append(ref_res['CA'])
                ref_atoms.append(ref_res['C'])
    # Do the same for the sample structure
    for sample_chain in patron_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in np.array([1]):
                sample_atoms.append(sample_res['N'])
                sample_atoms.append(sample_res['CA'])
                sample_atoms.append(sample_res['C'])
    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(patron_model.get_atoms())
    # Print RMSD:
    print(super_imposer.rms)
    io = Bio.PDB.PDBIO()
    io.set_structure(patron_model)
    io.save("newD-pro.pdb")
    ChangePDB(pdb, "newD-pro.pdb", Pro2Swap)


def ParseCA(pdb):
    cpt = 0
    with open(pdb, "r") as filin:
        for line in filin:
            if len(line) < 16:
                continue
            if line[12:16].split()[0] == "CA":
                cpt += 1
    return str(cpt)


def ChangeChirality(pdb, residues, tleap, acpype, outputs="./", subfold="./", cyclic=False, seq=0):
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
        if cyclic:
            filin.write("\nbond mol.1.N mol."+str(len(seq))+".C")
        for i in residues:
            filin.write("\nselect mol."+str(i)+".CA")
        filin.write("\nflip mol")
        filin.write("\nsaveamberparm mol pept_amber.prmtop pept_amber.inpcrd")
        filin.write("\nsavepdb mol "+pdb)
        filin.write("\nquit")
        print("generate amber's topology")
    Popen("tleap -f chirality.leap", shell=True).wait()
    Popen(acpype+" -x pept_amber.inpcrd -p pept_amber.prmtop", shell=True).wait()
    NewFile = ReadTopFile("pept_amber_GMX.top")
    WriteNewFile("pept_amber_GMX.top", NewFile)
    Popen("rm pept_amber.inpcrd pept_amber.prmtop pept_amber_GMX.top ", shell=True).wait()
    Popen("rm pept-H.pdb pept-good.pdb pept_amber.pdb", shell=True).wait()
    if CorrectFirstRes("pept_amber_GMX.gro"):
            CorrectDihedral("pept_amber_GMX_corrected.top")
    Popen("mv pept_amber_GMX.gro peptide.gro", shell=True).wait()
    Popen("mv pept_amber_GMX_corrected.top peptide.top", shell=True).wait()
    os.chdir(outputs)
    return "peptide.top"

def removeH(structure):
    newPDB =""
    with open(structure, "r") as pdb:
        for line in pdb:
            if line.split()[0] != "ATOM":
                newPDB += line
            elif line[12:16].find("H") >= 0:
                #print(line[12:16])
                continue
            else:
                newPDB += line
    with open("pept-H.pdb", "w") as filout:
    #with open("test.pdb", "w") as filout:
        filout.write(newPDB)

def MadeTleapScript4Cyclic(tleap, output, residue, cyx):
    with open(output, "w") as filin:
        filin.write("source "+tleap)
        filin.write("\nset default PBradii bondi")
        filin.write("\nclearpdbresmap")
        filin.write("\nmol = loadpdb pept-good.pdb")
        print("\nbond mol.1.N mol."+residue+".C")
        filin.write("\nbond mol.1.N mol."+residue.strip(" ")+".C")
        if len(cyx) > 0:
            filin.write("\nbond mol."+str(cyx[0])+".SG mol."+str(cyx[1])+".SG")
        filin.write("\nsaveamberparm mol pept_amber.prmtop pept_amber.inpcrd")
        filin.write("\nsavepdb mol pept_amber.pdb")
        filin.write("\nquit")
    print("generate amber's topology")


def MadeTleapScript4Linear(tleap, output, cyx, structure):
    with open(output, "w") as filin:
        filin.write("source "+tleap)
        filin.write("\nset default PBradii bondi")
        filin.write("\nclearpdbresmap")
        filin.write("\nmol = loadpdb "+structure)
        print("\nbond mol."+str(cyx[0])+".SG mol."+str(cyx[1])+".SG")
        filin.write("\nbond mol."+str(cyx[0])+".SG mol."+str(cyx[1])+".SG")
        filin.write("\nsaveamberparm mol pept_amber.prmtop pept_amber.inpcrd")
        filin.write("\nsavepdb mol pept_amber.pdb")
        filin.write("\nquit")
    print("generate amber's topology")

def cyxInPDB(structure):
    """
    Return residue's number of all the CYX
    """
    cyx = []
    with open(structure, "r") as filin:
        for line in filin:
            if line.split()[0] != "ATOM":
                continue
            elif line.split()[2] == "CA" and line.split()[3] == "CYX":
                cyx.append(int(line[22:26]))
    return cyx


def tleapCaped(tleap, structure):
    with open("linear.leap", "w") as filin:
        filin.write("source "+tleap)
        filin.write("\nset default PBradii bondi")
        filin.write("\nclearpdbresmap")
        filin.write("\nmol = loadpdb "+structure)
        filin.write("\nsavepdb mol peptide.pdb")
        filin.write("\nquit")
    print("generate amber's topology")


def tmp_caped(newPDB, CH3, OC1, OC2):
    with open("tmp_caped.pdb", "w") as filout:
        lines = newPDB.split("\n")
        for i in range(len(lines)):
            if i == CH3[0]-2:
                #print(CH3[1])
                filout.write("ATOM      0  CH3 ACE{0}".format(CH3[1][21:]))
                #print(lines[i])
            elif i == OC1[0]-1:
                filout.write("{0}O  {1}".format(OC1[1][:13], OC1[1][16:]))
                #print(OC1[1])
                continue
            elif i == OC2[0]-1:
                #print("NME")
                filout.write("ATOM      0  N   NME{0}".format(OC2[1][20:]))
                continue
            else:
                filout.write(lines[i]+"\n")



def removeH4linear(tleap, structure):
    newPDB =""
    CH3 = []
    OC1 = []
    OC2 = []
    flag = True
    cpt_line = 0
    with open(structure, "r") as pdb:
        for line in pdb:
            if line[12:16].find("H1") >= 0 and flag is True:
                CH3 = [cpt_line, line]
                cpt_line += 1
                flag = False
                print(line)
            if line.split()[0] != "ATOM":
                newPDB += line
                cpt_line += 1
            elif line[12:16].find("H") >= 0:
                #print(line[12:16])
                continue
            else:
                if line[12:16].find("OC1") >= 0:
                    OC1 = [cpt_line, line]
                elif line[12:16].find("OC2") >= 0:
                    OC2 = [cpt_line ,line]
                newPDB += line
                cpt_line += 1
    tmp_caped(newPDB, CH3, OC1, OC2)
    tleapCaped(tleap, "tmp_caped.pdb")
    Popen("tleap -f linear.leap", shell=True).wait()
    #Remove TER to maintain peptide's integrity
    Popen("sed -i \"/TER*/d\" peptide.pdb", shell=True).wait()


def makeTopology(structure, peptide, gmx, tleap, acpype, forcefield, output = "./", debug = False):
    """
    Generate topopoly files for gromacs from intial structure
    Arguments:
        _structure: PDB file
        _peptide: Flag to specify if the peptide is cyclic or not
        _gmx: path for gromacs
        _tleap: path for leap
        _debug: Keep all output
    """
    os.chdir(output)
    cyx = cyxInPDB(structure)
    if peptide is False:
        print("No cyclic peptide")
        #check if there is disulfure bonds
        if len(cyx) > 0:
            print("ADD disulfure bond !!!")
            MadeTleapScript4Linear(tleap, "script.leap", cyx, structure)
            Popen("tleap -f script.leap", shell=True).wait()
            Popen(acpype+" -x pept_amber.inpcrd -p pept_amber.prmtop", shell=True).wait()
            NewFile = ReadTopFile("pept_amber_GMX.top")
            WriteNewFile("pept_amber_GMX.top", NewFile)
            Popen("mv pept_amber_GMX.gro peptide.gro", shell=True).wait()
            Popen("mv pept_amber_GMX_corrected.top peptide.top", shell=True).wait()
            Popen("rm pept_amber.inpcrd pept_amber.prmtop pept_amber_GMX.top pept-H.pdb pept-good.pdb pept_amber.pdb", shell=True).wait()
        #if not. Generate topology files with gromacs
        else:
            cmd = gmx+" pdb2gmx -p peptide.top -ignh yes -ff "+forcefield+" -water none \
-o "+structure+" -f "+structure
            #Popen("rm *.gro peptide.top *#", shell=True).wait()
            Popen(cmd, shell=True).wait()
            removeH4linear(tleap, structure)
            cmd = gmx+" pdb2gmx -p peptide.top -ignh yes -ff "+forcefield+" -water none \
-o peptide.gro -f peptide.pdb"
            Popen(cmd, shell=True).wait()
            #Popen("rm peptide.pdb *# tmp_caped.pdb", shell=True).wait()
    else:
        print("Cyclic peptide")
        #First step we remove hhydrogen
        #Popen("grep -v 'H' "+structure+" > pept-H.pdb", shell=True).wait()
        removeH(structure)
        #Atoms's number beging at 1
        Popen(gmx+" editconf -f pept-H.pdb -o pept-good.pdb -resnr 1", shell=True).wait()
        residue = ParseCA("pept-good.pdb")
        #Create leap script to generate amber topology files
        #gromacs can not make correct topology file for cyclic peptide)
        print("generate script for amber")
        MadeTleapScript4Cyclic(tleap, "script.leap", residue, cyx)
        Popen("tleap -f script.leap", shell=True).wait()
        #delete unacessary files
        #Popen("rm pept-H.pdb pept-good.pdb pept_amber.pdb", shell=True).wait()
        Popen(acpype+" -x pept_amber.inpcrd -p pept_amber.prmtop", shell=True).wait()
        #If the first residu is a proline, no need to correct topology file
        NewFile = ReadTopFile("pept_amber_GMX.top")
        WriteNewFile("pept_amber_GMX.top", NewFile)
        if debug is False:
            Popen("rm pept_amber.inpcrd pept_amber.prmtop pept_amber_GMX.top ", shell=True).wait()
            Popen("rm pept-H.pdb pept-good.pdb pept_amber.pdb", shell=True).wait()
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
