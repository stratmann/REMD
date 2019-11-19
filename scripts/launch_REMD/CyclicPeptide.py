#!/usr/bin/env python3
# coding: utf8

import argparse
import os
import glob
from subprocess import Popen, PIPE
import sys
import shutil




code = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
"F", "P", "S", "T", "W", "Y", "V", "X"] #X corresponding to cystein with
#                                       with disulfure bridge


def parseSeq(fichier, path):
    """
    Read text and take amino acid sequences
    Check if the amino acid exit or not
    Arguments:
        fichier: file countaining amino acid sequences
        path: where output are stored
    """
    path = path+'seqLibrary/'
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    seqences = {} #dico countaining all the sequence
    flag_seq = True
    with open(fichier, 'r') as lines:
        for line in lines:
            if line[0] == "#":
                continue
            seqname = ""
            seq_pept = line.split('\n')[0] #only the sequence
            if flag_seq:
                if len(seq_pept) < 4:
                    print("Peptide too short\n")
                    sys.exit()
                flag_seq = False
            for aa in seq_pept:
                #Check if the amino acid are correct
                if aa.upper() in code:
                    seqname += aa
                else:
                    print("Sorry letter {0} is not a amino acid".format(aa))
                    sys.exit()
            fo = open(path + seqname + ".txt", 'w')
            #print("seqence {0}: {1}\n".format(i, seqname.upper()))
            fo.write(seqname.upper().replace("X", "C") + "\n") #For cycstein
            fo.close()

def cys2modify(peptide):
    """
    Read amino acid sequence and return position of X
    """
    print("Sequence : {0}".format(peptide))
    cyx = [] #residues to modify
    for i in range(len(peptide)):
        if peptide[i].upper() == "X":
            cyx.append(i+1) #sequence beging to 1
    return cyx

def cys2cyx(structure, cyx):
    lines = ""
    with open(structure, "r") as filin:
        for line in filin:
            if line.split()[0] != "ATOM":
                lines += line
                continue
            elif int(line.split()[5]) in cyx:
                lines += line.replace("CYS","CYX")
            else:
                lines += line

    with open(structure, "w") as filin:
        filin.write(lines)

def MakePeptideGreatAgain(scwrl, seq, cyclicPept, output="./"):
    """
    Arguments:
        _seq: files
    Return:
        _file name
    """
    lenPeptide = str(len(seq.split("/")[-1][:-4]))
    if cyclicPept is True:
        refbackbone ="/home/REMD/src/scwrl3_lin/CyclicPeptide/"+lenPeptide+".pdb"
    else:
        refbackbone ="/home/REMD/src/scwrl3_lin/LinearPeptide/"+lenPeptide+".pdb"
    peptide = seq.split("/")[-1][:-4]
    #tmp_pept = seq.replace("x","c").replace("X","C")
    Popen(scwrl+" -i "+refbackbone+" -s "+seq+" -o "+output+"ref"+peptide+".pdb", shell=True).wait()
    print(scwrl+" -i "+refbackbone+" -s "+seq+" -o "+output+"ref"+peptide+".pdb")
    if os.path.exists(output+"ref"+peptide+".pdb") is True:
        print("peptide structure was made.")
        structure = output+"ref"+peptide+".pdb"
        #Modify pdb to add disulfure bonds (CYS to CYS)
        if peptide.upper().find("X") >= 0:
            cyx = cys2modify(peptide) #list of residues to modify
            cys2cyx(structure, cyx)
        return structure
    else:
        print("peptide structure was not made. :'('")
        sys.exit(0)




if __name__ == '__main__':
    scwrl="/bigdata/jaysen/these/github/Docker/REMD/src/scwrl3_lin/./scwrl3"
    Popen("rm -rf ./seqLibrary", shell=True).wait()
    parser = argparse.ArgumentParser(description='topology fie (GRO,PDB)')
    parser.add_argument('-seq', action="store", dest="s", type=str, help="sequence file")
    parser.add_argument('-cyclic', action="store", dest="cyclic", type=bool, default=False, help="flag for cyclic peptide")
    arg = parser.parse_args()

    #Read the file to create the sequence and the pdb
    parseSeq(arg.s)
    for elmt in glob.glob("./seqLibrary/*.txt"):
        newfolder=elmt.split("/")[-1][:-4]
        #if not os.path.exists('./peptide/'+newfolder):
        #    os.makedirs('./peptide/'+newfolder)
        MakePeptideGreatAgain(scwrl, elmt, arg.cyclic)
