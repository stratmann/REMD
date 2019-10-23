#!/usr/bin/env python3
# coding: utf8

import argparse
import os
import glob
from subprocess import Popen, PIPE
import sys
import shutil


code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'}


def parseSeq(fichier, path):
    path = path+'seqLibrary/'
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    seqences = {} #dico countaining all the sequence
    flag_seq = True
    with open(fichier, 'r') as lines:
        for line in lines:
            items = line.split()
            #check if the first arguments is the amino acid's number
            if flag_seq:
                if items[0] != "1":
                    print("Format not good :(\nThe first id is not 1")
                    sys.exit()
                flag_seq = False
            if len(seqences) == 0:
                for i in range(1,len(items)):
                    seqences[i] = [items[i]]
            else:
                for i in range(1,len(items)):
                    seqences[i].append(items[i])
    #Convert code 3 letters to 1 letters
    for i in seqences.keys():
        seqname = ""
        for amino_acid in seqences[i]:
            if amino_acid.isupper() and amino_acid in code:
                seqname += code[amino_acid]
            elif amino_acid.islower() and amino_acid.upper() in code:
                seqname += code[amino_acid.upper()].lower()
                #seqname += code[amino_acid.upper()]
                #D- amino acid not yet implemented
            else:
                print("Amino acid not know\n")
        fo = open(path + seqname + ".txt", 'w')
        print("seqence {0}: {1}\n".format(i, seqname.upper()))
        fo.write(seqname.upper() + "\n")
        fo.close()


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
    Popen(scwrl+" -i "+refbackbone+" -s "+seq+" -o "+output+"ref"+peptide+".pdb", shell=True).wait()
    print(scwrl+" -i "+refbackbone+" -s "+seq+" -o "+output+"ref"+peptide+".pdb")
    if os.path.exists(output+"ref"+peptide+".pdb") is True:
        print("peptide structure was made.")
        return output+"ref"+peptide+".pdb"
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
