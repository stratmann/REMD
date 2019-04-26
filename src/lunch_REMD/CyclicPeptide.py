#!/usr/bin/env python
# coding: utf8

import argparse
import os
import glob
from subprocess import Popen, PIPE

#########################################################
#       Function to create peptide's structure file     #
#   Work as well for cyclic peptide and regular peptide #
#########################################################


code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'}


def parseSeq(fichier):
    if not os.path.exists('./seqLibrary'):
        os.makedirs('./seqLibrary')
    seq = []
    seqList = []
    fi = open(fichier, 'r')
    for line in fi.readlines():
        items = line.split()
        seq.append(items[1:])
    fi.close()
    for res1 in seq[0]:
        for res2 in seq[1]:
            for res3 in seq[2]:
                for res4 in seq[3]:
                    for res5 in seq[4]:
                        seqname = code[res1] + code[res2] + code[res3] + code[res4] + code[res5]
                        fo = open("./seqLibrary/" + seqname + ".txt", 'w')
                        fo.write(seqname + "\n")
                        fo.close()


def MakePeptideGreatAgain(scwrl, seq, cyclicPept, output="./"):
    """
    Arguments:
        _seq: files
    """
    lenPeptide = str(len(seq.split("/")[-1][:-4]))

    if cyclicPept is True:
        refbackbone ="../src/scwrl3_lin/BackboneReference/CyclicPeptide/"+lenPeptide+".pdb"
    else:
        refbackbone ="../src/scwrl3_lin/BackboneReference/LinearPeptide/"+lenPeptide+".pdb"
    Popen(scwrl+" -i "+refbackbone+" -s "+seq+" -o ref.pdb", shell=True).wait()
    if os.path.exists('ref.pdb') is True:
        print "peptide structure was made."
    else:
        print "peptide structure was not made. :'('"
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
