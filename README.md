# Author: Jaysen SAWMYNADEN Guillaume POSTIC
#                               LICENCE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the [GNU General Public License](http://www.gnu.org/licenses/)
along with this program.

# REMD


## Automated Replica Exchange Molecular Dynamics protocol

##                               Principle

The script is a automated Replica Exchange Molecular Dynamics protocol.



##                               Install and execute
To use it, first install all the requiered module:
```
argparse
subprocess
sys
re
```
By default, the user can provide a sequence file, the script will create a reference structure (PDB file)
and a topology file. 
To lunch the REMD users must provide also temperature reference (for each replic).
Then a minimization, equilibration and md production with gromacs with the force field Amber 96 are made
```
 main.py -seq seq.txt -temperature "300 318 337.97 358.81 380.85 404.27 429.12 455.50"
 
```
 Arguments bellow are optional
```
 -g: reference PDB used to create topology file (default None)
 -p: topology file (default None)
 -cyclic: flag for cyclic peptide (default False)
 -temperature: temperature reference
 -time: time simulation in ps (default 200000) 
 -o: output path filenamec (default current folder)
 -log: Log file's name (file.log)
```

##                              Example
You can execute the script to create a PDB and topology files:
```
python main.py -seq seq.txt
```
If you want to create a cyclic peptide you have to use the flag -cyclic
```
python main.py -seq seq.txt -cyclic True
```
To lunch REMD you must provide temperature for replica:
```
python main.py -seq seq.txt -temperature "300 318 337.97 358.81 380.85 404.27 429.12 455.50"
```
It is possible to use and other force field by providing a PDB and topology files:
```
python main.py -seq seq.txt -temperature "300 318 337.97 358.81 380.85 404.27 429.12 455.50" -g ref.pdb -p topology.top
```
##                               Output

2 files:
* Log file (arguments used and other info)
* PDB file contain the structures
Length unit is (and MUST be) ALWAYS in angstrom whatever your MD engine

#                               Warning

The script use tleap (amberTools 18), acpype, GROMACS 5.1.2 and scwrl 3.0. Be sure the softwares are installed

##                              UPCOMING FEATURES

* PEP 8
* Ramachadran plot for each residues
* Clustering
