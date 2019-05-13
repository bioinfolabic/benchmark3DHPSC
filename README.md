# Biological Sequence for the Protein Structure Problem with Integer Programming using 3D-HP-SC model

Requeriments:

1) Python 2.7

2) biopython (https://biopython.org/)

3) Java

4) IBM CPLEX

5) Linux 14.04 LTS

6) Rasmol 2.7.5.2

# Protein Sequence Dataset

The codes used for extracting the biological sequence from PDB and for the conversion to the HP sequence are contained in the EXT_CONV folder.

The program developed to extract the information and convert to the HP sequence is in the pdb3.py file, the pseudo-code is shown in the figure below. The hydrophobicity classification list is contained in the file alberts.txt, however, it can be modified if the user needs.

![alt text](https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/EXT_CONV/alg.png )


# How to execute the code:

$ java -Djava.library.path=<CPLEX_PATH> -jar Protein.jar <INPUT_FILE> 1 3 1 <INITIAL_STRUCUTRE> > <OUTPUT_FILE>

# Dictionary:

- CPLEX_PATH: Path of the cplex library, in this work we use version 10.4.1.

- Protein.jar: Executable program.

- INPUT_FILE: File with hydrophobic-polar sequence information.

- INITIAL_STRUCUTRE: It is possible to insert an initial structure for the optimization, otherwise, it can be replaced by the string "no".

- OUTPUT_FILE (optional): In this example, we save the information obtained from the I/O to a file.

# Example

java -Djava.library.path=/home/users/ILOG/CPLEX_Studio125/cplex/bin/x86-64_sles10_4.1 -jar Protein.jar "INPUT/protein_data_2GB1.txt" 1 3 1 "no" > output_2GB1.txt

# Related Works
```
@article{nunes2016integer,
  title={An integer programming model for protein structure prediction using the 3{D}--{HP} side chain model},
  author={Nunes, Luiz Fernando and Galv{\~a}o, Lauro Cesar and Lopes, Heitor Silv{\'e}rio and Moscato, Pablo and Berretta, Regina},
  journal={Discrete Applied Mathematics},
  volume={198},
  number={1},
  pages={206--214},
  year={2016}
}
```
