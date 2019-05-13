# Biological Sequence for the Protein Structure Problem with Integer Programming using 3D-HP-SC model

Requeriments:

1) Python

2) Java

3) IBM CPLEX

5) Linux 14.04 LTS

6) Rasmol 2.7.5.2

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
