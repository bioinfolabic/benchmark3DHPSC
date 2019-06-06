# Biological Sequence for the Protein Structure Problem with Integer Programming using 3D-HP-SC model

Requeriments:

1) Python 2.7

2) biopython (https://biopython.org/)

3) Java

4) IBM CPLEX (https://ibm.onthehub.com)

5) Linux 14.04 LTS or higher

6) Rasmol 2.7.5.2 (optional for drawing the proteins structure predicted)

# Protein Sequence Dataset

The codes used for extracting the biological sequence from PDB and for the conversion to the HP sequence are contained in the EXT_CONV folder.

The program developed to extract the information and convert to the HP sequence is in the pdb3.py file, the pseudo-code is shown in the figure below. The hydrophobicity classification list is contained in the file alberts.txt, however, it can be modified if the user needs.


<div id="fig:subfigures" class="subfigures" data-caption="Caption for figure" class="center">
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/EXT_CONV/alg.png" width="350" height="350">  
 
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/EXT_CONV/alberts.png" width="350" height="350">  

</div>


The table below shows the 19 sequences used for the Protein Structure Prediction problem.

<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/sequences.png" width="400" height="400" class="center">  



# How to execute the code:

$ java -Djava.library.path=<CPLEX_PATH> -jar Protein.jar <INPUT_FILE> 1 3 1 <INITIAL_STRUCUTRE> > <OUTPUT_FILE>
*Certifique-se de ter instalado as dependências necessárias. 
# Dictionary:

- CPLEX_PATH: Path of the cplex library, in this work we use version 10.4.1.

- Protein.jar: Executable program.

- INPUT_FILE: File with hydrophobic-polar sequence information.

- INITIAL_STRUCUTRE: It is possible to insert an initial structure for the optimization, otherwise, it can be replaced by the string "no".

- OUTPUT_FILE (optional): In this example, we save the information obtained from the I/O to a file.

# Example

java -Djava.library.path=/home/users/ILOG/CPLEX_Studio125/cplex/bin/x86-64_sles10_4.1 -jar Protein.jar "INPUT/protein_data_2GB1.txt" 1 3 1 "no" > output_2GB1.txt

*make sure that you have been installed the dependencies required.

# Results

The following are some proteins structure predicted (1DPK, 2IWJ, 1JBL, 1M23, 2IWJ, 2LKE, 2NVJ and 2BOY)  by the integer programming method with the 3D-HP-SC model. Subsequently, in the table, the results of the benchmarks are presented, showing the PDB ID of the protein, sequence size, hydrophobic-hydrophobic contact quantity, hydrophobic number of elements in the sequence, the computational time required for the prediction and predicted structure.

<div id="fig:subfigures" class="subfigures" data-caption="Caption for figure">
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/1DPK.png" width="200" height="200">  
 
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/2IWJ.png" width="200" height="200">  

<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/1JBL.png" width="200" height="200">  
 
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/1M23.png" width="200" height="200">  

\

<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/2IWJ.png" width="200" height="200">  
 
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/2LKE.png" width="200" height="200">  

<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/2NVJ.png" width="200" height="200">  
 
<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/B0Y.png" width="200" height="200">  
</div>


<img src="https://github.com/bioinfolabic/benchmark3DHPSC/blob/master/IMG/results.png" width="600" height="400" class="center">  

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
