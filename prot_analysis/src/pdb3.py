"""
=========================================================
        FEDERAL UNIVERSITY OF TECHNOLOGY - PARANA        
=========================================================

   NOTES:
---------------------------------------------------------
   
   Hydrophobic / Non-Polar / A = 0
   Hydrophilic / Polar / B = 1

---------------------------------------------------------
"""

import re
import numpy as np
import os
import Bio
from Bio.PDB import *
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

parser = MMCIFParser()

protein = input('Protein PDB code >> ')

###  Download the mmCif of the molecule if not in the system  ###
sv=PDBList(server='https://files.rcsb.org/', pdb=None, obsolete_pdb=True, verbose=True)
sv.retrieve_pdb_file(protein, obsolete=False, pdir="proteins/"+protein, file_format="mmCif", overwrite=False)



###  Add the protein to the proteins.txt  ###
with open("data_input/proteins.txt",mode = 'a+') as fo:
	with open("data_input/proteins.txt",mode ='r+') as fr:
		for exists in fr: 	
			if(re.search(protein, exists)):
				flag=1
				break
			else: 
				flag=0
				
	if flag==0:
		fo.write(protein+"%s \n")
	fr.close()
fo.close()



structure = parser.get_structure(protein, "proteins/"+protein+"/"+protein+".cif")
mmcif_dict = MMCIF2Dict("proteins/"+protein+"/"+protein+".cif")
sequence = mmcif_dict['_entity_poly_seq.mon_id']
len_seq = len(sequence)
io=PDBIO()

print ('Number of amino acids ', len_seq)


with open("proteins/"+protein+"/seq_"+protein+".txt", "w+") as fo:
	for i in range(len_seq):
		text = str(sequence[i])
		if i != len_seq-1:
			text = text + '\n'
		fo.write(text)
fo.close()

###  Create a file .pdb withe the atoms  ###
for chain in structure.get_chains():
	io.set_structure(chain)
	io.save("proteins/"+protein+"/"+protein+'_'+chain.get_id() + ".pdb")
a_p="proteins/"+protein+"/"+protein+'_'+chain.get_id() + ".pdb"

###  Locate and separete the CAs  ###
def CA(pdbFileName):
	fr = open(pdbFileName, 'r')
	fw = open("proteins/"+protein+"/"+protein+"_CA.txt", 'w')

	for found in fr:
		if(re.search(r'^ATOM\s+\d+\s+CA+', found)):
			fw.write(found)

	fw.close()
	fr.close()

CA(a_p)

###  Reading the hydrophobicity scale  ###
scale_file = input('Hydrophobicity scale file >> ')
with open("data_input/"+scale_file+'.txt', 'r') as fi:
	data = fi.read()
fi.close()

lines = data.strip().split('\n')
scale = dict()
for line in lines:
	amino, hydro = line.split(':')
	scale[amino] = hydro.strip()


###  Writing the AB sequence (as 1s and 0s)  ###
with open("proteins/"+protein+"/seq_10_"+protein+".txt", 'w+') as myfile:  
	with open("proteins/"+protein+"/seq_ab_"+protein+".txt", 'w+') as fo:
		fo.write('Number of residues of the protein\n'+str(len_seq)+'\n')
	#	fo.write('Dimension of one of the sides of the cube\n5\n') # PODE TIRAR?
		fo.write('AB -10 sequence\n') 
		for i in range(len_seq):
			try:
				n = float(scale[sequence[i]])
				print (n)
				if n >= 0:
					ab = '0\n'
				else:
					ab = '1\n'
				fo.write(ab)
				myfile;write(ab)
			except ValueError:
				if scale[sequence[i]] == 'A':
					ab = '0\n'
				else:
					ab = '1\n'
				fo.write(ab)
				myfile.write(ab)
		fo.write('break line')
	fo.close()
myfile.close()



###  Write the cordenaates of the CAs Hidrophibics   ###

j=0
i=0
n=0
seq_hh=[]
num_matrix=[]
with open("proteins/"+protein+"/posicao_com_numeracao_"+protein+".txt", 'w+') as fa: 
	with open("proteins/"+protein+"/seq_10_"+protein+".txt", 'r') as fr: 
		with open("proteins/"+protein+"/"+protein+"_CA.txt", 'r') as myfile: 	
			with open("proteins/"+protein+"/"+protein+"_pos_CA_h.txt", 'w+') as fw: 
				for seq in fr:
					seq_hh.append(seq.strip())

				pattern = re.compile(r'(?:ATOM\s+\d+\s+CA+\s+\D+\s+\D+\s+\d+)+(\s+[\D\d]+\s+[\D\d]+\s+[\D\d])+(?:\d)');
				aux=[]
				for fo in myfile:
					for found in re.findall(pattern,fo):
						aux.append(found.strip().split()[0:3])					
						if (seq_hh[j]=='0'):
							fw.write(aux[j][0]+" "+aux[j][1]+" "+aux[j][2]+"\n")					
							fa.write("atomo:"+str(j)+" coordenadas:"+aux[j][0]+" "+aux[j][1]+" "+aux[j][2]+"\n")					
							num_matrix.append(j)
						j+=1
				with open("proteins/"+protein+"/"+protein+"_CA_euc_map.txt", 'w+') as fo: 
					for i in range(len(num_matrix)):
						for n in range(len(num_matrix)):
							fo.write("["+str(num_matrix[i])+"]["+str(num_matrix[n])+"] ")
						fo.write("\n")

			fw.close()
		myfile.close()
	fr.close()
fa.close()

print("Success!\n")