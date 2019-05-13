"""
=========================================================
        FEDERAL UNIVERSITY OF TECHNOLOGY - PARANA        
=========================================================

   NOTES:
---------------------------------------------------------
   
   Hydrophobic / Non-Polar / H = 0
   Hydrophilic / Polar / P = 1

---------------------------------------------------------
"""

from Bio.PDB import *
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np


parser = MMCIFParser()

protein = raw_input('Protein PDB code >> ')
structure = parser.get_structure(protein, protein+'.cif')
mmcif_dict = MMCIF2Dict(protein+'.cif')
sequence = mmcif_dict['_entity_poly_seq.mon_id']
len_seq = len(sequence)

print 'Number of amino acids ', len_seq

with open('seq_'+protein+'.txt', "wb") as fo:
	for i in range(len_seq):
		text = str(sequence[i])
		if i != len_seq-1:
			text = text + '\n'
		fo.write(text)
fo.close()


###  Reading the hydrophobicity scale  ###
scale_file = raw_input('Hydrophobicity scale file >> ')
with open(scale, 'r') as fi:
	data = fi.read()
fi.close()

lines = data.strip().split('\n')
scale = dict()
for line in lines:
	amino, hydro = line.split(':')
	scale[amino] = hydro.strip()


###  Writing the hp sequence (as 1s and 0s)  ###
with open('seq_hp_'+protein+'.txt', 'wb') as fo:
	fo.write('Number of residues of the protein\n'+str(len_seq)+'\n')
#	fo.write('Dimension of one of the sides of the cube\n5\n') # PODE TIRAR?
	fo.write('hp sequence\n') # COMO ESTAMOS CONVERTENDO PARA 0 OU 1, MUDAMOS ESTE NOME?
	for i in range(len_seq):
		try:
			n = float(scale[sequence[i]])
			print n
			if n >= 0:
				hp = '0\n'
			else:
				hp = '1\n'
			fo.write(hp)
		except ValueError:
			if scale[sequence[i]] == 'A':
				hp = '0\n'
			else:
				hp = '1\n'
			fo.write(hp)
	fo.write('break line') # ISTO ESTAVA NO CODIGO ANTIGO, PODE TIRAR?
fo.close()
