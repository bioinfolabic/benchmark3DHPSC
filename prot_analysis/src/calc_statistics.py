"""
=========================================================
        FEDERAL UNIVERSITY OF TECHNOLOGY - PARANA        
=========================================================


---------------------------------------------------------
"""

import re 
import numpy as np
from scipy.spatial import distance
import os

print ("Initiating statistics...")
try:
    os.mkdir('data_output')
except OSError as e:
	del e


n=0
i=0
all_dist_euc=[]

###  Open protein list  ###
with open("data_input/proteins.txt",'r') as fr:

###  Open here will be saved all the euclidean distances  ###
	with open("data_output/dist_euc_all_prot.txt",'w+') as fi:

### Get coordenates  of hydrofobics amino acids's carbon alpha  ###
		for prot in fr:
			protein=prot.strip()

			with open("proteins/"+protein+"/"+protein+"_pos_CA_h.txt", 'r') as myfile: 
				aux=[]
				for fo in myfile:
					aux.append(fo.strip().split()[0:3])					
					i+=1
			myfile.close()				
			ca_pos=np.array(aux,dtype=np.float64)

###  Create file where will be saved each protein's euclidean distance among itself's hydrophibics residues  ###
			with open("proteins/"+protein+"/"+protein+"_euc_dist.txt",'w+') as fw:
				euc_dists=np.zeros((len(ca_pos),len(ca_pos)))

				for aux_euc in range(len(ca_pos)):
					for aux_euc2 in range(len(ca_pos)):
						euc_dists[aux_euc][aux_euc2]=distance.euclidean(ca_pos[aux_euc],ca_pos[aux_euc2])
						fw.write(str(euc_dists[aux_euc][aux_euc2])+" ")
					fw.write("\n")

				for aux_all_euc in range(len(ca_pos)):
					for aux_all_euc2 in range(aux_all_euc+1,len(ca_pos)):
						fi.write(str(euc_dists[aux_all_euc][aux_all_euc2])+" ")	
						all_dist_euc.append(euc_dists[aux_all_euc][aux_all_euc2])
						n+=1		


			fw.close()
	fi.close()
fr.close()

euc_dists=np.array(all_dist_euc,dtype=np.float64)

euc_dists.sort()

with open("data_output/all_euc_dist_asc.txt",'w+') as fw:
	fw.write(str(euc_dists))	
fw.close()


###  Calcule the euclidean distance mean  ####
mean_euc_dist=(np.sum(euc_dists))/len(euc_dists)


###  Get maximum , second to maximum and minime distances  ###
max_dist=np.max(euc_dists)

sec_max_dist=euc_dists[-2]

min_dist=np.min(euc_dists[np.nonzero(euc_dists)])


###  Calcule the Standard deviation  ###
desv_padrao=0
var=0
j=0

for aux_desv in range(len(euc_dists)):
		var= var + (euc_dists[aux_desv]-mean_euc_dist)**2	
		j+=1		

desv_padrao=(var/(j))**(0.5)

with open("data_output/statistics.txt",'w+') as fw:
	fw.write("Mean of all euclidean distances:"+str(mean_euc_dist)+"\n")
	fw.write("Maximum euclidean distance:"+str(max_dist)+"\n")
	fw.write("Maximum euclidean distance:"+str(sec_max_dist)+"\n")
	fw.write("Minimum euclidean distance:"+str(min_dist)+"\n")
	fw.write("Standard deviation of all euclidean distances:"+str(desv_padrao)+"\n")


fw.close()
print ("Success!\n")
