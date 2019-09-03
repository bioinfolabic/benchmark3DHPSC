"""
=========================================================
        FEDERAL UNIVERSITY OF TECHNOLOGY - PARANA        
=========================================================

   NOTES:
---------------------------------------------------------
   
   The pdb files are in a format accept by Rasmol

---------------------------------------------------------
"""

import re 
import numpy as np
from scipy.spatial import distance

print ("Initiating 3d model files...")


with open("data_input/proteins.txt",'r') as fp:
	 
	for prot in fp:
		protein=prot.strip()

	
### Patters to extraction of data  ###
		pattern  = re.compile(r'(?:xiv)+(\d+\D+\d+)')
		pattern2 = re.compile(r'(?:xjw)+(\d+\D+\d+)')
		pattern3 = re.compile(r'(?:yiw)+(\d+\D+\d+)')
		pattern4 = re.compile(r'(?:yjv)+(\d+\D+\d+)')
		pattern5 = re.compile(r'(?:^n\s+\D)+(\s+\d+)')
		pattern6 = re.compile(r'(?:hh)+(\d+\D+\d+\D+\d+\D+\d+)')

		bb_size=[]### backbone  ###
		sc_hh=[]### The id of the hydrophobics sidechains  ###
		flag=0
		x=0
		y=0
		z=0
		n=0
		m=0
		error_flag=0
###  Create the 3d model to help on vizualization  ###


		space_3d=np.zeros((5,5,5))


		nmrca_posa_posn=[]###  Save the data of each carbon alpha  ###
		nmrsc_posa_posn=[]###  Save the data of each sidechain  ###
		sc_contact=[]### Save the contact between the sidechans  ###

		try:
			with open("data_input/results-IP/output_protein_data_"+protein.upper()+".txt",'r') as fr:
				for fo in fr:
					x=0
					y=0
					z=0


###  Get the protein's size  ###
					if (flag==0):
						bb_size=str(re.findall(pattern5,fo)).replace("[","").replace("]","").replace("'","").replace("'","").strip()
						flag+=1


###  Get the hydrophobics sidechains  ####
					if(re.search("He:", fo)):
						sc_hh.append(fo.replace("He:","").split())

					if(re.search("Ho", fo)):
						sc_hh.append(fo.replace("Ho:","").split())

	###  Calcule and save the carbon alpha position in the 3d model,marked as "1" ###
					for found in re.findall(pattern,fo):
						aux=int((found.replace(","," ").split())[1])						
						while aux>=25:
							aux-=25
							x+=1
						while aux>=5:
							aux-=5
							y+=1				
						z=aux
						
						space_3d[x][y][z]=1

						nmrca_posa_posn.append(found.replace(","," ").strip()+" "+str(space_3d[x][y][z])+" x"+str(x)+" y"+str(y)+" z"+str(z))

						error_flag+=1
						n+=1
						continue

					for found in re.findall(pattern2,fo):
						aux=int((found.replace(","," ").split())[1])						
						while aux>=25:
							aux-=25
							x+=1
						while aux>=5:
							aux-=5
							y+=1				
						z=aux

						space_3d[x][y][z]=1

						nmrca_posa_posn.append(found.replace(","," ").strip()+" "+str(space_3d[x][y][z])+" x"+str(x)+" y"+str(y)+" z"+str(z))

						error_flag+=1
						n+=1
						continue


###  Calcule and save the sidechains position in the 3d model,marking hydrophobics as "2" and polar as "3" ###

					for found in re.findall(pattern3,fo):
						aux=int((found.replace(","," ").split())[1])						
						aux2=((found.replace(","," ").split())[0])						

						while aux>=25:
							aux-=25
							x+=1
						while aux>=5:
							aux-=5
							y+=1				
						z=aux

						
						if aux2 in (sc_hh[0]):
							space_3d[x][y][z]=2
						elif aux2 in (sc_hh[1]):
							space_3d[x][y][z]=2
						else:
							space_3d[x][y][z]=3

						nmrsc_posa_posn.append(found.replace(","," ").strip()+" "+str(space_3d[x][y][z])+" x"+str(x)+" y"+str(y)+" z"+str(z))

						error_flag+=1
						m+=1	
						continue


					for found in re.findall(pattern4,fo):
						aux=int((found.replace(","," ").split())[1])						
						aux2=((found.replace(","," ").split())[0])						

						while aux>=25:
							aux-=25
							x+=1
						while aux>=5:
							aux-=5
							y+=1				
						z=aux
						
						if aux2 in (sc_hh[0]):
							space_3d[x][y][z]=2
						elif aux2 in(sc_hh[1]):
							space_3d[x][y][z]=2
						else:
							space_3d[x][y][z]=3

						nmrsc_posa_posn.append(found.replace(","," ").strip()+" "+str(space_3d[x][y][z])+" x"+str(x)+" y"+str(y)+" z"+str(z))

						error_flag+=1
						m+=1		
						continue

					for found in re.findall(pattern6,fo):
						aux=found.replace(","," ").split()
						sc_contact.append(aux[0]+" "+aux[2])		

				fr.close()			
		except IOError:
			print (u'File of '+protein+' protein not found')
			continue

		if(error_flag==0):
			print('The file of '+protein+' protein is incomplete')
			continue

###  Convert the hydrophobics sidechains list in a int array  ###
		sc_hh=np.concatenate((sc_hh[1],sc_hh[0]))
		sc_hh=np.array(sc_hh,dtype=int)
		sc_hh=np.sort(sc_hh)

		aux1=[]
		aux2=[]
		i=0
		for i in range(len(nmrsc_posa_posn)):
			aux1.append(nmrca_posa_posn[i].split())
			aux2.append(nmrsc_posa_posn[i].split())
			aux1[i][0]=int(aux1[i][0])
			aux2[i][0]=int(aux2[i][0])

		aux1.sort()
		aux2.sort()


###  Process to format the string into pdb format  ###
		def format_atom_cd(i, atom, residue,n_residue, x_vertic, y_vertic, z_vertic,atom_min):
		    return (f'{"ATOM": <4}{i: >7}{"  ": >2}{atom: <4}{residue: <4}{"A":<1}{n_residue: >4}{x_vertic: >12.3f}{y_vertic: >8.3f}{z_vertic: >8.3f}{"  1.00  0.00": >11}{atom_min: >12}')
		def bond_atom(atom_a,atom_b,pula_linha):
		    return (f'{"CONECT": <6}{atom_a: >5}{atom_b: >5}{pula_linha: >1}')

		i=0
		n=3
	
		with open("proteins/"+protein+"/coordenadas_rasmol_"+protein+".pdb",'w+') as fw:
			for i in range(len(nmrsc_posa_posn)):
				i_atom=aux1[i][0]*2+1
				if (aux2[i][2]=="2.0"):
					ca=format_atom_cd( i_atom , "CA", "GLY" , i+1 , n*float(aux1[i][3].replace("x","")), n*float(aux1[i][4].replace("y","")), n*float(aux1[i][5].replace("z","")), "C\n")
					#h=format_atom_cd( i_atom , "H", "GLY" , i+1 , n*float(aux2[i][3].replace("x","")), n*float(aux2[i][4].replace("y","")), n*float(aux2[i][5].replace("z","")), "H\n")			
					sc=format_atom_cd( i_atom+1 , "O", "GLY" , i+1 , n*float(aux2[i][3].replace("x","")), n*float(aux2[i][4].replace("y","")), n*float(aux2[i][5].replace("z","")), "O\n")
				elif(aux2[i][2]=="3.0"):
					ca=format_atom_cd( i_atom , "CA", "LYS" , i+1 , n*float(aux1[i][3].replace("x","")), n*float(aux1[i][4].replace("y","")), n*float(aux1[i][5].replace("z","")), "C\n")
					#h=format_atom_cd( i_atom , "H", "LYS" , i+1 , n*float(aux2[i][3].replace("x","")), n*float(aux2[i][4].replace("y","")), n*float(aux2[i][5].replace("z","")), "H\n")			
					sc=format_atom_cd( i_atom+1 , "N", "LYS" , i+1 , n*float(aux2[i][3].replace("x","")), n*float(aux2[i][4].replace("y","")), n*float(aux2[i][5].replace("z","")), "N\n")
				else:
					print ("ERRO CRIACAO PDB")
					ERRO
					break

				fw.write(ca)
				fw.write(sc)

###  Create the bonds of the protein  ###	
			for i in range(1,len(nmrca_posa_posn)*2,2): 
				bd_at=bond_atom(i,i+2,"\n")
				fw.write(bd_at)
				bd_at=bond_atom(i,i+1,"\n")
				fw.write(bd_at)
			fw.write("END")
			fw.close()




fp.close()

print("Success!\n")