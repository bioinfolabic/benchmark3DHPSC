"""
=========================================================
        FEDERAL UNIVERSITY OF TECHNOLOGY - PARANA        
=========================================================


---------------------------------------------------------
"""

import re 
import numpy as np
import csv
from scipy.spatial import distance

print ("Initiating positions and score...")



### Initiate the csv files where the data is going to be saved  ##
c = csv.writer(open("data_output/protein_matchs.csv", "w+"),delimiter=";")
p = csv.writer(open("data_output/protein_f1scores.csv", "w+"),delimiter=";")
rec=csv.writer(open("data_output/protein_recall.csv", "w+"),delimiter=";")
precs= csv.writer(open("data_output/protein_precision.csv", "w+"),delimiter=";")


###  Csv file that will be read  ###
arquivo = open('data_input/pred.csv','r')
linhas = csv.reader(arquivo,delimiter=";")

###  Get data from .csv  ###
data=[]
for linha in linhas:
	data.append((str(linha[0].lower())+" "+str(linha[2])).split())

c_dxs = csv.writer(open("data_output/Dists_x_Score.csv", "w+"))
c_dxs.writerow(["Distance;Score"])


flag1=0
flag2=0
n_values=[]
n_total_fs=[]
n_total_prcs=[]
n_total_rc=[]


###  Travel through the obtained proteins and the range of threshold desired  ###
with open("data_input/proteins.txt",'r') as fr:

	for prot in fr:
		n=3.8
		n_matchs=[]
		n_fs=[]
		n_prcs=[]
		n_rc=[]
		
		while n<=9.5:
			
			protein=prot.strip()
			hh=[]
			score=0

###  Read the csv file e get the score  ###
			ec_d1=[]
			with open("proteins/"+protein+"/"+protein+"_euc_dist.txt", 'r') as fo:
				for dists in fo:
					ec_d1.append(dists.strip().split())		
				aux1=np.array(ec_d1,dtype=np.float64)
			fo.close()

			for i in range(len(aux1)):
				for j in range(len(aux1)):
					if (aux1[i][j]>n):
						aux1[i][j]=0
						
					if (aux1[i][j]==aux1[j][i]):
						aux1[j][i]=0
						
			hh_count=np.count_nonzero(aux1)



			for i in range(len(data)):
				if data[i][0]==protein:
					break;
			score+=(hh_count-float(data[i][1]))
			
			c_dxs.writerow([str(n)+";"+str(score)])

			i=0
			j=0

###  Get the number of contacts in the predict model  ###
			try:
				with open("data_input/results-IP/output_protein_data_"+protein.upper()+".txt",'r') as fr2:
					for aux_r in fr2:
						if(re.search(r'^hh', aux_r)):
							aux_r=aux_r.replace("hh","")
							aux_r=aux_r.replace(","," ")
							aux_r=aux_r.strip().split()
							if int(aux_r[0])<int(aux_r[2]):
								hh.append("["+aux_r[0]+"]["+aux_r[2]+"]")
							else:
								hh.append("["+aux_r[2]+"]["+aux_r[0]+"]")
				hh_len=len(hh)
				hh=np.array(hh)
				fr2.close()	
			
			except IOError:
				print (u'File of '+protein+' protein not found.')
				n+=100
				continue

###  Get the number of contacts for the threshold in the real model  ###
			with open("proteins/"+protein+"/"+protein+"_CA_euc_map.txt", 'r') as fa:
				ec_d=[]
				num=[]
				for aux_n in fa:
					num.append(aux_n.strip().split())
			fa.close()
			
			with open("proteins/"+protein+"/"+protein+"_euc_dist.txt", 'r') as fo:
				for dists in fo:
					ec_d.append(dists.strip().split())		
			fo.close()
			
			aux=np.array(ec_d,dtype=np.float64)
			for i in range(len(aux)):
				for j in range(len(aux)):
					if (aux[i][j]>n):
						aux[i][j]=0
						num[i][j]=0
					if (aux[i][j]==aux[j][i]):
						aux[j][i]=0
						num[j][i]=0
			a=np.array(num)
			a=a.flatten()
			if(np.count_nonzero(a)!=0):
				a=a[a!='0']				


###  Get the intersction between the predict and real models  ###
			its=np.intersect1d(hh,a)


###  Get a specif number in the range  ###
### 		if str(n)=="6.6": 
###  Get the F1 Score,Precision and Recall  ###
			tp=float(len(its))   			##True positive
			fn=float(len(a)-len(its))		##False Negative
			fp=float((hh_len-len(its)))		##False positive
			
			fs=(((2*tp)/(2*tp + fp + fn)))

			precision=tp/(tp+fp)

			recall=tp/(tp+fn)


			matchs=len(its)
			n_matchs.append(str(matchs))
			n_fs.append(str(fs))
			n_prcs.append(str(precision))
			n_rc.append(str(recall))

			if flag1==0:
				n_values.append(str(n))

			if str(n)=="9.4":
				n=9.5
			else:
				n+=0.2


###  Flag to go to the next protein  ###
		if n>100:
			continue


		for i in range(len(n_values)):
			if flag2==0:
				n_total_fs.append(0)
				n_total_prcs.append(0)
				n_total_rc.append(0)
			n_total_fs[i]+=float(n_fs[i])
			n_total_rc[i]+=float(n_rc[i])
			n_total_prcs[i]+=float(n_prcs[i])

###  Write into the csv file  ### 
		if flag2==0:
			c.writerow(["Protein"]+n_values)
			p.writerow(["Protein"]+n_values)
			rec.writerow(["Protein"]+n_values)
			precs.writerow(["Protein"]+n_values)			
			flag2+=1
			flag1+=1

		c.writerow([protein]+n_matchs)
		p.writerow([protein]+n_fs)
		rec.writerow([protein]+n_rc)
		precs.writerow([protein]+n_prcs)			

	for i in range(len(n_values)):
		n_total_fs[i]=str(n_total_fs[i])
		n_total_rc[i]=str(n_total_rc[i])
		n_total_prcs[i]=str(n_total_prcs[i])

	p.writerow(["Total:"]+n_total_fs)
	rec.writerow(["Total:"]+n_total_rc)
	precs.writerow(["Total:"]+n_total_prcs)


fr.close()

print("Success!\n")
# Bom range 3.8  - 9.5