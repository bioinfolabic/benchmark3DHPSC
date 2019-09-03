
pdb3 - input 

affinity hydrophilic scale table, exemple : alberts.txt



pdb3 - output

proteins.txt
"protein".ciff
seq_"protein".txt
"protein"_A.pdb
"protein"_CA.txt
seq_10_"protein".txt
seq_ab_"protein".txt
posicao_com_numeracao_"protein".txt
"protein"_pos_CA_h.txt

##########################################################################

calc_statistcs.py - input

proteins.txt
"protein"_pos_CA_h.txt



calc_statistcs.py - output

dist_euc_all_prot.txt
"protein"_euc_dist.txt
all_euc_dist_asc.txt
statistcs.txt

###########################################################################

ver_pos_score.py - input

pred.csv
proteins.txt
"protein"_euc_dist.txt
output_protein_data_"protein".txt
"protein"_CA_euc_map.txt



ver_pos_score.py - output

protein_matchs.csv
protein_f1scores.csv
protein_recall.csv
protein_precision.csv
Dists_x_Score.csv


############################################################################

create_3d_model - input

proteins.txt
output_protein_data_"protein".txt



create_3d_model - input

coordenadas_rasmol_"protein".pdb










