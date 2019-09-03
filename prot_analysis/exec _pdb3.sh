#!/bin/bash

echo "
=========================================================
        Federal University of Technology - Paraná 
=========================================================

---------------------------------------------------------
"
table=alberts
flag='0'
o_pa='('
c_pa=')'

while [ flag='0' ];
do
read  -p 'Want add a protein? '$o_pa'y/n'$c_pa'  ' var



if [  $var = 'y' ] && [ $flag==0 ] ;
then

read -p 'Protein name '$o_pa'lowercase'$c_pa':  ' protein
read -p 'Table of preference '$o_pa'lowercase'$c_pa':  ' table
python3 src/pdb3.py <<EOF
$protein
$table
EOF



elif [ $var = 'n' ];
then 

while read p
do

protein=$p
python3 src/pdb3.py <<EOF
$protein
$table
EOF
echo ""
done < data_input/proteins.txt
break 




else 
echo "Invalid command"

fi



done

