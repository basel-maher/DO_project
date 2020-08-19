#!/bin/bash


i=$SLURM_ARRAY_TASK_ID
tis_name=$(sed -n "${i}{p;}" < /scratch/bma8ne/GTEx_v7/coloc_BAN/src/tis_supersplit_v7)

echo $tis_name

fileOut=/scratch/bma8ne/GTEx_v7/coloc_gefos_BAN/eqtl_200k/${tis_name}_200k.txt


while read rsid chr pos id_19 ens symbol; do awk -v "fileOut=$fileOut" -v "tis_name=$tis_name" -v "ens=$ens" -v "pos=$pos" -v "chr=$chr" -v "rsid=$rsid" -v "symbol=$symbol" -v "id_19=$id_19" -F ':|\t|_' ' function abs(v) {return v < 0 ? -v : v} (ens == substr($1,1,15) && chr == $2 && abs($3) <= abs(pos) + 200000 && abs($3) >= abs(pos) - 200000) {print $0, rsid, chr, pos, ens, symbol, id_19 >> fileOut}' /scratch/bma8ne/GTEx_v7/coloc_BAN/split/${tis_name}; done < /scratch/bma8ne/GTEx_v7/coloc_gefos_BAN/src/lead_merged
