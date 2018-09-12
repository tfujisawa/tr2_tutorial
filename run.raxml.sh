#! /bin/bash

for sq in "$@"
do
	#echo "raxml -m GTRGAMMA -T 2 -f a -p 12345 -s $sq -# 100 -x 100 -n ${sq//.fasta} -o Agkistrodon_contortrix_1,Agkistrodon_contortrix_2,Agkistrodon_piscivorus_1,Agkistrodon_piscivorus_2"
	raxmlHPC -m GTRGAMMA -T 2 -f a -p 12345 -s $sq -# 100 -x 100 -n ${sq//.fasta} -o Agkistrodon_contortrix_1,Agkistrodon_contortrix_2,Agkistrodon_piscivorus_1,Agkistrodon_piscivorus_2
done

