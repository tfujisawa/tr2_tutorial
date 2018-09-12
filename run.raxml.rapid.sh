#! /bin/bash
for sq in "$@"
do
	raxmlHPC -m GTRGAMMA -T 2 -f d -p 12345 -s $sq -n ${sq//.fasta} -o  Agkistrodon_contortrix_1,Agkistrodon_contortrix_2,Agkistrodon_piscivorus_1,Agkistrodon_piscivorus_2
done
