# set up kallisto
kallisto index -i na1000 na1000_cds.fa

# get the data
#fasterq-dump SRR5772180
#fasterq-dump SRR5772181
#
#gzip SRR5772180.fastq
#gzip SRR5772181.fastq

kallisto quant --single -o wt -i na1000 -t 6 -l 200 -s 25 SRR5772180.fastq.gz
kallisto quant --single -o gapR_dep -i na1000 -t 6 -l 200 -s 25 SRR5772181.fastq.gz

# now merge tpm across replicate
cat wt/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, $5}' > wt_tpm.txt
cat gapR_dep/abundance.tsv  | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, $5}' > gapR_dep_tpm.txt

# get the genes above our tpm threshold
n_genes=$(( `cat wt_tpm.txt | wc -l` - 1 ))
n_thresh=$(( $n_genes / 5 ))

grep -v target_id wt_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_wt.txt
grep -v target_id gapR_dep_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_gapR_dep.txt

# set up the target list
grep -f ../good_terms c_crescentus_na1000.tab | grep -v -f ../bad_terms | awk -F "\t" '{print $21}' | sed s/";"/"\n"/g | grep "\S" > target_list.txt

# now get our real targets
grep -h -f target_list.txt high_exp_wt.txt high_exp_gapR_dep.txt | awk '{print $1}' | cut -c 21-31 > nap_candidates_refseq.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > wt_nap_hits.txt
grep -f nap_candidates_refseq.txt c_crescentus_na1000.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> wt_nap_hits.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > gapR_dep_nap_hits.txt
grep -f nap_candidates_refseq.txt c_crescentus_na1000.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> gapR_dep_nap_hits.txt

# and merge this back in with the expression data
python ../add_expr_dat.py wt_nap_hits.txt nap_hits_test.txt wt_tpm.txt
python ../add_expr_dat.py nap_hits_test.txt nap_hits_full.txt gapR_dep_tpm.txt

