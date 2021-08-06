# set up kallisto
kallisto index -i Ec_mg1655 GCF_000005845.2_ASM584v2_cds_from_genomic.fna

kallisto quant -o logphase_R1 -i Ec_mg1655 -t 12 glu_l_rna_rep1/glu_l_rna_rep1_R1.fq.gz glu_l_rna_rep1/glu_l_rna_rep1_R2.fq.gz
kallisto quant -o logphase_R2 -i Ec_mg1655 -t 12 glu_l_rna_rep2/glu_l_rna_rep2_R1.fq.gz glu_l_rna_rep2/glu_l_rna_rep2_R2.fq.gz
kallisto quant -o statphase_R1 -i Ec_mg1655 -t 12 glu_s2_rna_rep1/glu_s2_rna_rep1_R1.fq.gz glu_s2_rna_rep1/glu_s2_rna_rep1_R2.fq.gz
kallisto quant -o statphase_R2 -i Ec_mg1655 -t 12 glu_s2_rna_rep2/glu_s2_rna_rep2_R1.fq.gz glu_s2_rna_rep2/glu_s2_rna_rep2_R2.fq.gz


# now merge tpm across replicates
paste logphase_R1/abundance.tsv logphase_R2/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10) /2}' > control_average_tpm.txt
paste statphase_R1/abundance.tsv statphase_R2/abundance.tsv  | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10) /2}' > stat_average_tpm.txt

# get the genes above our tpm threshold
n_genes=$(( `cat control_average_tpm.txt | wc -l` - 1 ))
n_thresh=$(( $n_genes / 5 ))

grep -v target_id control_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_control.txt
grep -v target_id stat_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_stat.txt

# set up the target list
grep -f ../good_terms mg1655_proteome.tab | grep -v -f ../bad_terms | awk -F "\t" '{print $21}' | sed s/";"/"\n"/g | grep "\S" > target_list.txt

# now get our real targets
grep -h -f target_list.txt high_exp_control.txt high_exp_stat.txt | awk '{print $1}' | cut -c 21-31 > nap_candidates_refseq.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > nap_hits.txt
grep -f nap_candidates_refseq.txt mg1655_proteome.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> nap_hits.txt

# and merge this back in with the expression data
python ../add_expr_dat.py nap_hits.txt nap_hits_test.txt control_average_tpm.txt
python ../add_expr_dat.py nap_hits_test.txt nap_hits_full.txt stat_average_tpm.txt



