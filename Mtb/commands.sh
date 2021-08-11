# set up kallisto
kallisto index -i m_tb m_tuberculosis_H37Rv_NC_000962_cds.fa

# download data
#fasterq-dump SRR9042978
#fasterq-dump SRR9042979
#fasterq-dump SRR9042980
#
#gzip SRR9042978.fastq
#gzip SRR9042979.fastq
#gzip SRR9042980.fastq

kallisto quant --single -o wt_rep1 -i m_tb -t 6 -l 200 -s 25 SRR9042978.fastq.gz
kallisto quant --single -o wt_rep2 -i m_tb -t 6 -l 200 -s 25 SRR9042979.fastq.gz
kallisto quant --single -o wt_rep3 -i m_tb -t 6 -l 200 -s 25 SRR9042980.fastq.gz

# now merge tpm across replicates
paste wt_rep1/abundance.tsv wt_rep2/abundance.tsv wt_rep2/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10 + $15) /3}' > wt_average_tpm.txt

# get the genes above our tpm threshold
n_genes=$(( `cat wt_average_tpm.txt | wc -l` - 1 ))
n_thresh=$(( $n_genes / 5 ))

grep -v target_id wt_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_wt.txt

# set up the target list
grep -f ../good_terms m_tb_ATCC_25618.tab | grep -v -f ../bad_terms | awk -F "\t" '{print $21}' | sed s/";"/"\n"/g | grep "\S" > target_list.txt

# now get our real targets
grep -h -f target_list.txt high_exp_wt.txt | awk '{print $1}' | cut -c 21-29 > nap_candidates_refseq.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > nap_hits.txt
grep -f nap_candidates_refseq.txt m_tb_ATCC_25618.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> nap_hits.txt

# and merge this back in with the expression data
python ../add_expr_dat.py nap_hits.txt nap_hits_full.txt wt_average_tpm.txt

