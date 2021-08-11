# get the sequencing reads

for acc in SRR11998391 SRR11998392 SRR11998393  SRR11998412 SRR11998413 SRR11998414 ; do
	fasterq-dump $acc
    gzip ${acc}.fastq
done

# set up kallisto
kallisto index -i N_men GCF_000008805.1_ASM880v1_cds_from_genomic.fna

kallisto quant -o control_r1 -i N_men -t 6 --single -l 250 -s 75 SRR11998391.fastq.gz 
kallisto quant -o control_r2 -i N_men -t 6 --single -l 250 -s 75 SRR11998392.fastq.gz 
kallisto quant -o control_r3 -i N_men -t 6 --single -l 250 -s 75 SRR11998393.fastq.gz 
kallisto quant -o stat_r1 -i N_men -t 6 --single -l 250 -s 75 SRR11998412.fastq.gz 
kallisto quant -o stat_r2 -i N_men -t 6 --single -l 250 -s 75 SRR11998413.fastq.gz 
kallisto quant -o stat_r3 -i N_men -t 6 --single -l 250 -s 75 SRR11998414.fastq.gz 

# now merge tpm across replicates
paste control_r1/abundance.tsv control_r2/abundance.tsv control_r3/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10 + $15) /3}' > control_average_tpm.txt
paste stat_r1/abundance.tsv stat_r2/abundance.tsv stat_r3/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10 + $15) /3}' > stat_average_tpm.txt

# get the genes above our tpm threshold
n_genes=$(( `cat control_average_tpm.txt | wc -l` - 1 ))
n_thresh=$(( $n_genes / 5 ))

grep -v target_id control_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_control.txt
grep -v target_id stat_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_stat.txt

# set up the target list
grep -f ../good_terms Nm_proteome.tab | grep -v -f ../bad_terms | awk -F "\t" '{print $21}' | sed s/";"/"\n"/g | grep "\S" > target_list.txt

# now get our real targets
grep -h -f target_list.txt high_exp_control.txt high_exp_stat.txt | awk '{print $1}' | cut -c 21-34 > nap_candidates_refseq.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > nap_hits.txt
grep -f nap_candidates_refseq.txt Nm_proteome.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> nap_hits.txt               

python ../add_expr_dat.py nap_hits.txt nap_hits_test.txt control_average_tpm.txt
python ../add_expr_dat.py nap_hits_test.txt nap_hits_full.txt stat_average_tpm.txt

