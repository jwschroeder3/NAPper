## set up kallisto
kallisto index -i bsub_168 bsub_168_NC_000964_cds.fa

# download data and run kallisto
for cond in 8h 24h; do
    echo $cond
    while IFS= read -r line; do
        acc=$(echo $line | cut -f1 -d' ')
        fasterq-dump $acc
        gzip ${acc}.fastq
        cond=$(echo $line | cut -f2 -d' ')
        rep=$(echo $line | cut -f3 -d' ')
        kallisto quant --single -o "${cond}_rep${rep}" -i bsub_168 \
            -t 6 -l 200 -s 25 ${acc}.fastq.gz
    done  < <(grep "${cond}" bsub_sample_info.txt)
done

# now merge tpm across replicates
paste 8h_rep1/abundance.tsv 8h_rep2/abundance.tsv 8h_rep3/abundance.tsv 8h_rep4/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10 + $15 + $20) /4}' > 8h_average_tpm.txt
paste 24h_rep1/abundance.tsv 24h_rep2/abundance.tsv 24h_rep3/abundance.tsv 24h_rep4/abundance.tsv | awk 'BEGIN {OFS="\t"}; NR==1 {print $1, $5}; NR>1 {print $1, ($5 + $10 + $15 + $20) /4}' > 24h_average_tpm.txt

# get the genes above our tpm threshold
n_genes=$(( `cat 8h_average_tpm.txt | wc -l` - 1 ))
n_thresh=$(( $n_genes / 5 ))

grep -v target_id 8h_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_8h.txt
grep -v target_id 24h_average_tpm.txt | sort -g --key=2 | tail -n $n_thresh > high_exp_24h.txt

# set up the target list
grep -f ../good_terms b_subtilis_UP000001570.tab | grep -v -f ../bad_terms | awk -F "\t" '{print $21}' | sed s/";"/"\n"/g | grep "\S" > target_list.txt

# now get our real targets
grep -h -f target_list.txt high_exp_8h.txt high_exp_24h.txt | awk '{print $1}' | cut -c 21-31 > union_nap_candidates_refseq.txt
echo "uniprot_id;protein;gene;refseq" | sed s/";"/"\t"/g > union_nap_hits.txt
grep -f union_nap_candidates_refseq.txt b_subtilis_UP000001570.tab | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1, $4, $10, $21}' >> union_nap_hits.txt

# and merge this back in with the expression data
python ../add_expr_dat.py union_nap_hits.txt nap_hits_test.txt 8h_average_tpm.txt
python ../add_expr_dat.py nap_hits_test.txt nap_hits_full.txt 24h_average_tpm.txt

