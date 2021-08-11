# NAPper
Utilities and examples for screening proteomes for potential NAPs using transcriptomics data

Each directory, Bs, Cc, Ec, Mtb, Nm, Pa, contains uniprot proteome files with GO term annotations and a file called `commands.sh` containing the commands that can be run to screen proteomes for potential NAPs.

Within each of the Bs, Cc, Ec, Mtb, Nm, and Pa directories, the final output of the screen is contained in the file `nap_hits_full.txt`. The R code in the src directory can be used to calculate the log2(TPM) from the TPM values in each file and to save the final results for each bacterium in this repository.

The file `Bs_results_category_supplemented.tsv` contains the *Bacillus subtilis* results with each protein's category supplemented as an extra column. The R code in the src directory can also be used to make plots comparing expression levels between the two conditions in the *B. subtilis* data we used.
