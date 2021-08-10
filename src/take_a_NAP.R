library(tidyverse)
srcdir = 'src'
source(file.path(srcdir,'helper_functions.R'))

bact_direcs = c(
  'Bs', # Bacillus subtilis
  'Mtb', # Mycobacterium tuberculosis
  'Ec', # Escherichia coli
  'Pa', # Pseudomonas aeruginosa
  'Nm', # Nisseria meningitidis
  'Cc' # Caulobacter crescentus
)

data_tib = tibble(direc = bact_direcs) %>%
  mutate(hits = purrr::map(file.path(bact_direcs, 'nap_hits_full.txt'), read_tsv))

(data_tib %>%
  filter(direc=="Bs") %>%
  .$hits)[[1]] %>%
  mutate(logTPM_8h = log2(`8h_average_tpm.txt`),
         logTPM_24h = log2(`24h_average_tpm.txt`)) %>%
  write_tsv(file="Bs_results.txt")

(data_tib %>%
  filter(direc=="Ec") %>%
  .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  write_tsv(file="Ec_results.txt")

(data_tib %>%
    filter(direc=="Nm") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  write_tsv(file="Nm_results.txt")

(data_tib %>%
    filter(direc=="Pa") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  write_tsv(file="Pa_results.txt")

(data_tib %>%
    filter(direc=="Cc") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_wt = log2(`wt_tpm.txt`),
         logTPM_gapR_depleted = log2(`gapR_dep_tpm.txt`)) %>%
  write_tsv(file="Cc_results.txt")

(data_tib %>%
    filter(direc=="Mtb") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_wt = log2(`wt_average_tpm.txt`)) %>%
  write_tsv(file="Mtb_results.txt")
