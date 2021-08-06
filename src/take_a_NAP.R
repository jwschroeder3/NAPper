library(tidyverse)
source('helper_functions.R')

main_dir = '/home/jeremy/src/take_a_nap/transcriptomics'
# vector to iterate over the directories I have within main_dir
bact_direcs = c(
  'b_subtilis',
  #     's_pyogenes',
  'm_tb',
  'c_crescentus'
)

# set up initial tibble with file names and use purrr::map to
#   read in the kallisto-generated output for each SRR run
tib = tibble(
    bugs = bact_direcs,
    dir = paste(main_dir, bugs, sep='/'),
    srr_dir = purrr::map(dir, get_srr_dirs),
  ) %>%
  tidyr::unnest(srr_dir) %>%
  mutate(
    srr = purrr::map_chr(srr_dir, grab_srr),
    fname = paste(srr_dir, "abundance.tsv", sep='/'),
    data = purrr::map(fname, read_delim, delim='\t')
  )

# read in tables with each dataset's run information
#   we can use the info here to filter just the 
#   runs we want
bsub_sample_info = read_csv('b_sub_sample_info.csv')
# keep the 8-hour and 24-hour timepoints from Zizi's work
bsub_soi = bsub_sample_info %>%
  filter(condition %in% c('8h', '24h'), replicate == 1)

cc_sample_info = read_csv('c_crescentus_sample_info.csv')
# keep just WT and GapR-depleted RNA-seq runs from Monica's work
cc_soi = cc_sample_info %>%
  filter(condition %in% c('WT glucose 6h', 'GapR depletion glucose 6h'))

mtb_sample_info = read_csv('m_tb_sample_info.csv')
# keep just the wild type data from Christina Stallings' lab's work
mtb_soi = mtb_sample_info %>%
  filter(condition == "WT carD", replicate == 1)

all_run_info = bind_rows(bsub_sample_info, cc_sample_info, mtb_sample_info)

# place run ids of interest into a single vector
srrs_of_interest = c(bsub_soi$SRR, cc_soi$SRR, mtb_soi$SRR)
srrs_of_interest

# keep only the runs of interest
filtered_tib = tib %>%
  filter(srr %in% srrs_of_interest)

# these files were downloaded from uniprot
#   rename the crucial column for convenience
bsub_uniprot = read_delim(
    'uniprot_refs/b_subtilis_UP000001570.tab',delim='\t'
  ) %>%
  mutate(bug = "B. subtilis") %>%
  rename(
    locus_tag = `Gene names  (ordered locus )`,
    gene = `Gene names  (primary )`,
    protein = `Protein names`,
    GO_terms = `Gene ontology IDs`
  )
cc_uniprot = read_delim(
    'uniprot_refs/c_crescentus_na1000.tab',delim='\t'
  ) %>%
  mutate(bug = "C. crescentus") %>%
  rename(
    locus_tag = `Gene names  (ordered locus )`,
    gene = `Gene names  (primary )`,
    protein = `Protein names`,
    GO_terms = `Gene ontology IDs`
  )
mtb_uniprot = read_delim(
    'uniprot_refs/m_tb_ATCC_25618.tab',delim='\t'
  ) %>%
  mutate(bug = "M. tuberculosis") %>%
  rename(
    locus_tag = `Gene names  (ordered locus )`,
    gene = `Gene names  (primary )`,
    protein = `Protein names`,
    GO_terms = `Gene ontology IDs`
  )
# throw all this info into a single tibble now
uniprot_tib = bind_rows(bsub_uniprot, cc_uniprot, mtb_uniprot)

# read in lookup tables for associating protein IDs with gene locus tags
bsub_lut = read_delim(
    'b_sub_protID_locus_lut.tsv',
    delim = '\t',
    col_names = c("prot_id","locus_tag")
  ) %>%
  mutate(bug = "B. subtilis")
cc_lut = read_delim(
    'c_crescentus_protID_locus_lut.tsv',
    delim = '\t',
    col_names = c("prot_id","locus_tag")
  ) %>%
  mutate(bug = "C. crescentus")
mtb_lut = read_delim(
    'm_tb_protID_locus_lut.tsv',
    delim = '\t',
    col_names = c("prot_id","locus_tag")
  ) %>%
  mutate(bug = "M. tuberculosis")
# put 'em all into one tibble
lut = bind_rows(bsub_lut, cc_lut, mtb_lut)

expanded_tib = filtered_tib %>%
  mutate(
    # calculate percentiles
    percentiles = purrr::map(data, calc_tpm_percentile),
    # add protein IDs
    prot_id = purrr::map(data, extract_prot_id)
  ) %>%
  # unravel it all into a long tibble, rather than a nested one
  unnest(cols=c(data,percentiles,prot_id)) %>%
  # finally, join the lookup table up to the tibble
  left_join(lut, by="prot_id") %>%
  left_join(all_run_info %>% dplyr::select(-bug), by=c("srr"="SRR"))

# so at this point, I have to remove the "_" from the B sub locus tags
#   since older tags have BSU_XXXXX, and newer ones or BSUXXXXX
new_tags = NA
for (i in 1:nrow(expanded_tib)) {
  tag = expanded_tib$locus_tag[i]
  if (grepl("^BSU", tag)) {
    new_tags[i] = str_replace(tag, "_", "")
  } else {
    new_tags[i] = tag
  }
}

expanded_tib = expanded_tib %>%
  mutate(locus_tag = new_tags)

with_prot_info = expanded_tib %>%
  # join the crucial uniprot information to this tibble
  left_join(
    uniprot_tib %>% dplyr::select(
      protein,
      gene,
      locus_tag,
      GO_terms
    ),
    by="locus_tag"
  )

bsub_test_naps = c(
  'rok',
  'hupA', # uniprot's annotated name for hbs
  'lrpC',
  'abrB'
)

with_prot_info %>%
  filter(
    bugs == "b_subtilis",
    gene %in% bsub_test_naps
  ) %>%
  arrange(gene)

DNA_binding_term = "GO:0003677"
NOT_terms = c(
  "GO:0003887", # DNA-directed DNA polymerase activity
  "GO:0016987", # Sigma factor activity
  "GO:0006412", # Translation
  "GO:0003899", # DNA-directed 5'-3' RNA polymerase activity
  "GO:0006354", # DNA-templated transcription, elongation
  "GO:0006261", # DNA-dependent DNA replication
  "GO:0003896" # DNA primase activity
)

# the filter_NAPs function comes from helper_functions.R
#   it joins the uniprot info onto the first input tibble
#   then it filters GO terms and gene expression percentile
DNA_binders = filter_NAPs(
  expanded_tib,
  uniprot_tib,
  DNA_binding_term,
  NOT_terms,
  perc_threshold = 0
)

nested_binders = DNA_binders %>%
  tidyr::nest(data = c(locus_tag, protein, gene, GO_terms, percentiles, condition))

wide_nested_binders = nested_binders %>%
  mutate(
    wide_data = purrr::map(data, purr_widen, names=condition, values=percentiles),
    thing = str_replace(bug, ". ", "_"),
    out_fname = paste0(thing, "_binders.tsv")
  ) %>%
  dplyr::select(-c(data,thing))

mapply(
  function(.x, .y) write_delim(x=.x, file=.y, delim='\t'),
  wide_nested_naps$wide_data,
  wide_nested_naps$out_fname
)

potential_naps = filter_NAPs(
  expanded_tib,
  uniprot_tib,
  DNA_binding_term,
  NOT_terms,
  perc_threshold = 80
)

nested_naps = potential_naps %>%
  tidyr::nest(data = c(locus_tag, protein, gene, GO_terms, percentiles, condition))

wide_nested_naps = nested_naps %>%
  mutate(
    wide_data = purrr::map(data, purr_widen, names=condition, values=percentiles),
    thing = str_replace(bug, ". ", "_"),
    out_fname = paste0(thing, "_potential_NAPs.tsv")
  ) %>%
  dplyr::select(-c(data,thing))

mapply(
  function(.x, .y) write_delim(x=.x, file=.y, delim='\t'),
  wide_nested_naps$wide_data,
  wide_nested_naps$out_fname
)

#---- B. subtilis
bsub_binders = wide_nested_binders %>%
  filter(bug == "B. subtilis") %>%
  unnest(cols=wide_data)

p = bsub_binders %>%
  ggplot(aes(x = `8h`, y = `24h`, label = gene, color=ifelse((`8h`>=80 | `24h`>=80), 'r', 'b'))) +
  geom_point() +
  scale_color_manual(name="Potential NAP", labels=c("Yes","No"), breaks=c('r','b'), values=c("red","#999999")) +
  labs(y="24 hours in biofilm developemnt", x="8 hours in biofilm development")
p
ggplotly(p)


#---- C crescentus
cc_binders = wide_nested_binders %>%
  filter(bug == "C. crescentus") %>%
  unnest(cols=wide_data)

p = cc_binders %>%
  ggplot(aes(x = `WT glucose 6h`, y = `GapR depletion glucose 6h`, label = gene, color=ifelse((`WT glucose 6h`>=80 | `GapR depletion glucose 6h`>=80), 'r', 'b'))) +
  geom_point() +
  scale_color_manual(name="Potential NAP", labels=c("Yes","No"), breaks=c('r','b'), values=c("red","#999999")) +
  labs(y="GapR depletion", x="WT")
p
ggplotly(p)

#---- M tb
# mtb_binders = wide_nested_binders %>%
#   filter(bug == "M. tuberculosis") %>%
#   unnest(cols=wide_data)
# 
# p = mtb_binders %>%
#   ggplot(aes(y=`WT carD`, x=bug, label = gene, color=ifelse(`WT carD`>=80, 'r', 'b'))) +
#   # geom_boxplot() +
#   geom_point(position="jitter") +
#   scale_color_manual(name="Potential NAP", labels=c("Yes","No"), breaks=c('r','b'), values=c("red","#999999")) +
#   labs(y="Expression percentile", x="Mtb expression")
# p
# ggplotly(p)
