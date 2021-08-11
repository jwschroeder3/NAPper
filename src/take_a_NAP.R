library(tidyverse)
library(ggrepel)

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
  arrange(desc(logTPM_8h)) %>%
  write_tsv(file="Bs_results.txt")

(data_tib %>%
  filter(direc=="Ec") %>%
  .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  arrange(desc(logTPM_control)) %>%
  write_tsv(file="Ec_results.txt")

(data_tib %>%
    filter(direc=="Nm") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  arrange(desc(logTPM_control)) %>%
  write_tsv(file="Nm_results.txt")

(data_tib %>%
    filter(direc=="Pa") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_control = log2(`control_average_tpm.txt`),
         logTPM_stat = log2(`stat_average_tpm.txt`)) %>%
  arrange(desc(logTPM_control)) %>%
  write_tsv(file="Pa_results.txt")

(data_tib %>%
    filter(direc=="Cc") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_wt = log2(`wt_tpm.txt`),
         logTPM_gapR_depleted = log2(`gapR_dep_tpm.txt`)) %>%
  arrange(desc(logTPM_wt)) %>%
  write_tsv(file="Cc_results.txt")

(data_tib %>%
    filter(direc=="Mtb") %>%
    .$hits)[[1]] %>%
  mutate(logTPM_wt = log2(`wt_average_tpm.txt`)) %>%
  arrange(desc(logTPM_wt)) %>%
  write_tsv(file="Mtb_results.txt")

#---- plot results

plot_data = read_tsv('Bs_results_category_supplemented.tsv')

exp_data_8h = read_tsv('Bs/8h_average_tpm.txt')
perc_8h = log2(quantile(exp_data_8h$tpm, c(0.8,0.9)))
exp_data_24h = read_tsv('Bs/24h_average_tpm.txt')
perc_24h = log2(quantile(exp_data_24h$tpm, c(0.8,0.9)))

xmin = 7.3
xmax = 9.7
ymin = 7.8
ymax = 9.1

p = plot_data %>%
  ggplot() +
  geom_vline(xintercept=perc_24h, linetype='dashed', color='red') +
  geom_hline(yintercept=perc_8h, linetype='dashed', color='red') +
  geom_point(aes(x=logTPM_24h, y=logTPM_8h, color=category)) +
  scale_color_discrete(name="Category", breaks=c("GR","LR","NAP","PNAP","Other"), labels=c("Global regulator", "Local regulator", "NAP", "Potential NAP", "Other")) +
  geom_text_repel(aes(x=logTPM_24h, y=logTPM_8h, label=gene), color='black', fontface='italic') +
  geom_rect(data=tibble(xl=xmin,xu=xmax,yl=ymin,yu=ymax), aes(xmin=xl, xmax=xu, ymin=yl, ymax=yu), fill="NA", color="black") +
  theme_linedraw() +
  theme(text=element_text(size=12)) +
  labs(x="log2(TPM) 24 hours into biofilm development", y="log2(TPM) 8 hours into biofilm development")
p
p %>% ggsave(filename="full_plot.pdf")

inset_p = plot_data %>%
  filter(logTPM_24h > xmin,
         logTPM_24h < xmax,
         logTPM_8h > ymin,
         logTPM_8h < ymax) %>%
  ggplot(aes(x=logTPM_24h, y=logTPM_8h, color=category)) +
  geom_vline(xintercept=perc_24h, linetype='dashed', color='red') +
  geom_hline(yintercept=perc_8h, linetype='dashed', color='red') +
  geom_point() +
  scale_color_discrete(name="Category", breaks=c("GR","LR","NAP","PNAP","Other"), labels=c("Global regulator", "Local regulator", "NAP", "Potential NAP", "Other")) +
  geom_text_repel(aes(label=gene), color='black', fontface='italic') +
  theme_linedraw() +
  theme(text=element_text(size=12)) +
  labs(x="log2(TPM) 24 hours into biofilm development", y="log2(TPM) 8 hours into biofilm development")
inset_p
inset_p %>% ggsave(filename="inset_plot.pdf")
